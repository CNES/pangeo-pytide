import argparse
import collections
import datetime
import logging
import math
import os
import pathlib
import re
import dask.array
import dask.distributed
import netCDF4
import numpy
import xarray
import pytide

# Pattern to decode the time contains in filenames
PATTERN = re.compile(r'\w+_t(\d+)\.nc').search

# Start date of the simulation
T0 = pytide.core.timestamp(datetime.datetime(2011, 9, 10))

# Size units
UNITS = [('octets', 0), ('ko', 0), ('Mo', 1), ('Go', 2), ('To', 2), ('Po', 2)]

# Definition of the calculation period of the analysis (the spin-up period is
# not included).
START_DATE = '2011-11-13'
END_DATE = '2012-11-12'


def t_axis(dirname):
    """Get the time axis representing the time series"""
    ds = xarray.open_zarr(dirname)
    return ds.dtime.values


def compute_nodal_corrections(client, waves, time_series):
    t = time_series.astype(numpy.float64) * 1e-9
    f, v0u = waves.compute_nodal_corrections(t)
    return (dask.array.from_delayed(client.scatter(f, broadcast=True),
                                    shape=f.shape,
                                    dtype=f.dtype),
            dask.array.from_delayed(client.scatter(v0u, broadcast=True),
                                    shape=v0u.shape,
                                    dtype=v0u.dtype))


def create_result(template,
                  target,
                  variable,
                  wave_table,
                  chunk=None,
                  dtype=numpy.float32):
    src = netCDF4.Dataset(template, mode="r")
    tgt = netCDF4.Dataset(target, mode="w")

    if chunk is None:
        chunk = None

    square_size = len(src.variables["j"][chunk])

    try:
        tgt.setncatts(src.__dict__)
        for dimname, dim in src.dimensions.items():
            if dimname == "time":
                continue
            size = len(dim)
            if dimname in ['j', 'i']:
                size = square_size
            tgt.createDimension(dimname, size)
        tgt.createDimension("wave", len(wave_table))

        var = tgt.createVariable("wave", str, ('wave', ))
        var.setncattr("long_name", "Tidal constituents")
        for idx, item in enumerate(wave_table):
            var[idx] = item.name()

        def variable_properties(varname):
            ncvar = src.variables[varname]
            if hasattr(ncvar, '_FillValue'):
                fill_value = ncvar._FillValue
            else:
                fill_value = None

            attributes = ncvar.__dict__
            if '_FillValue' in attributes:
                del attributes['_FillValue']

            return ncvar, fill_value, attributes

        # Copy variables.
        for varname in src.variables:
            if varname in [variable, 'iters', 'time']:
                continue

            ncvar, fill_value, attributes = variable_properties(varname)
            var = tgt.createVariable(varname,
                                     ncvar.dtype,
                                     ncvar.dimensions,
                                     fill_value=fill_value,
                                     zlib=True,
                                     complevel=6)
            var.setncatts(attributes)
            var[:] = ncvar[chunk] if varname in ['i', 'j'] else ncvar[:]

        # Create variables
        ncvar, _, _ = variable_properties(variable)
        var = tgt.createVariable(f"{variable}_real",
                                 dtype, ('face', 'wave', 'j', 'i'),
                                 fill_value=numpy.nan,
                                 zlib=True,
                                 complevel=9)
        var.setncatts(attributes)

        var = tgt.createVariable(f"{variable}_imag",
                                 dtype, ('face', 'wave', 'j', 'i'),
                                 fill_value=numpy.nan,
                                 zlib=True,
                                 complevel=9)
        var.setncatts(attributes)

    finally:
        tgt.close()
        src.close()


def write_one_face(path, waves, face, variable):
    if isinstance(waves, numpy.ma.MaskedArray):
        waves = waves.data
    with netCDF4.Dataset(path, mode="a") as ds:
        var = ds.variables[f'{variable}_real']
        var[face, :, :, :] = waves.real
        var = ds.variables[f'{variable}_imag']
        var[face, :, :, :] = waves.imag


def load_faces(dirname, face, variable, period, chunk):
    """Load a face from the time series"""
    ds = xarray.open_zarr(dirname)
    ds = ds.transpose("face", "j", "i", "time")
    return ds.isel(face=face, j=chunk, i=chunk, time=period)[variable].data


def pretty_size(size):
    """Display a size readable by a human"""
    if size in [0, 1]:
        return f"{size} octet"

    exponent = min(int(math.log(size, 1000)), len(UNITS) - 1)
    quotient = float(size) / 1000**exponent
    unit, num_decimals = UNITS[exponent]
    format_string = '{:.%sf} {}' % (num_decimals)
    return format_string.format(quotient, unit)


def dask_array_properties(da):
    """Display the dask array size handled"""
    chunks_size = 1
    for item in da.chunks:
        chunks_size *= item[0]
    chunks_size *= da.dtype.itemsize
    nblocks = numpy.array(da.numblocks).prod()
    return "array size=" + pretty_size(
        da.nbytes) + f" (#{nblocks} * " + pretty_size(
            chunks_size) + "); chunk size=" + str(da.chunksize)


def dask_array_rechunk(da, nblocks, axis=2):
    """TODO rechunk"""
    chunks = []
    div = int(math.sqrt(nblocks))
    for index, item in enumerate(da.chunks):
        chunks.append(numpy.array(item).sum() * (div if index == axis else 1))
    return tuple(item // div for index, item in enumerate(chunks))


def _apply_along_axis(arr, func1d, func1d_axis, func1d_args, func1d_kwargs):
    """Wrap apply_along_axis"""
    return numpy.apply_along_axis(func1d, func1d_axis, arr, *func1d_args,
                                  **func1d_kwargs)


def apply_along_axis(func1d, axis, arr, *args, **kwargs):
    """Apply the harmonic analysis to 1-D slices along the given axis."""
    arr = dask.array.core.asarray(arr)

    # Validate and normalize axis.
    arr.shape[axis]
    axis = len(arr.shape[:axis])

    # Rechunk so that analyze is applied over the full axis.
    arr = arr.rechunk(arr.chunks[:axis] + (arr.shape[axis:axis + 1], ) +
                      arr.chunks[axis + 1:])

    # Test out some data with the function.
    test_data = numpy.ones(args[0].shape[1], dtype=arr.dtype)
    test_result = numpy.array(func1d(test_data, *args, **kwargs))

    # Map analyze over the data to get the result
    # Adds other axes as needed.
    result = arr.map_blocks(
        _apply_along_axis,
        name=dask.utils.funcname(func1d) + '-along-axis',
        dtype=test_result.dtype,
        chunks=(arr.chunks[:axis] + test_result.shape + arr.chunks[axis + 1:]),
        drop_axis=axis,
        new_axis=list(range(axis, axis + test_result.ndim, 1)),
        func1d=func1d,
        func1d_axis=axis,
        func1d_args=args,
        func1d_kwargs=kwargs,
    )

    return result


def valid_directory(s):
    if not pathlib.Path(s).is_dir():
        raise argparse.ArgumentTypeError(
            "The directory does not exist or is not readable: " + s)
    return s


def valid_slice(s):
    if s == "None":
        return None
    try:
        return int(s)
    except ValueError:
        raise argparse.ArgumentTypeError("invalid slice value:" + s)


def valid_date(s):
    try:
        date = datetime.datetime.strptime("%Y-%m-%dT%H:%M:%S")
        return numpy.datetime64(date)
    except ValueError:
        raise argparse.ArgumentTypeError("invalid date value:" + s)


def usage():
    parser = argparse.ArgumentParser(
        description="Harmonic analysis for MIT/GCM")
    parser.add_argument("dirname",
                        help="Zarr directory containing MIT/GCM data",
                        type=valid_directory)
    parser.add_argument("template",
                        help="Template of the NetCDF file to be taken into "
                        "account",
                        type=argparse.FileType(mode='r'))
    parser.add_argument("variable",
                        help="Variable of the dataset to be processed")
    parser.add_argument("result", help="Path to the NetCDF File to create")
    parser.add_argument("--log",
                        metavar='PATH',
                        help="path to the logbook to use",
                        type=argparse.FileType("w"))
    parser.add_argument("--face",
                        help="face to process",
                        metavar='INT',
                        type=int,
                        nargs="+",
                        default=list(range(13)))
    parser.add_argument("--time_chunk",
                        help="Slice of the time axis to process",
                        metavar='INT',
                        nargs=3,
                        type=valid_slice,
                        default=None)
    parser.add_argument("--start_date",
                        help="Start date to process",
                        metavar='DATE',
                        nargs=1,
                        type=valid_date,
                        default=START_DATE)
    parser.add_argument("--end_date",
                        help="End date to process",
                        metavar='DATE',
                        nargs=1,
                        type=valid_date,
                        default=END_DATE)
    parser.add_argument("--var_chunk",
                        help="Slice of the variable to process",
                        metavar='INT',
                        nargs=3,
                        type=valid_slice,
                        default=None)
    parser.add_argument("--scheduler_file",
                        help="Path to a file with scheduler information",
                        metavar='PATH',
                        type=argparse.FileType("r"),
                        default=None)
    parser.add_argument("--local-cluster",
                        help="Use a dask local cluster for testing purpose",
                        action="store_true")
    consitutents = pytide.WaveTable.known_constituents()
    parser.add_argument("--tidal_constituents",
                        help="List of tidal waves to be studied. "
                        "Choose from the following consitutents: " +
                        ", ".join(consitutents),
                        metavar='WAVE',
                        nargs="+",
                        choices=consitutents,
                        default=[])
    parser.add_argument(
        "--nblocks",
        help="Number of grid divisions for the calculation of the harmonic"
        "analysis.",
        type=int,
        default=200)
    args = parser.parse_args()
    if not args.local_cluster:
        if args.scheduler_file is None:
            args.scheduler_file = str(pathlib.Path.home().joinpath(
                pathlib.Path("scheduler.json")))
        argparse.FileType('r')(args.scheduler_file)

    def setup_chunk(chunk):
        if chunk is None:
            return slice(0, None)
        return slice(*tuple(chunk))

    args.time_chunk = setup_chunk(args.time_chunk)
    args.var_chunk = setup_chunk(args.var_chunk)

    return args


def setup_logging(filename):
    """Setup the logging system"""
    kwargs = dict(format='%(asctime)s: %(message)s', level=logging.INFO)
    if filename is not None:
        filename.close()
        kwargs["filename"] = filename.name
    logging.basicConfig(**kwargs)


def main():
    args = usage()
    setup_logging(args.log)

    client = dask.distributed.Client(
        dask.distributed.LocalCluster(threads_per_worker=1)
    ) if args.local_cluster else dask.distributed.Client(
        scheduler_file=args.scheduler_file)
    logging.info(client)

    # Reading the list of files and associated dates to
    # be processed.
    time_series = t_axis(args.dirname)
    period = (time_series >= args.start_date) & (time_series <= args.end_date)
    logging.info("number of files to process %d", len(time_series[period]))
    logging.info("period [%s, %s]", time_series[period].min(),
                 time_series[period].max())

    wave_table = pytide.WaveTable(args.tidal_constituents)
    logging.info("%d tidal constituents to be analysed", len(wave_table))

    f, v0u = compute_nodal_corrections(client, wave_table, time_series[period])

    if not os.path.exists(args.result):
        # Create the result file
        logging.info("create the result file %r", args.result)
        create_result(args.template.name,
                      args.result,
                      args.variable,
                      wave_table,
                      chunk=args.var_chunk)
    else:
        # otherwise, we test the access to the file before continuing.
        with netCDF4.Dataset(args.result):
            pass

    del wave_table
    del time_series

    # Load the time series
    for face in args.face:
        logging.info("processing face %d", face)
        ds = load_faces(args.dirname,
                        face,
                        args.variable,
                        period,
                        chunk=args.var_chunk)
        logging.info("loaded %s", dask_array_properties(ds))
        ds = ds.rechunk(dask_array_rechunk(ds, args.nblocks))
        logging.info("fragmented %s", dask_array_properties(ds))
        future = apply_along_axis(pytide.WaveTable.harmonic_analysis, 2, ds,
                                  *(f, v0u))
        result = future.compute()
        result = numpy.transpose(result, [2, 0, 1])
        logging.info("write face #%d", face)
        write_one_face(args.result, result, face, args.variable)
        logging.info("calculation completed for face #%d", face)
        client.cancel(ds)

    logging.info("calculation done")


if __name__ == "__main__":
    main()
