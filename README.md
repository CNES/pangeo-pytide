[![conda](https://anaconda.org/conda-forge/pytide/badges/installer/conda.svg?service=github)](https://www.anaconda.com/distribution/)
[![platforms](https://anaconda.org/conda-forge/pytide/badges/platforms.svg?service=github)](https://anaconda.org/conda-forge/pytide)
[![latest-release-date](https://anaconda.org/conda-forge/pytide/badges/latest_release_date.svg?service=github)](https://github.com/CNES/pangeo-pytide/commits/master)
[![license](https://anaconda.org/conda-forge/pytide/badges/license.svg?service=github)](https://opensource.org/licenses/BSD-3-Clause)


# pangeo-pytide

## About

`pytide` allows to analyze the tidal constituents of a time series from a
[harmonic
analysis](https://pangeo-pytide.readthedocs.io/en/latest/pytide.html#pytide.WaveTable.harmonic_analysis).
The definition of tidal constants and astronomical arguments is taken from
[FES2014 tidal prediction
software](https://bitbucket.org/cnes_aviso/fes/src/master/).

It was developed to analyze the *MIT/GCM LLC4320* model. The script
["mit_gcm_detiding.py"](https://github.com/CNES/pangeo-pytide/blob/master/src/scripts/mit_gcm_detiding.py)
used to perform this analysis is distributed with this distribution.

## How To Install

[Anaconda](https://anaconda.org) is a free and open-source distribution of the
Python programming language for scientific computing, that aims to simplify
package management and deployment. Package versions are managed by the package
management system conda.

The first step is to install the anaconda distribution. The installation manual
for this software is detailed
[here](https://docs.anaconda.com/anaconda/install/).

To install the software using conda simply execute the following command:

    conda install pangeo-pytide -c conda-forge

This command will install the software and the necessary dependencies. More
information is available on the syntax of this command on the [related
documentation](https://conda.io/projects/conda/en/latest/commands/install.html).

> If you want to build the package yourself, you can find more information on
> the [help page](https://pangeo-pytide.readthedocs.io/en/latest/setup.html) of
> the project.

## Quick Tutorial

The distribution contains a time series
[fes_tide_time_series.nc](https://github.com/CNES/pangeo-pytide/blob/master/tests/dataset/fes_tide_time_series.nc)
that will be used in this help.

The first step is to read this time series using the NetCDF4 library, for
example:

    import netCDF4
    import pytide

    with netCDF4.Dataset("tests/dataset/fes_tide_time_series.nc") as dataset:
        time = dataset['time'][:] * 1e-6    # microseconds to epoch
        h = dataset['ocean'][:] * 1e-2      # cm to m

Then, we will create an instance of a `pytide.WaveTable` object:

    wt = pytide.WaveTable()

By default, all components known by this object are loaded into memory. The
list of components known by this object can be retrieved using the
[pytide.WaveTable.known_constituents](https://pangeo-pytide.readthedocs.io/en/latest/pytide.html#pytide.WaveTable.known_constituents)
method.

If you want to restrict the analysis to only a few components, you must provide
a list to the constructor in order to specify the waves to be analyzed.

    wt = pytide.WaveTable(["M2", "K1", "O1", "P1", "Q1", "S1"])

The
[pytide.WaveTable.constituents](https://pangeo-pytide.readthedocs.io/en/latest/pytide.html#pytide.WaveTable.constituents)
method allows to retrieve the list of waves defined during the construction of
the object.

The different nodal corrections are then calculated from the time series to be
analyzed:

    f, vu = wt.compute_nodal_corrections(time)
    f, vu = f.T, vu.T  # The matrices must be transposed.

These coefficients are used by [harmonic
analysis](https://pangeo-pytide.readthedocs.io/en/latest/pytide.html#pytide.WaveTable.harmonic_analysis)
to determine the properties of the different tidal waves defined during the
construction of the instance.

    w = wt.harmonic_analysis(h, f, vu)

This result can then be used to determine a tidal height for the analyzed time
series:

    hp = wt.tide(time, w)
