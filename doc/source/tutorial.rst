Tutorial
--------

The distribution contains a time series ``fes_tide_time_series.nc`` that will
be used in this help. This file is located in the ``tests/dataset`` directory
at the root of the project.

The first step is to read this time series using the NetCDF4 library, for
example:

.. code:: python

    import netCDF4
    import pytide

    with netCDF4.Dataset("tests/dataset/fes_tide_time_series.nc") as dataset:
        time = dataset['time'][:] * 1e-6    # microseconds to epoch
        h = dataset['ocean'][:] * 1e-2      # cm to m


Then, we will create an instance of a :py:class:`pytide.WaveTable` object:

.. code:: python

    wt = pytide.WaveTable()

By default, all components known by this object are loaded into memory. The
list of components known by this object can be retrieved using the
:py:meth:`pytide.WaveTable.known_constituents` method.

If you want to restrict the analysis to only a few components, you must provide
a list to the constructor in order to specify the waves to be analyzed.

.. code:: python

    wt = pytide.WaveTable(["M2", "K1", "O1", "P1", "Q1", "S1"])

The :meth:`pytide.WaveTable.constituents` method allows to retrieve the list of
waves defined during the construction of the object.

The different nodal corrections are then calculated from the time series to be
analyzed:

.. code:: python

    f, vu = wt.compute_nodal_corrections(time)
    f, vu = f.T, vu.T  # The matrices must be transposed.

These coefficients are used by :meth:`harmonic analysis
<pytide.WaveTable.harmonic_analysis>` to determine the properties of the
different tidal waves defined during the construction of the instance.

.. code:: python

    w = wt.harmonic_analysis(h, f, vu)

This result can then be used to determine a tidal height for the analyzed time
series:


.. code:: python

    hp = wt.tide(time, w)
