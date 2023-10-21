# Copyright (c) 2022 CNES
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Tidal constituents analysis
###########################
"""
from typing import TYPE_CHECKING, Dict, Iterator, List, Tuple, Union
import datetime
import distutils.version
import sys

import numpy

if TYPE_CHECKING and distutils.version.LooseVersion(
        numpy.__version__) >= distutils.version.LooseVersion(
            '1.20') and sys.version_info >= (3, 9):
    import numpy.typing
    NDArrayDateTime64 = numpy.typing.NDArray[numpy.datetime64]
    NDArrayComplex128 = numpy.typing.NDArray[numpy.complex128]
    NDArrayFloat64 = numpy.typing.NDArray[numpy.float64]
else:
    NDArrayDateTime64 = numpy.ndarray
    NDArrayComplex128 = numpy.ndarray
    NDArrayFloat64 = numpy.ndarray
from . import core, version

__version__ = version.release()
__date__ = version.date()


class AstronomicAngle(core.AstronomicAngle):
    """Initialize some astronomic data useful for nodal corrections.

    Args:
        date (datetime.datetime, optional): Desired UTC time
    """


class Wave(core.Wave):
    """Tidal wave properties."""


class WaveTable(core.WaveTable):
    """Properties of tidal constituents."""

    def __repr__(self) -> str:
        constituents = self.constituents()
        if len(constituents) > 9:
            constituents = constituents[:4] + ['...'] + constituents[-4:]

        return '{}.{}({})'.format(self.__class__.__module__,
                                  self.__class__.__name__,
                                  ', '.join(constituents))

    def freq(self) -> NDArrayFloat64:
        """Gets the waves frequencies in radians per seconds."""
        return numpy.array([wave.freq for wave in self], dtype=numpy.float64)

    def constituents(self) -> List[str]:
        """Gets the wave constituents handled by this instance."""
        return [wave.name() for wave in self]

    def compute_nodal_corrections(
            self, date: datetime.datetime) -> core.AstronomicAngle:
        """Compute nodal corrections.

        Args:
            date (datetime.datetime): Desired date

        Return:
            core.AstronomicAngle: The astronomic angle, indicating the date on
            which the tide is to be calculated.
        """
        return super().compute_nodal_corrections(date)

    def compute_nodal_modulations(
        self, dates: Union[List[datetime.datetime], NDArrayDateTime64]
    ) -> Tuple[NDArrayFloat64, NDArrayFloat64]:
        """Compute nodal modulations for amplitude and phase.

        Args:
            dates (list, numpy.ndarray): Desired dates

        Return:
            tuple: the nodal correction for amplitude, v (greenwich argument)
            + u (nodal correction for phase)
        """
        if isinstance(dates, list) and all(
                isinstance(item, datetime.datetime) for item in dates):
            epoch = numpy.array([core.timestamp(item) for item in dates]) * 1e6
            return super().compute_nodal_modulations(
                epoch.astype('datetime64[us]'))
        # The method throws an error if the dates are not datetime64
        return super().compute_nodal_modulations(dates)  # type: ignore

    @staticmethod
    def harmonic_analysis(h: NDArrayFloat64,
                          f: NDArrayFloat64,
                          vu: NDArrayFloat64,
                          dtype=None) -> NDArrayComplex128:
        """Harmonic Analysis.

        The harmonic analysis method consists in expressing the ocean tidal
        variations as a sum of independent constituents accordingly to the
        tidal potential spectrum. Then the sea surface elevation at a point
        :math:`(x, y)` and time :math:`t` can be expressed as a linear sum as
        follow:

        .. math::

            S_{ap} = S_{0}(x, y) + \\sum_{k=0}^n f_{k}(t)S_{k}(x, y)
            \\times cos [\\omega_{k}t + {v}_{k}(t) + u_{k}(t) - G_{k}(x,y)]

        where:

            * :math:`n` is the number of constituents,
            * :math:`S_{0}(x, y)` is the mean sea level,
            * :math:`S_{k}(x, y)` is the amplitude of the constituent of index
              :math:`k`,
            * :math:`G_{k}(x, y)` is the phase lag relative to Greenwich time,
            * :math:`w_{k}` is the angular frequency of the constituent of
              index :math:`k`,
            * :math:`v_{k}` is the astronomical argument at time :math:`t`,
            * :math:`f_{k}(t)` is the nodal correction coefficient applied to
              the amplitude of the constituent of index :math:`k`,
            * :math:`u_{k}(t)` is the nodal correction coefficient applied to
              the phase of the constituent of index :math:`k`.

        The a priori analysis spectrum includes the most important astronomical
        constituents in the Darwin development, completed by Shureman in 1958,
        and many non-linear waves. The definition of tidal constants and
        astronomical arguments is taken from FES2014 tidal prediction software
        and a complete definition of waves is also available in Shureman
        (1958). This spectrum is the most commonly used for harmonic analysis
        due the simplification given by the nodal correction concept (:math:`f`
        and :math:`u` coefficients above) which allows dealing with slow
        motions of the lunar ascending node and reducing the number of
        constituents in the tidal spectrum.

        More details about this harmonic analysis method can be found in
        Ponchaut et al. 1999.

        Args:
            h (numpy.ndarray): Sea level.
            f (numpy.ndarray): Nodal correction coefficient applied to the
                amplitude of the constituents analyzed.
            vu (numpy.ndarray): Astronomical argument at time :math:`t` + the
                nodal correction coefficient applied to the phase of the
                constituents analyzed
            dtype (numpy.dtype, optional): Data type of the complex numbers

        Returns:
            numpy.ndarray: The complex number representing the different
            reconstructed waves.
        """
        if dtype is None:
            return core.WaveTable.harmonic_analysis(h, f, vu)
        return core.WaveTable.harmonic_analysis(h, f, vu).astype(dtype)

    @staticmethod
    def select_waves_for_analysis(duration, n_periods=2) -> Iterator[str]:
        """Returns the list of tidal waves such that their period is more than
        twice the duration of the time series analyzed."""
        self = core.WaveTable()
        for wave in self:
            period = (numpy.pi * 2) / wave.freq / 86400
            if period * n_periods < duration:
                yield wave.name()


class WaveDict(WaveTable):
    """Manages the tidal wave table as a dictionary."""

    def freq(self):
        """Gets the waves frequencies in radians per seconds."""
        return {wave.name(): wave.freq for wave in self}

    def harmonic_analysis(self,
                          h: NDArrayFloat64,
                          f: NDArrayFloat64,
                          vu: NDArrayFloat64,
                          dtype=None) -> Dict[str, NDArrayComplex128]:
        """Harmonic Analysis.

        Args:
            h (numpy.ndarray): Sea level.
            f (numpy.ndarray): Nodal correction coefficient applied to the
                amplitude of the constituents analyzed.
            vu (numpy.ndarray): Astronomical argument at time :math:`t` + the
                nodal correction coefficient applied to the phase of the
                constituents analyzed
            dtype (numpy.dtype, optional): Data type of the complex number

        Returns:
            dict: A mapping between the wave name and its complex number
            representing it.

        .. seealso::

            :py:meth:`WaveTable.harmonic_analysis`
        """
        analysis = super().harmonic_analysis(h, f, vu, dtype=dtype)
        return {
            constituent: coefficient
            for constituent, coefficient in zip(self.constituents(), analysis)
        }

    def tide_from_tide_series(
            self, dates: NDArrayDateTime64,
            wave: Dict[str, numpy.complex128]) -> NDArrayFloat64:
        """Calculates the tide of a given time series.

        Args:
            dates (numpy.ndarray): time series data
            wave (dict): Tidal wave properties.

        Returns:
            numpy.ndarray: The tide calculated for the time series provided.
        """
        if len(wave) != len(self):
            raise ValueError('wave must contain as many items as tidal '
                             'constituents loaded')
        wave_properties = numpy.array([wave[item] for item in self])
        return super().tide_from_tide_series(dates, wave_properties)
