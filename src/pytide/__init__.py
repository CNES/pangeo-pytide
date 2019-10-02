# Copyright (c) 2019 CNES
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Tidal constituents analysis
###########################
"""
import numpy
from . import core


class AstronomicAngle(core.AstronomicAngle):
    """Initialize some astronomic data useful for nodal corrections.

Args:
  epoch (float, optional): Desired UTC time
"""
    pass


class Wave(core.Wave):
    """Tidal wave properties"""
    pass


class WaveTable(core.WaveTable):
    """Properties of tidal constituents"""

    def __repr__(self):
        return "%s.%s(%s)" % (self.__class__.__module__,
                              self.__class__.__name__,
                              ', '.join(self.constituents()))

    def freq(self, as_dict=True):
        """Gets the waves frequencies in radians per seconds"""
        if as_dict:
            return {wave.name(): wave.freq for wave in self}
        return numpy.array([wave.freq for wave in self], dtype=numpy.float64)

    def constituents(self):
        """Gets the wave constituents handled by this instance"""
        return [wave.name() for wave in self]

    def compute_nodal_corrections(self, time):
        f, vu = super(WaveTable, self).compute_nodal_corrections(time)
        return f.T, vu.T

    #@staticmethod
    def harmonic_analysis(self, h, f=None, vu=None, dtype=None, as_dict=True):
        """Harmonic Analysis

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
        and a complete definition of waves is also available in Shureman (1958).
        This spectrum is the most commonly used for harmonic analysis due the
        simplification given by the nodal correction concept (:math:`f` and
        :math:`u` coefficients above) which allows dealing with slow motions of
        the lunar ascending node and reducing the number of constituents in the
        tidal spectrum.

        More details about this harmonic analysis method can be found in
        Ponchaut et al. 1999.

        Args:
            h (numpy.ndarray): Sea level.
            f (numpy.ndarray): Nodal correction coefficient applied to the
                amplitude of the constituents analyzed.
            vu (numpy.ndarray): Astronomical argument at time :math:`t` + the
                nodal correction coefficient applied to the phase of the
                constituents analyzed

        Returns:
            numpy.ndarray: The complex number representing the different
            reconstructed waves.
        """
        if as_dict:
            return {constituent: coefficient for constituent, coefficient in
                    zip(self.constituents(),
                        core.WaveTable.harmonic_analysis(h, f, vu))}
        if dtype is None:
            return core.WaveTable.harmonic_analysis(h, f, vu)
        return core.WaveTable.harmonic_analysis.astype(dtype)

    @staticmethod
    def select_waves_for_analysis(duration, n_periods=2):
        """Returns the list of tidal waves such that their period is more than
        twice the duration of the time series analyzed.
        """
        self = core.WaveTable()
        for wave in self:
            period = (numpy.pi * 2) / wave.freq / 86400
            if period * n_periods < duration:
                yield wave.name()

    def tide_from_time_series(self, time, w):
        if isinstance(w,dict):
            _w = list(w.values())
        else:
            _w = w
        return super(WaveTable, self).tide_from_tide_series(time, _w)
