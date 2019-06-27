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

    def freq(self):
        """Gets the waves frequencies in radians per seconds"""
        return numpy.array([wave.freq for wave in self], dtype=numpy.float64)

    def constituents(self):
        """Gets the wave constituents handled by this instance"""
        return [wave.name() for wave in self]

    @staticmethod
    def harmonic_analysis(h, f, v0u, dtype=None):
        """Perform tidal harmonic analysis for a given signal h.

        f and v0u are the nodal factor calculated for the time series studied.

        This method returns the real and imaginary part of the
        constituents defined in nodal factors handled
        """
        if dtype is None:
            return core.WaveTable.harmonic_analysis(h, f, v0u)
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
