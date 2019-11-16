# Copyright (c) 2019 CNES
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
import datetime
import os
import math
import numpy
import unittest
import netCDF4
from pytide import core


class AstronomicAngle(unittest.TestCase):
    def test_init(self):
        aa = core.AstronomicAngle(core.timestamp(datetime.datetime(2000, 1,
                                                                   1)))
        self.assertTrue(isinstance(aa, core.AstronomicAngle))
        self.assertAlmostEqual(aa.h, 4.886452089967941, delta=1e-6)
        self.assertAlmostEqual(aa.n, 2.182860931126595, delta=1e-6)
        self.assertAlmostEqual(aa.nu, 0.20722218671046477, delta=1e-6)
        self.assertAlmostEqual(aa.nuprim, 0.13806065629468897, delta=1e-6)
        self.assertAlmostEqual(aa.nusec, 0.13226438100551682, delta=1e-6)
        self.assertAlmostEqual(aa.p, 1.4537576754171415, delta=1e-6)
        self.assertAlmostEqual(aa.p1, 4.938242223271549, delta=1e-6)
        self.assertAlmostEqual(aa.r, 0.1010709894525481, delta=1e-6)
        self.assertAlmostEqual(aa.s, 3.6956255851976114, delta=1e-6)
        self.assertAlmostEqual(aa.t, 3.1415926536073755, delta=1e-6)
        self.assertAlmostEqual(aa.x1ra, 1.1723206438502318, delta=1e-6)
        self.assertAlmostEqual(aa.xi, 0.1920359426758722, delta=1e-6)


class Wave(unittest.TestCase):
    def test_init(self):
        with self.assertRaises(TypeError):
            core.Wave()


class WaveTable(unittest.TestCase):
    DATASET = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "dataset", "fes_tide_time_series.nc")

    def test_init(self):
        wt = core.WaveTable()
        self.assertEqual(len(wt), 67)
        self.assertEqual(len([item for item in wt]), 67)
        self.assertEqual(wt.wave("M2"), wt.wave(core.Wave.Ident.kM2))
        self.assertNotEqual(wt.wave("M2"), wt.wave(core.Wave.Ident.kK1))
        self.assertTrue(wt.wave("__M2__") is None)
        self.assertListEqual(sorted(wt.known_constituents()),
                             sorted([item.name() for item in wt]))
        for item in wt:
            self.assertEqual(item.ident,
                             getattr(core.Wave.Ident, "k" + item.name()))

        wt = core.WaveTable(["M2", "K1", "O1", "P1", "Q1", "S1"])
        self.assertEqual(len(wt), 6)
        self.assertListEqual(sorted([item.name() for item in wt]),
                             sorted(["M2", "K1", "O1", "P1", "Q1", "S1"]))

    def test_wave(self):
        aa = core.AstronomicAngle(core.timestamp(datetime.datetime(2000, 1,
                                                                   1)))
        wt = core.WaveTable(["M2"])
        wave = wt.wave("M2")
        self.assertAlmostEqual(wave.freq * 86400,
                               12.140833182614747,
                               delta=1e-6)
        self.assertEqual(wave.type, wave.TidalType.kShortPeriod)

    def test_analysis(self):
        with netCDF4.Dataset(self.DATASET) as dataset:
            time = dataset['time'][:] * 1e-6
            h = dataset['ocean'][:] * 1e-2

        wt = core.WaveTable()
        f, vu = wt.compute_nodal_corrections(time)
        w = wt.harmonic_analysis(h, f, vu)
        delta = h - wt.tide_from_tide_series(time, w)

        self.assertAlmostEqual(delta.mean(), 0, delta=1e-16)
        self.assertAlmostEqual(delta.std(), 0, delta=1e-12)


if __name__ == "__main__":
    unittest.main()
