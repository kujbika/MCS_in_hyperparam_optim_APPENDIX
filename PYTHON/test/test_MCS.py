"""Tests for `trfl` backtester service"""
"""you can run it in terminal: python -m unittest tests.Backtester -v"""

import unittest
import os
import numpy as np
from main.MCS import bootstrap_sample as bs


class TestMCSCalculations(unittest.TestCase):

    def test_bootstrap(self):
        """test when the posi changes sign"""
        data = np.random.rand(1,10)
        a = bs(data, 3, 1000)
        print("lol")


