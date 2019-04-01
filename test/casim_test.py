""" :module casim_test: Test module for casim."""

# Import class to be tested.
from casim import casim

from collections import namedtuple
from test_utilities import _remove_test_files
from io import StringIO

import sys
import unittest

class casim_test(unittest.TestCase):
    """ :class: Test class for the casim """

    @classmethod
    def setUpClass(cls):
        """ Setup the test class. """

        # Setup a list of test files.
        cls._static_test_files = []

    @classmethod
    def tearDownClass(cls):
        """ Tear down the test class. """

        _remove_test_files(cls._static_test_files)

    def setUp (self):
        """ Setup the test instance. """

        # Setup list of test files to be removed immediately after each test method.
        self._test_files = []

    def tearDown (self):
        """ Tear down the test instance. """
        _remove_test_files(self._test_files)

    def test_10x10_seed_1(self):
        """ Run a test case with 100x100 cells and prng seed 1. """


        arguments = namedtuple('arguments', ('seed'))
        arguments.seed = 1

        # Capture stdout.
        old_stdout = sys.stdout
        stream = StringIO()
        sys.stdout = stream

        # Run the simulation.
        casim.main(arguments)

        # Read stdout.
        sim_out = stream.getvalue()

        # Reset stdout.
        sys.stdout = old_stdout

        # Remove last line (timing report).
        sim_out = sim_out.split("\n")

        self.assertIn("mut container updated [(0, 0), (0, 1), (1, 2), (1, 3), (3, 4), (3, 5), (2, 6)]", sim_out)

        self.assertIn(" [0 0 0 0 0 5 7 0 0 0]", sim_out)

if __name__ == "__main__":

    # Setup the test suite.
    suite = unittest.makeSuite(casim_test, 'test')

    # Run the suite.
    result = unittest.TextTestRunner(verbosity=2).run(suite)

    # Report test results.
    if result.wasSuccessful():
        print('---> All tests passed. <---')
        sys.exit(0)

    sys.exit(1)

