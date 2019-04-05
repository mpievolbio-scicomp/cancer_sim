""" :module casim_test: Test module for casim."""

# Import class to be tested.
from casim import casim
from casim.casim import CancerSimulator, CancerSimulatorParameters, check_set_number

from collections import namedtuple
from test_utilities import _remove_test_files
from io import StringIO
import numpy
import sys
import unittest

class CancerSimulatorParametersTest(unittest.TestCase):
    """ :class: Test class for the CancerSimulator """

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

    def test_default_constructor (self):
        """ Test initialization without arguments. """

        parameters = CancerSimulatorParameters()

        self.assertEqual(parameters.matrix_size,                          10  )
        self.assertEqual(parameters.number_of_generations,                 2  )
        self.assertEqual(parameters.division_probability,                  1  )
        self.assertEqual(parameters.advantageous_division_probability,     1  )
        self.assertEqual(parameters.death_probability,                     0  )
        self.assertEqual(parameters.fitness_advantage_death_probability,   0.0)
        self.assertEqual(parameters.mutation_rate,                         0.8)
        self.assertEqual(parameters.advantageous_mutation_probability,     1  )
        self.assertEqual(parameters.mutations_per_division,                1  )
        self.assertEqual(parameters.time_of_advantageous_mutation,     50000  )
        self.assertEqual(parameters.number_of_clonal,                     1   )

    def test_shaped_constructor (self):
        """ Test initialization with arguments. """

        parameters = CancerSimulatorParameters(
                                matrix_size =                           20  ,
                                number_of_generations =                 10  ,
                                division_probability =                   0.5,
                                advantageous_division_probability =      0.3,
                                death_probability =                      0.1,
                                fitness_advantage_death_probability =    0.4,
                                mutation_rate =                          0.2,
                                advantageous_mutation_probability =      0.8,
                                mutations_per_division =                10  ,
                                time_of_advantageous_mutation =      30000  ,
                                number_of_clonal =                       2  ,
                                                )

        self.assertEqual(parameters.matrix_size,                          20  )
        self.assertEqual(parameters.number_of_generations,                10  )
        self.assertEqual(parameters.division_probability,                  0.5)
        self.assertEqual(parameters.advantageous_division_probability,     0.3)
        self.assertEqual(parameters.death_probability,                     0.1)
        self.assertEqual(parameters.fitness_advantage_death_probability,   0.4)
        self.assertEqual(parameters.mutation_rate,                         0.2)
        self.assertEqual(parameters.advantageous_mutation_probability,     0.8)
        self.assertEqual(parameters.mutations_per_division,               10  )
        self.assertEqual(parameters.time_of_advantageous_mutation,     30000  )
        self.assertEqual(parameters.number_of_clonal,                      2  )


    def test_check_set_number(self):
        """ Test the numer checking utility. """

        self.assertEqual(1, check_set_number(1, int))
        self.assertEqual(1, check_set_number(1, int, None, 0, 10))
        self.assertEqual(1, check_set_number(1.0, int, None, 0, 10))

        self.assertRaises(TypeError, check_set_number, numpy.array([1.,2.]), float)

    def test_run(self):
        """ Test running a simulation. """

        default_parameters = CancerSimulatorParameters()

        cancer_sim = CancerSimulator(default_parameters)

        cancer_sim.run()

class CancerSimulatorTest(unittest.TestCase):
    """ :class: Test class for the CancerSimulator """

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

    def test_default_constructor (self):
        """ Test the construction of the Simulator without arguments. """

        casim = CancerSimulator()

        self.assertIsInstance(casim, CancerSimulator)


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

    unittest.main()

    ## Setup the test suite.
    #suite = unittest.makeSuite(casim_test, 'test')

    ## Run the suite.
    #result = unittest.TextTestRunner(verbosity=2).run(suite)

    ## Report test results.
    #if result.wasSuccessful():
        #print('---> All tests passed. <---')
        #sys.exit(0)

    sys.exit(1)

