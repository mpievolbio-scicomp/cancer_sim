""" :module casim_test: Test module for casim."""

# Import class to be tested.
from casim import casim
from casim.casim import CancerSimulator, CancerSimulatorParameters, check_set_number, load_cancer_simulation, LOGGER

from collections import namedtuple
from test_utilities import _remove_test_files
from io import StringIO
import numpy
import logging
from subprocess import Popen
import re
import os
import sys
import unittest
from tempfile import mkdtemp


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
        self.assertEqual(parameters.advantageous_death_probability,   0.0)
        self.assertEqual(parameters.mutation_rate,                         0.8)
        self.assertEqual(parameters.advantageous_mutation_probability,     1  )
        self.assertEqual(parameters.mutations_per_division,                1  )
        self.assertEqual(parameters.time_of_advantageous_mutation,     50000  )
        self.assertEqual(parameters.number_of_clonal,                     1   )
        self.assertEqual(parameters.tumour_multiplicity,              'single')

    def test_shaped_constructor (self):
        """ Test initialization with arguments. """

        parameters = CancerSimulatorParameters(
                                matrix_size =                           20  ,
                                number_of_generations =                 10  ,
                                division_probability =                   0.5,
                                advantageous_division_probability =      0.3,
                                death_probability =                      0.1,
                                advantageous_death_probability =    0.4,
                                mutation_rate =                          0.2,
                                advantageous_mutation_probability =      0.8,
                                mutations_per_division =                10  ,
                                time_of_advantageous_mutation =      30000  ,
                                number_of_clonal =                       2  ,
                                tumour_multiplicity =                'single',
                                                )

        self.assertEqual(parameters.matrix_size,                          20  )
        self.assertEqual(parameters.number_of_generations,                10  )
        self.assertEqual(parameters.division_probability,                  0.5)
        self.assertEqual(parameters.advantageous_division_probability,     0.3)
        self.assertEqual(parameters.death_probability,                     0.1)
        self.assertEqual(parameters.advantageous_death_probability,   0.4)
        self.assertEqual(parameters.mutation_rate,                         0.2)
        self.assertEqual(parameters.advantageous_mutation_probability,     0.8)
        self.assertEqual(parameters.mutations_per_division,               10  )
        self.assertEqual(parameters.time_of_advantageous_mutation,     30000  )
        self.assertEqual(parameters.number_of_clonal,                      2  )
        self.assertEqual(parameters.tumour_multiplicity,               'single' )

    def test_check_set_number(self):
        """ Test the numer checking utility. """

        self.assertEqual(1, check_set_number(1, int))
        self.assertEqual(1, check_set_number(1, int, None, 0, 10))
        self.assertEqual(1, check_set_number(1.0, int, None, 0, 10))

        self.assertRaises(TypeError, check_set_number, numpy.array([1.,2.]), float)


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

        # Test that the Simulator cannot be constructed without parameters.
        with self.assertRaises(ValueError):
            casim = CancerSimulator()

    def test_reference_run_nondefaults(self):
        """ Test running a simulation with non-default parameters and check values. """

        default_parameters = CancerSimulatorParameters(
                                matrix_size =                           20  ,
                                number_of_generations =                 10  ,
                                division_probability =                   0.5,
                                advantageous_division_probability =      0.8,
                                death_probability =                      0.1,
                                advantageous_death_probability =    0.4,
                                mutation_rate =                          0.2,
                                advantageous_mutation_probability =      0.8,
                                mutations_per_division =                10  ,
                                time_of_advantageous_mutation =      30000  ,
                                number_of_clonal =                       2  ,
                               )


        cancer_sim = CancerSimulator(default_parameters, seed=1, outdir='/tmp/1')

        cancer_sim.run()

        # Get cell matrix and compare to reference result.
        matrix = cancer_sim._CancerSimulator__mtx

        print(matrix.todense())
        reference_matrix = numpy.array(
                                       [[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                                        [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                                        [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 55, 54,  0,  0,  0,  0,  0,  0,  0,  0],
                                        [ 0,  0,  0,  0,  0,  0, 62,  1, 52,  1,  1, 55,  1,  1,  1,  0,  0,  0,  0,  0],
                                        [ 0,  0,  0,  0,  0,  0,  0, 63, 38, 53,  1, 36,  1,  1,  1,  0, 15,  0,  0,  0],
                                        [ 0,  0,  0,  0,  0,  0,  1, 46,  1, 39,  1,  1, 37,  1,  1, 15,  0,  0,  0,  0],
                                        [ 0,  0,  0,  0,  0, 12,  1,  1, 47,  1,  1,  1,  1,  1, 15, 22,  0,  0,  0,  0],
                                        [ 0,  0,  0,  0,  0, 28, 12, 12,  1, 35,  1,  1,  1, 15, 22, 15, 42,  0,  0,  0],
                                        [ 0,  0,  0,  0, 28, 28, 12, 10, 12,  4,  1,  1, 15, 23, 43, 42, 48,  0,  0,  0],
                                        [ 0,  0,  0,  0, 56, 12, 28, 12, 10, 13,  5,  1,  1, 25, 18, 49, 18,  0,  0,  0],
                                        [ 0,  0,  0,  0, 50, 57,  6, 29, 27, 11,  1,  3,  2,  1,  1, 18,  1,  1,  0,  0],
                                        [ 0,  0,  0,  0, 31, 51, 31, 26, 16, 17,  1,  1,  2,  2,  2,  1,  1,  1,  0,  0],
                                        [ 0,  0,  0,  0, 31, 31, 30, 26,  1,  1,  1,  1,  1,  2, 33, 41, 40,  0,  0,  0],
                                        [ 0,  0,  0,  0, 30, 30,  1,  1,  1,  1,  1,  9,  1,  1, 32, 61,  0,  0,  0,  0],
                                        [ 0,  0,  0,  0,  0, 45, 45,  1, 20, 21,  1,  8,  9,  1,  1,  0, 60,  0,  0,  0],
                                        [ 0,  0,  0,  0,  0, 44,  1, 20, 21,  1,  8,  1,  8,  1,  1,  0,  0,  0,  0,  0],
                                        [ 0,  0,  0,  0,  0,  1, 20, 20,  1,  1,  8,  8, 59,  8,  1,  1,  0,  0,  0,  0],
                                        [ 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  8,  0, 58,  1,  0,  0,  0,  0,  0,  0],
                                        [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                                        [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0]]
                                       )

        self.assertEqual( numpy.linalg.norm(matrix - reference_matrix), 0)

        # Get mutation container and compare to reference result.
        mutation_container = cancer_sim._CancerSimulator__mut_container
        reference_mutation_container = [(0, 0), (0, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (1, 10), (1, 11), (7, 12), (7, 13), (1, 14), (1, 15), (1, 16), (1, 17), (1, 18), (1, 19), (1, 20), (1, 21), (14, 22), (14, 23), (19, 24), (19, 25), (6, 26), (6, 27), (6, 28), (6, 29), (6, 30), (6, 31), (2, 32), (2, 33), (1, 34), (1, 35), (1, 36), (1, 37), (1, 38), (1, 39), (2, 40), (2, 41), (24, 42), (24, 43), (1, 44), (1, 45), (34, 46), (34, 47), (1, 48), (1, 49), (6, 50), (6, 51), (1, 52), (1, 53), (1, 54), (1, 55), (31, 56), (31, 57), (1, 58), (1, 59), (33, 60), (33, 61), (1, 62)]

        for r,m in zip(reference_mutation_container, mutation_container):
            self.assertEqual(r[0], m[0])
            self.assertEqual(r[1], m[1])

    def test_reference_run_defaults(self):
        """ Test running a simulation with default parameters and check values. """

        default_parameters = CancerSimulatorParameters()

        cancer_sim = CancerSimulator(default_parameters, seed=1)

        cancer_sim.run()

        # Get cell matrix and compare to reference result.
        matrix = cancer_sim._CancerSimulator__mtx

        reference_matrix = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0, 6, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 5, 7, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 4, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                                       ]
                                      )

        self.assertEqual( numpy.linalg.norm(matrix - reference_matrix), 0)

        # Get mutation container and compare to reference result.
        mutation_container = cancer_sim._CancerSimulator__mut_container
        reference_mutation_container = [(0, 0),
                                        (0, 1),
                                        (1, 2),
                                        (1, 3),
                                        (3, 4),
                                        (3, 5),
                                        (2, 6),
                                       ]
        for r,m in zip(reference_mutation_container, mutation_container):
            self.assertEqual(r[0], m[0])
            self.assertEqual(r[1], m[1])

    def test_setup_io(self):
        """ Test the IO handling. """

        default_parameters = CancerSimulatorParameters()
        cancer_sim = CancerSimulator(default_parameters, seed=1)

        self.assertIsNone(cancer_sim.outdir)
        self.assertIsNone(cancer_sim._CancerSimulator__seeddir)
        self.assertIsNone(cancer_sim._CancerSimulator__logdir)
        self.assertIsNone(cancer_sim._CancerSimulator__simdir)

        # Create an empty dir.
        tmpdir = mkdtemp()
        self._test_files.append(tmpdir)

        # This should work.
        cancer_sim.outdir = tmpdir

        # Check value.
        self.assertEqual(cancer_sim.outdir, tmpdir)

        # Get seed dir.
        seeddir = os.path.join(tmpdir, 'cancer_%d' % cancer_sim._CancerSimulator__seed)
        self.assertEqual(cancer_sim._CancerSimulator__seeddir, seeddir)
        self.assertTrue(os.path.isdir(cancer_sim._CancerSimulator__seeddir))

        # Check all subdirectories are correctly named and exist.

        self.assertEqual(cancer_sim._CancerSimulator__logdir,
                         os.path.join(seeddir, 'log'))
        self.assertTrue(os.path.isdir(cancer_sim._CancerSimulator__logdir))

        self.assertEqual(cancer_sim._CancerSimulator__simdir,
                         os.path.join(seeddir, 'simOutput'))
        self.assertTrue(os.path.isdir(cancer_sim._CancerSimulator__simdir))

        # Check export_tumour flag
        self.assertTrue(cancer_sim._CancerSimulator__export_tumour)

        # But not twice.
        with self.assertRaises(IOError) as exc:
            cancer_sim.outdir = tmpdir

    def test_serialize(self):
        """ The the serialization of the entire object. """
        parameters = CancerSimulatorParameters()
        cancer_sim = CancerSimulator(parameters, seed=1, outdir=mkdtemp())

        # dump before run.
        cancer_sim.dump()

        # Reload
        loaded_simulation = load_cancer_simulation(cancer_sim.dumpfile)

        self.assertIsInstance(loaded_simulation, CancerSimulator)

        # Check parameters.
        loaded_parameters = loaded_simulation.parameters

        self.assertEqual(loaded_parameters.number_of_generations,                parameters.number_of_generations)
        self.assertEqual(loaded_parameters.matrix_size,                          parameters.matrix_size)
        self.assertEqual(loaded_parameters.number_of_generations,                parameters.number_of_generations)
        self.assertEqual(loaded_parameters.division_probability,                 parameters.division_probability)
        self.assertEqual(loaded_parameters.advantageous_division_probability,    parameters.advantageous_division_probability)
        self.assertEqual(loaded_parameters.death_probability,                    parameters.death_probability)
        self.assertEqual(loaded_parameters.advantageous_death_probability,  parameters.advantageous_death_probability)
        self.assertEqual(loaded_parameters.mutation_rate,                        parameters.mutation_rate)
        self.assertEqual(loaded_parameters.advantageous_mutation_probability,    parameters.advantageous_mutation_probability)
        self.assertEqual(loaded_parameters.mutations_per_division,               parameters.mutations_per_division)
        self.assertEqual(loaded_parameters.time_of_advantageous_mutation,        parameters.time_of_advantageous_mutation)
        self.assertEqual(loaded_parameters.number_of_clonal,                     parameters.number_of_clonal)
        self.assertEqual(loaded_parameters.tumour_multiplicity,              parameters.tumour_multiplicity)

        # Check we can run.
        loaded_simulation.run()

        # dump again.
        loaded_simulation.dump()

        # Load again
        loaded_again_simulation = load_cancer_simulation(loaded_simulation.dumpfile)

        # Run once more.
        loaded_again_simulation.run()

    def test_export_tumour_matrix(self):
        """ Test exporting the tumour matrix. """

        parameters = CancerSimulatorParameters()
        cancer_sim = CancerSimulator(parameters, seed=1, outdir = mkdtemp())
        cancer_sim.run()

        # Check files where created.
        listing = os.listdir(cancer_sim._CancerSimulator__simdir)
        for f in ['mtx.p', 'mut_container.p', 'death_list.p', 'mtx_VAF.txt']:
            self.assertIn(f, listing)


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

    def test_cli(self):
        """ Test the command line interface. """
        # Setup command.
        python = "python"
        module = casim.__file__

        # Run with seed only.
        args = ['1']

        proc = Popen([python, module] + args)
        proc.wait()
        self.assertEqual(proc.returncode, 0)

        # Run with positional argument.
        outdir = 'cancer_sim_output'
        self._test_files.append(outdir)
        args += ['-o', outdir ]
        proc = Popen([python, module] + args)
        proc.wait()
        self.assertEqual(proc.returncode, 0)

        # run with positional argument (long version).
        args = ['2', '--outdir', outdir]
        proc = Popen([python, module] + args)
        proc.wait()
        self.assertEqual(proc.returncode, 0)

        # run with positional argument (long version).
        args = ['3', '--outdir', outdir, '-vv']
        proc = Popen([python, module] + args)
        proc.wait()
        self.assertEqual(proc.returncode, 0)

    def test_10x10_seed_1(self):
        """ Run a test case with 10x10 cells and prng seed 1. """

        arguments = namedtuple('arguments', ('seed', 'outdir', 'loglevel'))
        arguments.seed = 1
        arguments.outdir='cancer_sim_out'
        arguments.loglevel = 2
        self._test_files.append(arguments.outdir)

        # Capture stdout.
        stream = StringIO()
        log = LOGGER
        for handler in log.handlers:
           log.removeHandler(handler)
        myhandler = logging.StreamHandler(stream)
        myhandler.setLevel(logging.DEBUG)
        log.addHandler(myhandler)

        # Run the simulation.
        casim.main(arguments)

        # Flush log.
        myhandler.flush()

        # Read stdout.
        sim_out = stream.getvalue()

        # Reset stdout.
        log.removeHandler(myhandler)
        handler.close()

        mut_container_regex = re.compile(r"1 \[\(1, 4.0\), \(2, 2.0\), \(3, 2.0\), \(4, 1.0\), \(5, 1.0\), \(6, 1.0\), \(7, 1.0\)\]")
        # self.assertRegex(sim_out, mut_container_regex)


if __name__ == "__main__":

    unittest.main()
