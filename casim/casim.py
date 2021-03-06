# -*- coding: utf-8 -*-
#!/usr/bin/env python3
__author__ = 'Luka Opasic, MD'
__email__ = 'opasic@evolbio.mpg.de'
__version__ = '1.1.0'


from argparse import ArgumentParser
from operator import itemgetter
from random import shuffle
from scipy.sparse import lil_matrix
from time import sleep, time
from timeit import default_timer as timer
import dill
import itertools
import logging
import math
import matplotlib.pyplot as plt
import numpy
import os
import pickle
import random as prng
import sys
from importlib.util import spec_from_file_location, module_from_spec

np = numpy

LOG_FORMAT_STR = '%(asctime)s %(levelname)s: %(message)s'
logging.basicConfig(
    level=logging.INFO,
    format=LOG_FORMAT_STR,
    handlers=[
        logging.StreamHandler(sys.stderr)
    ]
)

class CancerSimulatorParameters(object):
    """
        :class CancerSimulatorParameters: Represents the parameters for a cancer simulation.
    """

    def __init__(self,
                 matrix_size=None,
                 number_of_generations=None,
                 division_probability=None,
                 adv_mutant_division_probability=None,
                 death_probability=None,
                 adv_mutant_death_probability=None,
                 mutation_probability=None,
                 adv_mutant_mutation_probability=None,
                 number_of_mutations_per_division=None,
                 adv_mutation_wait_time=None,
                 number_of_initial_mutations=None,
                 tumour_multiplicity=None,
                 read_depth=None,
                 sampling_fraction=None,
                 plot_tumour_growth=None,
                 export_tumour=None,
                 sampling_positions=None,
                ):
        """
        Construct a new CancerSimulationParameters object.

        :param matrix_size: The size of the (square) grid in each dimension. 
        :type  matrix_size: int (matrix_size > 0)

        :param number_of_generations: The number of generations to simulate.
        :type  number_of_generations: int (number_of_generations > 0)

        :param division_probability: The probability for a cell division to occur during one generation.
        :type  division_probability: float (0.0 <= division_probability <= 1.0)

        :param adv_mutant_division_probability: The probability for the division of a cell with advantageous mutation to occur during one generation.
        :type  adv_mutant_division_probability: float (0.0 <= adv_mutant_division_probability <= 1.0)

        :param death_probability: The probability for a cell to die during one generation.
        :type  death_probability: float (0.0 <= death_probability <= 1.0)

        :param adv_mutant_death_probability: The probability for a cell with advantageous mutation to die during one generation.
        :type  adv_mutant_death_probability: float (0.0 <= adv_mutant_death_probability <= 1.0)

        :param mutation_probability: The probability of mutation.
        :type  mutation_probability: float (0.0 <= mutation_probability <= 1.0)

        :param adv_mutant_mutation_probability: The probability for a mutation to occur during one generation in a cell with adv. mutation.
        :type  adv_mutant_mutation_probability: float (0.0 <= adv_mutant_mutation_probability <= 1.0)

        :param number_of_mutations_per_division: The number of mutations per division
        :type  number_of_mutations_per_division: int (0 < number_of_mutations_per_division)

        :param adv_mutation_wait_time: The number of generations into the simulation after which the advantageous mutation is inserted.
        :type  adv_mutation_wait_time: int (adv_mutation_wait_time > 0)

        :param number_of_initial_mutations: Number of mutations present in first cancer cell.
        :type  number_of_initial_mutations: int (number_of_initial_mutations >= 0)

        :param tumour_multiplicity: Run in single or double tumour mode (i.e. consider growth of one single tumour or two tumours simultaneously). Possible values: "single", "double".
        :type  tumour_multiplicity: str

        :param read_depth: The sequencing read depth (read length * number of reads / genome length). Default: 100.
        :type  read_depth: int (read_depth >= 0)

        :param sampling_fraction: The fraction of cells to include in a sample. Default: 0.
        :type  sampling_fraction: float  (0 <= sampling_fraction <= 1)

        :param sampling_positions: The positions of cells to include in a sample. Default: Random position.
        :type  sampling_positions: List (or array) of tuples of ints. E.g.  ([10,20], [2,31]).

        :param plot_tumour_growth: Render graph of the tumour size as function of time. Default: True.
        :type plot_tumour_growth: bool

        :param export_tumour: Dump the tumour data to file. Default: True.
        :type export_tumour: bool

        """
        # Store parameters on the object.
        self.matrix_size = matrix_size
        self.number_of_generations = number_of_generations
        self.division_probability = division_probability
        self.adv_mutant_division_probability = adv_mutant_division_probability
        self.death_probability = death_probability
        self.adv_mutant_death_probability = adv_mutant_death_probability
        self.mutation_probability = mutation_probability
        self.adv_mutant_mutation_probability = adv_mutant_mutation_probability
        self.number_of_mutations_per_division = number_of_mutations_per_division
        self.adv_mutation_wait_time = adv_mutation_wait_time
        self.number_of_initial_mutations = number_of_initial_mutations
        self.tumour_multiplicity = tumour_multiplicity
        self.read_depth = read_depth
        self.sampling_fraction = sampling_fraction
        self.sampling_positions = sampling_positions
        self.plot_tumour_growth = plot_tumour_growth
        self.export_tumour = export_tumour
    
        
    @property
    def matrix_size(self):
        return self.__matrix_size
    @matrix_size.setter
    def matrix_size(self, val):
        self.__matrix_size = check_set_number(val, int, 10, 1, None )

    @property
    def number_of_generations(self):
        return self.__number_of_generations
    @number_of_generations.setter
    def number_of_generations(self, val):
        self.__number_of_generations = check_set_number(val, int, 2, 1, None)

    @property
    def division_probability(self):
        return self.__division_probability
    @division_probability.setter
    def division_probability(self, val):
        self.__division_probability = check_set_number(val, float, 1, 0.0, 1.0)

    @property
    def adv_mutant_division_probability(self):
        return self.__adv_mutant_division_probability
    @adv_mutant_division_probability.setter
    def adv_mutant_division_probability(self, val):
        self.__adv_mutant_division_probability = check_set_number(val, float, 1, 0.0, 1.0)

    @property
    def death_probability(self):
        return self.__death_probability
    @death_probability.setter
    def death_probability(self, val):
        self.__death_probability = check_set_number(val, float, 0, 0.0, 1.0)

    @property
    def adv_mutant_death_probability(self):
        return self.__adv_mutant_death_probability
    @adv_mutant_death_probability.setter
    def adv_mutant_death_probability(self, val):
        self.__adv_mutant_death_probability = check_set_number(val, float, 0.0, 0.0, 1.0)

    @property
    def mutation_probability(self):
        return self.__mutation_probability
    @mutation_probability.setter
    def mutation_probability(self, val):
        self.__mutation_probability = check_set_number(val, float, 0.8, 0.0, 1.0)

    @property
    def adv_mutant_mutation_probability(self):
        return self.__adv_mutant_mutation_probability
    @adv_mutant_mutation_probability.setter
    def adv_mutant_mutation_probability(self, val):
        self.__adv_mutant_mutation_probability = check_set_number(val, float, 1.0, 0.0, 1.0)

    @property
    def number_of_mutations_per_division(self):
        return self.__number_of_mutations_per_division
    @number_of_mutations_per_division.setter
    def number_of_mutations_per_division(self, val):
        self.__number_of_mutations_per_division = check_set_number(val, int, 1, 0)

    @property
    def adv_mutation_wait_time(self):
        return self.__adv_mutation_wait_time
    @adv_mutation_wait_time.setter
    def adv_mutation_wait_time(self, val):
        self.__adv_mutation_wait_time = check_set_number(val, int, 50000, 0)

    @property
    def number_of_initial_mutations(self):
        return self.__number_of_initial_mutations
    @number_of_initial_mutations.setter
    def number_of_initial_mutations(self, val):
        self.__number_of_initial_mutations = check_set_number(val, int, 1, 0)

    @property
    def tumour_multiplicity(self):
        return self.__tumour_multiplicity
    @tumour_multiplicity.setter
    def tumour_multiplicity(self,val):
        if val is None:
            val = 'single'
        if not isinstance(val, str):

            raise TypeError("Wrong type for parameter 'tumour_multiplicity'. Expected str, got %s" % type(val))

        if val not in ["single", "double"]:
            raise ValueError("Only 'single' and 'double' are allowed values for parameter 'tumour_multiplicity'.")

        self.__tumour_multiplicity = val

    @property
    def read_depth(self):
        return self.__read_depth
    @read_depth.setter
    def read_depth(self, val):
        self.__read_depth = check_set_number(val, int, 100, 0, None)

    @property
    def sampling_fraction(self):
        return self.__sampling_fraction
    @sampling_fraction.setter
    def sampling_fraction(self, val):
        self.__sampling_fraction = check_set_number(val, float, 0.0, 0.0, 1.0)

    @property
    def sampling_positions(self):
        return self.__sampling_positions
    @sampling_positions.setter
    def sampling_positions(self, val):
        if val is not None:
            for pos in val:
                if not hasattr(val, "__iter__"):
                    raise TypeError("Sampling positions must be list of tuples.")
                if len(pos) != 2:
                    raise ValueError("Sampling positions must be list of 2-tuples (x,y coordinates).")
                for xy in pos:
                    if not isinstance(xy, int):
                        raise TypeError("Sampling position must be integer")
                    if xy < 0 or xy > self.matrix_size:
                        raise ValueError("Sampling position must be positive integer not larger than the matrix size.")

        self.__sampling_positions = val
            

    @property
    def plot_tumour_growth(self):
        return self.__plot_tumour_growth
    @plot_tumour_growth.setter
    def plot_tumour_growth(self, val):
        if val is None:
            val = True
        try:
            val = bool(val)
        except:
            raise TypeError("Incompatible type: Expected bool, got {}.".format(type(val)))
        
        self.__plot_tumour_growth = val

    @property
    def export_tumour(self):
        return self.__export_tumour
    @export_tumour.setter
    def export_tumour(self, val):
        if val is None:
            val = True
        try:
            val = bool(val)
        except:
            raise TypeError("Incompatible type: Expected bool, got {}.".format(type(val)))
        
        self.__export_tumour = val



class CancerSimulator(object):
    """
        :class CancerSimulator: Represents the Monte-Carlo simulation of cancer tumour growth on a 2D grid.
    """

    def __init__(self, parameters=None,
                       seed=None,
                       outdir=None,
                       ):
        """
        Construct a new CancerSimulation.

        :param parameters: The cancer simulation parameters
        :type  parameters: CancerSimulationParameters

        :param seed: The random seed.
        :type  seed: int

        :param outdir: The directory where simulation data is saved. Default: "casim_out/" in the current working directory.
        :type  outdir: (str || path-like object)
        """

        if parameters is None:
            raise ValueError("No parameters given, simulation cannot execute.")
        self.parameters = parameters

        # Setup internal variables.
        self.__mtx = lil_matrix((self.parameters.matrix_size, self.parameters.matrix_size), dtype=int)
        self.__mut_container = None
        self.__xaxis_histogram = None
        self.__biopsy_timing = None
        self.__beneficial_mutation = []
        self.__growth_plot_data = None
        self.__s = self.parameters.number_of_mutations_per_division


        # Keep track of how many steps where performed in previous run if this
        # is a reloaded run.
        self.__init_step = 0
        # Keep track of mutation count in previous run if this is a rerun.
        self.__mutation_counter = 1

        # Handle direct parameters.
        self.seed = seed
        self.outdir = outdir
        self.__ploidy=2
        self.__mut_multiplier=[self.__s]*100000
        
        if self.parameters.tumour_multiplicity == 'single':
            logging.info('Running in single tumour mode.')
            initLoc=(int(self.parameters.matrix_size*0.5),int(self.parameters.matrix_size*0.5))

            logging.info("First cell at %s.", str(initLoc))

            self.__mtx[initLoc]=1
            self.__mut_container=[(0, 0), (0, 1)]
            self.__pool=[initLoc]

        #start the pool of cancer cells by adding the initial cancer cell into it
        if self.parameters.tumour_multiplicity == 'double':
            logging.info('Running in sdsa mode.')

            ### COMMENT: Should these be given as parameters?
            distance_between_tumours=0.05
            initLoc=(int(self.parameters.matrix_size*0.45),int(self.parameters.matrix_size*0.5))
            secondinitLoc=(int(self.parameters.matrix_size*0.65),int(self.parameters.matrix_size*0.51))

            self.__mtx[initLoc]=1
            self.__mtx[secondinitLoc]=2
            self.__mut_container=[(0, 0), (0, 1), (0,2)]
            self.__pool=[initLoc, secondinitLoc]

        # create lists used in loops
        self.__growth_plot_data=[]


    @property
    def seed(self):
        return self.__seed
    @seed.setter
    def seed(self, val):
        """ Set the random seed for the simulation.
        :param val: The seed to set
        :type  val: int
        """

        # If not given: Set it to number of seconds since Jan. 1. 1970 (T0)
        if val is None:
            val = int(time())
        if not isinstance(val, int):
            raise TypeError("Wrong type for parameter 'seed'. Expected int, got %s." % type(val))
        if not val > 0:
            raise ValueError("The parameter 'seed' must a positive integer (int).")

        self.__seed = val

    @property
    def dumpfile(self):
        return self.__dumpfile

    @property
    def outdir(self):
        return self.__outdir
    @outdir.setter
    def outdir(self, val):
        """ Create the output directory if not existing. If simulation data is already present inside an existing directory, the simulation aborts. """
        self._setup_io(val)

    def _setup_io(self, outdir):
        """ """
        """ Setup the output directories.

        :param outdir: The directory under which all simulation output will be stored.
        :type  outdir: str
        :raises: IOError (Directory for this seed already exists)"""

        self.__outdir = './outdir'
        self.__seeddir = './seeddir'
        self.__logdir  = './logdir'
        self.__simdir = './simdir'

        if outdir is None:
            outdir = "casim_out"

        # Create top-level outdir.
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        self.__outdir = outdir

        # Check if directory for this seed already exists. Bail out if yes.
        seeddir = os.path.join(outdir, 'cancer_%d' % self.__seed)
        if os.path.exists(seeddir):
            raise IOError("The directory %s already exists. Cowardly refusing to overwrite. Please specify another seed or a different outdir." % seeddir)

        os.mkdir(seeddir)
        self.__seeddir = seeddir

        # Setup dump file
        self.__dumpfile = os.path.join(self.__seeddir, 'cancer_sim.py.dill')

        # Create subdirectories.
        logdir = os.path.join(seeddir, "log")
        simdir = os.path.join(seeddir, "simOutput")
        os.mkdir(logdir)
        os.mkdir(simdir)

        # Store on object.
        self.__logdir = logdir
        self.__simdir = simdir

        # Configure the logging filehandler.
        root_logger = logging.getLogger()
        fhandler = logging.FileHandler(os.path.join(logdir, "casim.log"))
        fhandler.setFormatter(root_logger.handlers[0].formatter)
        root_logger.addHandler(fhandler)

    def extend_sample(self, sample_center, sample_size):
        """ Takes a subset of cells from the tumour positioned around single input cell with specific coordinates. Output is a list of tuples of cells belonging to the sample.
        :param sample_center: coordinates of cell that will be center of the sample
        :type  sample: tuple

        :param sample_size: The size of the sample (fraction of total cells.)
        :type  sample_size: float
        """

        biopsy_size=math.ceil(sample_size*len(self.__pool))

        #look at z-tier neighbours around sample_center
        for z in range(1,len(self.__pool)):
            expanded_sample=[]

            #look at z-tier neighbours around sample_center
            for i in range(-z, z+1):
                for j in range(-z, z+1):
                    nc=(sample_center[0]+i, sample_center[1]+j)

                    #if surrounding cell in the pool add it to sample list
                    if nc in self.__pool:
                        expanded_sample.append(nc)

            # if chunk is larger than wanted percentage of total tumour
            if len(expanded_sample)>biopsy_size:

                #remove last value until desired chunk size
                while len(expanded_sample)>biopsy_size:
                    expanded_sample=expanded_sample[:-1]
                break

        if len(expanded_sample) == 0:
            logging.warning("""
Sample is empty. Consider enlarging the `sampling_fraction` parameter.
If that does not help, you may be sampling an empty region of the tumour matrix.
Inspect the tumour matrix data `mtx.p` in the output directory""")

        return expanded_sample


    def dump(self):
        """ Serialize the object. The current simulation will be stored in a machine readable
        format to <OUTDIR>/cancer_<SEED>/cancer_sim.py.dill, where <OUTDIR> is the specified output
        directory or (if the latter was not defined) a temporary directory."""

        with open(self.dumpfile, 'wb') as fp:
            dill.dump(self, fp)

    def run(self):
        """ Run the simulation.
        
        :return: 0 if the run finishes successfully.

        After a successful run, simulation output and log will be written to
        the output directory `<DIR>/cancer_<SEED>/simOutput` and
        `<DIR>/cancer_<SEED>/log`, respectively. Simulation output is split into
        several files:

        - `mtx_VAF.txt` is a datafile with three columns: `mutation_id` lists the index of
          each primary mutation, `additional_mut_id` indexes the subsequent mutations that occur in a cell of
          a given `mutation_id`; `frequency` is the frequency at which a given mutation occurs. 

        - `sample_out_XXX_YYY.txt` lists all mutations of the artificial sample
          taken from the whole tumour. Columns are identical to `mtx_VAF.txt`.
         
        - `wholeTumourVAFHistogram.pdf` contains a histogram plot of the
          mutation frequencies for the  whole tumour
        - `sampleHistogram_XXX_YYY.pdf` is the mutation frequency histogram for
          the sampled portion of the tumour. The two numbers XXX and YYY are the
          positional coordinates (grid indices) in the tumour matrix.

        - `mtx.p` is the serialized (aka "pickled") 2D tumour matrix in sparse
          matrix format.
        - `mut_container.p` is the serialized (aka "pickled") mutation list, a
          list of tuples `[t_i]`. Each tuple `t_i` consists of two values, `t_i =
          (c_i, m_i)`. The first element `c_i` is the cell number in which the i'th mutation
          occurs. The second element, `m_i`, is the mutation index `m_i=i`. 
        """ 
        # Setup square matrix.
        matrix_size=self.parameters.matrix_size

        self.__pre_run_log()
        logging.info('Tumour growth in progress.')

        start = timer()

        seed=self.__seed
        prng.seed(seed)
        numpy.random.seed(seed)

        #run growth function
        #output variable (true_vaf) is list of tuples with mutation id and frequency of mutation in the tumour [(mut_id, frequency),...]
        true_vaf=self.tumour_growth()

        # Export a graph containing change in tumour size over time
        if self.parameters.plot_tumour_growth:
            self.growth_plot()

        # Sampling
        # Setup list of coordinates that serve as center of sampling [(x,y)]
        samples_coordinates_list=self.__find_sample_coordinates()

        #iterate over each sample from the list of samples
        for center_cell_coordinates in samples_coordinates_list:
            #get sample of certain size
            extended_sample=self.extend_sample(center_cell_coordinates, sample_size=self.parameters.sampling_fraction)

            #extract mutation profiles of all cells found in the sample
            dna_from_sample=self.mutation_reconstruction(extended_sample)

            #count the number of detected mutations and calculate frequency of each mutation (getFrequencies=False gives count for each mutation)
            counted_sample=self.count_mutations(dna_from_sample, get_frequencies=True)

            if self.parameters.number_of_mutations_per_division==1 and self.parameters.number_of_initial_mutations==1:
                #export mutational profile of the sample
                self.export_sample(counted_sample, center_cell_coordinates)


            if self.parameters.number_of_mutations_per_division>1 or self.parameters.number_of_initial_mutations>1:
                #increases number of mutations in the tumour by factor from params.number_of_number_of_mutations_per_division
                increased_mut_number_sample=self.increase_mut_number(counted_sample)
                # Additional mutation serves to distinguish different mutations that occured
                # in the same cell at the same time.

                # Introduce sequencing noise, works only with increased number of mutations
                noisy_data=self.simulate_seq_depth(increased_mut_number_sample)
                self.export_sample(noisy_data, center_cell_coordinates)

                #creates and exports histogram of mutational frequencies
                self.export_histogram(noisy_data, center_cell_coordinates)

        end=timer()

        self.__post_run_log()

        logging.info("Consumed Wall time of this run: %f s.", end - start)

        return 0

    def __pre_run_log(self):
        message = ""

        logging.info("Ready to start CancerSim run with these parameters:")
        for k,v in self.parameters.__dict__.items():
            logging.info("%s = %s", k.split("__")[-1], v)
    
    def __post_run_log(self):
        logging.info("CancerSim run has finished.")
        logging.info("Simulation output written to: %s.", self.__simdir)
        logging.info("Log files written to: %s.""", self.__logdir)

    def export_histogram(self, sample_data, sample_coordinates):
        """ Create and export histogram of mutational frequencies (aka variant allelic frequencies)

        :param sample_data: List of mutations and their frequencies
        :type  sample_data: list

        :param sample_coordinates: coordinates of central sample cell
        :type  sample_coordinates: tuple (i,j) of cell indices
        """



        xaxis_histogram=np.arange(0.0,1,0.01)
        #setdetection limit of the mutation in the sample (depends on the sequencing machine and sequencing depth)
        detection_limit=0.05

        #plots all mutations with frequences above detection threshold
        fig, ax = plt.subplots()
        _ = ax.hist([s[1] for s in sample_data if s[1]>detection_limit], bins=xaxis_histogram)

        _ = ax.set_xlabel('Mutation frequency')
        _ = ax.set_ylabel('Number of mutations')

        #export VAF histogram of the whole tumour
        if sample_coordinates=='whole_tumour':
            figure_path = os.path.join(self.__simdir,'wholeTumourVAFHistogram.pdf')

        #export VAF histogram of sample
        else:
            figure_path = os.path.join(self.__simdir,'sampleHistogram_'+str(sample_coordinates[0])+'_'+str(sample_coordinates[1])+'.pdf')

        fig.savefig(figure_path)
        plt.close(fig)

    def export_sample(self, sample_data, sample_coordinates):
        """ Export (write to disk) frequencies of samples.

        :param sample_data: List of mutations and their frequencies
        :type  sample_data: list

        :param sample_coordinates: coordinates of central sample cell
        :type  sample_coordinates: tuple (i,j) of cell indices
        """

        if len(sample_data) == 0:
            return

        fname = os.path.join(self.__simdir, 'sample_out_'+str(sample_coordinates[0])+'_'+str(sample_coordinates[1])+'.txt')
        logging.info("Writing sampled tumour data to %s.", fname)

        with open(fname,'w') as sample_vaf_ex:
            if len(sample_data[0])==2:
                    sample_vaf_ex.write('mutation_id'+'\t'+'frequency'+'\n')
                    for i in sample_data:
                        sample_vaf_ex.write(str(i[0])+'\t'+str(i[1])+'\n')


            elif len(sample_data[0])==3:
                sample_vaf_ex.write('mutation_id'+'\t'+'additional_mut_id'+'\t'+'frequency'+'\n')
                for i in sample_data:
                    sample_vaf_ex.write(str(i[0])+'\t'+str(i[2])+'\t'+str(i[1])+'\n')

    def export_tumour_matrix(self, tumour_mut_data):
        """ Export (write to disk) the matrix of tumour cells.

        :param tumour_matrix: The tumour matrix to export
        :type  tumour_matrix: array like

        """
        if not self.parameters.export_tumour:
            return 

        fname = os.path.join(self.__simdir, 'mtx_VAF.txt')
        logging.info('Writing tumour profile to %s.', fname)
        
        # save VAF to text file
        with open(fname,'w') as vaf_ex:
            if len(tumour_mut_data[0])==2:
                vaf_ex.write('mutation_id'+'\t'+'frequency'+'\n')
                for i in tumour_mut_data:
                    vaf_ex.write(str(i[0])+'\t'+str(i[1])+'\n')

            if len(tumour_mut_data[0])==3:
                vaf_ex.write('mutation_id'+'\t'+'additional_mut_id'+'\t'+'frequency'+'\n')
                for i in tumour_mut_data:
                    vaf_ex.write(str(i[0])+'\t'+str(i[2])+'\t'+str(i[1])+'\n')

        # Pickle the data.
        fname = os.path.join(self.__simdir, 'mtx.p')
        logging.info('Writing simulation matrix to %s.', fname)
        with open(fname,'wb') as fp:
            pickle.dump(self.__mtx, fp)

        fname = os.path.join(self.__simdir, 'mut_container.p')
        logging.info('Writing mutation list to %s.', fname)
        with open(fname,'wb') as fp:
            pickle.dump(self.__mut_container, fp)

    def growth_plot(self):
    #Plots number of cancer cells over time and outputs it in .pdf
        if self.outdir is None:
            return

        fig, ax = plt.subplots()
        _ = ax.plot([x/2 for x in range(len(self.__growth_plot_data))], self.__growth_plot_data)

        _ = ax.set_xlabel('Division cycle')
        _ = ax.set_ylabel('Number of tumour cells')
        figure_path = os.path.join(self.__simdir,'growthCurve.pdf')
        fig.savefig(figure_path)
        plt.close(fig)

        logging.info("Growth curve graph written to %s.", figure_path)

    def count_mutations(self,mutation_list, get_frequencies):
        """ Count number each time mutation is detected in the sample

        :param mutation_list: mutation profiles of each cell in the sample
        :type  mutation_list: list of lists

        """
        mut_count=[]

        #flatten the list of mutations
        reduced=list(itertools.chain(*[j for j in mutation_list]))

        #count number of unique mutations in whole tumour at time step
        for i in set(reduced):
            mut_count.append((i, float(reduced.count(i))))

        #sort list of mutations based on the mutation id just in case they are not sorted
        mut_count=sorted(mut_count,key=itemgetter(0))
        mut_freq=[]
        if get_frequencies:
            for mutation in mut_count:
                #getting grequency of each mutation by dividing absolute number of detected mutations for each mutation with number of mutations "1" as this one is a proxy for sample size as every cancer cell has mutaiton "1"
                mut_freq.append((mutation[0],(mutation[1]/mut_count[0][1])/self.__ploidy))

            return mut_freq

        return mut_count

    def simulate_seq_depth(self, extended_vaf):
        """ Adds a beta binomial noise to sampled mutation frequencies

        :param extended_vaf: The list of cells to take a sample from.
        :type  extended_vaf: list
        """

        depth=np.random.poisson(self.parameters.read_depth, len(extended_vaf))

        AF=np.array([i[1] for i in extended_vaf])


        samp_alleles=np.random.binomial(depth, AF)

        VAF = samp_alleles/depth

        return [(extended_vaf[i][0], VAF[i], extended_vaf[i][2]) for i in range(len(extended_vaf)) if VAF[i]!=0]

    def increase_mut_number(self, original_mut_list):
        """ Scale up the number of mutations according to the 'number_of_initial_mutations' 'and number_of_mutations_per_division' parameter.

        :param solid_pre_vaf: The list of mutations to scale.
        :type  solid_pre_vaf: list

        """

        extended_mut_list=[]

        target_mut_solid=[]

        for i in original_mut_list:
            #first mutation duplicate N number of times
            # adding additional clonal mutations

            if i[0]==1:
                for j in range(self.parameters.number_of_initial_mutations):
                    extended_mut_list.append((i[0] , float(i[1]),j))

            else:
                # for all subsequent mutations duplicate number
                # of them based on poisson distribution in variable self.__mutMultiplier


                for j in range(self.__mut_multiplier[i[0]]):
                    extended_mut_list.append((i[0] , float(i[1]),j))

        # Return the multiplied mutations.

        return extended_mut_list

    def terminate_cell(self, cell, step):
        """ Kills cancer cell and removes it from the pool of cancer cells

        :param cell: cell chosen for termination
        :type cell: tuple (i,j) of cell indices.

        :param int step: The time step in the simulation
        """
        #removes cell from the pool
        self.__pool.remove(cell)

        #resets value of position on matrix to zero
        self.__mtx[cell]=0

    def death_step(self, step):

        """ Takes a group of random cells and kills them

        :param int step: The time step in the simulation

        """
        for cell in self.__pool:
            beneficial = self.__mtx[cell] in self.__beneficial_mutation
            r = prng.random()

            if (beneficial and r < self.parameters.adv_mutant_death_probability) or r < self.parameters.death_probability:
                self.terminate_cell(cell, step)


    def mutation_reconstruction(self,cells_to_reconstruct):
        """ Reconstructs list of mutations of individual cell by going thorough its ancestors.

        :param cell: Cell for which mutational profile will be recovered.
        :type cell: list of tuples [(i,j)] of cell indices.

        """
        # Return container.
        reconstructed = []

        # Map mutation count to origin (could save this step if elements in
        # mut_container where (c,o) instead of (o,c)).
        lookup_map = dict([(k,v) for v,k in self.__mut_container])

        # Loop over cell indices.
        for i in cells_to_reconstruct:

            # Get cell.
            cell = self.__mtx[i]
            logging.debug("Untangling cell %d.", cell)

            # Setup intermediate container.
            mut_prof=[]

            # Start with the first mutation of this cell.
            mc=self.__mut_container[cell]

            # Get mutation count
            m = mc[1]

            # Now go through the mutation container and trace back the history.
            while m>0:
                # Append current mutation count.
                mut_prof.append(m)

                # Get mutation origin of this count.
                m = lookup_map[m]

            # Store on return container in reverse order.
            reconstructed.append(mut_prof[::-1])

        return reconstructed

    def tumour_growth(self):
        """ Run the tumour growth simulation.  """

        # setup a counter to keep track of number of mutations that occur in this run.
        # take into account mutations from previous runs if rerun.
        if self.parameters.tumour_multiplicity == 'single':
            mutation_counter = self.__mutation_counter
        if self.parameters.tumour_multiplicity == 'double':
            mutation_counter = self.__mutation_counter + 1

        # Loop over time steps.
        for step in range(self.__init_step, self.__init_step+self.parameters.number_of_generations):
            # logging.debug("Cell matrix: \n%s", str(self.__mtx.todense()))
            logging.debug('%d/%d generation started', step+1, self.__init_step + self.parameters.number_of_generations + 1)

            # setup a temporary list to store the mutated cells in this iteration.
            temp_pool=[]

            # reshuffle the order of pool to avoid that cells with low number divide always first.
            shuffle(self.__pool)

            # logging.debug('list of cancer cells %s', str(self.__pool))

            # Loop over all cells in the pool.
            for cell in self.__pool:
                logging.debug('cell to divide %s', str(cell))

                # Get the existing neighboring cells.
                neigh=self.neighbours(cell)

                # first condition: if available neighbors
                if neigh:
                    # if cell has beneficial mutation.
                    beneficial = self.__mtx[cell] in self.__beneficial_mutation
                    r = prng.random()
                    if (beneficial and r < self.parameters.adv_mutant_division_probability) or r < self.parameters.division_probability:
                        mutation_counter = self.division(cell, beneficial, neigh, step, mutation_counter, temp_pool)


            # add new cancer cells to a pool of cells available for division next round
            [self.__pool.append(v) for v in temp_pool]

            self.__growth_plot_data.append(len(self.__pool))

            self.death_step(step)
            self.__growth_plot_data.append(len(self.__pool))

        logging.info("All generations finished. Starting tumour reconstruction.")


        # Update internal step counter if we dump and reload.
        self.__init_step = step+1
        self.__mutation_counter = mutation_counter

        # Reconstruct mutation history.
        reconstructed = self.mutation_reconstruction(self.__pool)

        logging.info("Reconstruction done,  get statistics.")

        mutation_counts=self.count_mutations(reconstructed, get_frequencies=True)

        if self.parameters.number_of_mutations_per_division==1 and self.parameters.number_of_initial_mutations==1:
            self.export_tumour_matrix(mutation_counts)
            return mutation_counts

        if self.parameters.number_of_mutations_per_division>1 or self.parameters.number_of_initial_mutations>1:

            increased_mut_number_tumour=self.increase_mut_number(mutation_counts)    #increases number of mutations in the tumour by factor from params.number_of_number_of_mutations_per_division

            noisy_data=self.simulate_seq_depth(increased_mut_number_tumour)       #introduce sequencing noise, works only with increased number of mutations

            self.export_tumour_matrix(noisy_data)
            center_cell_coordinates='whole_tumour'
            self.export_histogram(noisy_data, center_cell_coordinates)       #creates and exports histogram of mutational frequencies
            return noisy_data



    def neighbours(self, cell):
        """ Returns the nearest-neighbor cells around the given node.

        :param cell: The node for which to calculate the neighbors.
        :type  cell: tuple (i,j) of cell indices.

        """
        # make list of all surrounding nodes
        neighboursList=[
                (cell[0]-1, cell[1]+1),
                (cell[0]  , cell[1]+1),
                (cell[0]+1, cell[1]+1),
                (cell[0]-1, cell[1]  ),
                (cell[0]+1, cell[1]  ),
                (cell[0]-1, cell[1]-1),
                (cell[0]  , cell[1]-1),
                (cell[0]+1, cell[1]-1)]

        # return nodes that are not cancerous, do not contain mutation index
        return [y for y in neighboursList if self.__mtx[y]==0]

    def place_to_divide(self):
        """ Selects random unoccupied place on the matrix where cell will divide."""

        a = prng.randint(0,self.parameters.matrix_size-1)
        b = prng.randint(0,self.parameters.matrix_size-1)
        random_place_to_divide=(a,b)

        if self.__mtx[random_place_to_divide]==0:
            return a, b

        else:
            while self.__mtx[random_place_to_divide]!=0:
                a = prng.randint(0,self.parameters.matrix_size-1)
                b = prng.randint(0,self.parameters.matrix_size-1)
                random_place_to_divide=(a,b)

            return a, b

    def division(self, cell, beneficial, neighbors, step, mutation_counter, pool):
        """ Perform a cell division.

        :param tuple cell: The mother cell coordinates.
        :param bool beneficial: Flag to indicate if the cell carries the beneficial mutation.
        :param list neighbors: The neighboring cells.
        :param int step: The time step in the simulation
        :param int mutation_counter: The counter of mutations to be updated
        :param list pool: The (temporary) pool of cells.

        """
        # Draw a free neighbor.
        place_to_divide=prng.choice(neighbors)
        pool.append(place_to_divide)

        if beneficial:
            logging.info("Division of beneficial mutation carrier. Cell index = %s, mutation index = %d, place_to_divide=%s", str(cell), self.__mtx[cell], str(place_to_divide)) 

        mutation_counter = self.mutation(cell, neighbors, step, mutation_counter, pool,place_to_divide, beneficial)

        return mutation_counter

    def mutation(self, *args):
        """ Perform a mutation.

        :param cell: At which cell the mutation occurs
        :param neighbors: The neighboring cells
        :param mutation_counter: The current number of mutations, to be incremented.
        :param pool: The pool of all cells.
        :param place_to_divide: The position at which the mutation occurs.
        :param beneficial: Flag to control whether the mutation is beneficial or not.

        """

        cell, neighbors, step, mutation_counter, pool, place_to_divide, beneficial = args

        # Mutation.
        if prng.random()<self.parameters.mutation_probability:

            # Increment mutation counter.
            mutation_counter=mutation_counter+1

            # New cell gets the index number of largest number of mutation
            self.__mtx[place_to_divide]=len(self.__mut_container)
            self.__mut_container.append((self.__mut_container[self.__mtx[cell]][1], mutation_counter))

            # Log
            logging.debug('Neighbor cell has new index %d', self.__mtx[place_to_divide])
            logging.debug("%d, %d", self.__mut_container[self.__mtx[cell]][1], mutation_counter)
            # logging.debug('mut container updated: %s', str(self.__mut_container))

            if beneficial:
                self.__beneficial_mutation.append(int(self.__mtx[place_to_divide]))
                logging.info("Mutation  of beneficial mutation carrier. Cell index = %s, mutation index = %d, place to divide = %s", cell, self.__mtx[cell], str(place_to_divide))

            else:
                # Decide whether an advantageous mutation occurs.
                if prng.random()<self.parameters.adv_mutant_mutation_probability \
                        and len(self.__beneficial_mutation)==0 \
                        and step==self.parameters.adv_mutation_wait_time:
                    logging.info('New beneficial mutation: %d', int(self.__mtx[place_to_divide]))
                    self.__beneficial_mutation.append(int(self.__mtx[place_to_divide]))

            # Mother cell mutates
            mutation_counter=mutation_counter+1
            self.__mut_container.append((self.__mut_container[self.__mtx[cell]][1], mutation_counter))

            # Update mutation list.
            self.__mtx[cell]=len(self.__mut_container)-1

        # No new mutation.
        else:
            logging.debug('No new mutation in normal division, inheriting from parent')
            self.__mtx[place_to_divide]=self.__mtx[cell]

        return mutation_counter

    def __find_sample_coordinates(self):
        """ """
        """ Find the sample coordinates based on the list of coordinates given
        at startup. """

        # If no positions where given, find the center of the matrix.
        if self.parameters.sampling_positions is None:
            self.parameters.sampling_positions = [prng.choice(self.__pool)]

        return self.parameters.sampling_positions

def main(arguments):
    """ The entry point for the command line interface.

    :param arguments: The command line arguments for the cancer simulation tool.
    :type  arguments: Namespace

    """

    parameters = CancerSimulatorParameters()

    if os.path.isfile(arguments.params):

        spec = spec_from_file_location("params", arguments.params)
        params = module_from_spec(spec)
        spec.loader.exec_module(params)

        parameters = CancerSimulatorParameters(
                matrix_size=params.matrix_size,
                number_of_generations=params.number_of_generations,
                division_probability=params.division_probability,
                adv_mutant_division_probability=params.adv_mutant_division_probability,
                death_probability=params.death_probability,
                adv_mutant_death_probability=params.adv_mutant_death_probability,
                mutation_probability=params.mutation_probability,
                adv_mutant_mutation_probability=params.adv_mutant_mutation_probability,
                number_of_mutations_per_division=params.number_of_mutations_per_division,
                adv_mutation_wait_time=params.adv_mutation_wait_time,
                number_of_initial_mutations=params.number_of_initial_mutations,
                tumour_multiplicity=params.tumour_multiplicity,
                sampling_fraction=params.sampling_fraction,
                sampling_positions=params.sampling_positions,
                read_depth=params.read_depth,
                export_tumour=params.export_tumour,
                plot_tumour_growth=params.plot_tumour_growth,
                )

    # Set loglevel.
    loglevel = {0 : logging.WARNING,
                1 : logging.INFO,
                2 : logging.DEBUG,
                }

    if not arguments.loglevel in loglevel.keys():
        arguments.loglevel = 0

    logging.getLogger().setLevel(loglevel[arguments.loglevel])
    casim = CancerSimulator(parameters, seed=arguments.seed, outdir=arguments.outdir)

    return (casim.run())


def check_set_number(value, typ, default=None, minimum=None, maximum=None):
    """ Checks if a value is instance of type and lies within permissive_range if given. """

    if value is None:
        return default

    if not isinstance(value, typ):
        try:
            value = typ(value)
        except:
            raise TypeError("Incompatible type: Expected {0}, got {1}.".format(typ, type(value)))

    if minimum is not None:
        if value < minimum:
            raise ValueError("Value must be larger than {}.".format(minimum))

    if maximum is not None:
        if value > maximum:
            raise ValueError("Value must be smaller than {}.".format(maximum))

    return value


def load_cancer_simulation(dumpfile):
    """ Unpickle a cancer simulation from a dill generated dump.

    :param dumpfile: Path to the file that contains the dumped object.
    :type  dumpfile: str

    """

    with open(dumpfile, 'rb') as fp:
        obj = dill.load(fp)

    return obj


if __name__ == "__main__":
    # Entry point
    # Setup the command line parser.
    parser = ArgumentParser()

    # Seed parameter.
    parser.add_argument("-p",
                        "--params",
                        help="""Path to the python file holding the simulation
                        parameters. Defaults to `params.py` in the current working
                        directory. If no file is found, default parameters will
                        be chosen, see the API for CancerSimulatorParameters for
                        details.""",
                        default="params.py",
                        type=str,
                        metavar="PARAMS",
                        )

    parser.add_argument("-s",
                        "--seed",
                        help="The prng seed.",
                        type=int,
                        default=1,
                        metavar="SEED",
                        )

    parser.add_argument("-o",
                        "--outdir",
                        dest="outdir",
                        metavar="DIR",
                        default=None,
                        help="Directory where simulation data is saved.",
                        type=str,
                        )

    parser.add_argument("--verbose",
                        "-v",
                        dest="loglevel",
                        action='count',
                        default=0,
                        help="Increase the verbosity level by adding 'v's."
                        )

    # Parse the arguments.
    arguments = parser.parse_args()

    sys.exit(main(arguments))

