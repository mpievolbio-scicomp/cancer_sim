# -*- coding: utf-8 -*-
#!/usr/bin/env python3
__author__ = 'Luka Opasic, MD'
__email__ = 'opasic@evolbio.mpg.de'
__version__ = '0.0.1'


from argparse import ArgumentParser
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
from operator import itemgetter

np = numpy

# Configure logging.
LEVEL = logging.WARNING
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(LEVEL)
HANDLER = logging.StreamHandler()
HANDLER.setLevel(LEVEL)
FORMATTER=logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
HANDLER.setFormatter(FORMATTER)
LOGGER.addHandler(HANDLER)


class CancerSimulatorParameters(object):
    """
        :class CancerSimulatorParameters: Represents the parameters for a cancer simulation.
    """

    def __init__(self,
                 matrix_size                         = None,
                 number_of_generations               = None,
                 division_probability                = None,
                 advantageous_division_probability   = None,
                 death_probability                   = None,
                 advantageous_death_probability = None,
                 mutation_rate                       = None,
                 advantageous_mutation_probability   = None,
                 mutations_per_division              = None,
                 time_of_advantageous_mutation       = None,
                 number_of_clonal                    = None,
                 tumour_multiplicity                 = None,
                ):
        """
        Construct a new CancerSimulationParameters object.

        :param matrix_size: The size of the grid in each dimension.
        :type  matrix_size: int

        :param number_of_generations: The number of generations to simulate.
        :type  number_of_generations: int

        :param division_probability: The probability for a cell division to occur during one generation.
        :type  division_probability: float (0.0 <= division_probability <= 1.0)

        :param advantageous_division_probability: The probability for the division of a cell with advantageous mutation to occur during one generation.
        :type  advantageous_division_probability: float (0.0 <= division_probability <= 1.0)

        :param death_probability: The probability for a cell to die during one generation.
        :type  death_probability: float (0.0 <= division_probability <= 1.0)

        :param advantageous_death_probability: The probability for a cell with advantageous mutation to die during one generation.
        :type  advantageous_death_probability: float (0.0 <= division_probability <= 1.0)

        :param mutation_rate: The rate of mutation (probability per generation).
        :type  mutation_rate: float (0.0 <= division_probability <= 1.0)

        :param advantageous_mutation_probability: The rate for an advantageous mutation to occur during one generation.
        :type  advantageous_mutation_probability: float (0.0 <= division_probability <= 1.0)

        :param mutations_per_division: The number of mutations per division
        :type  mutations_per_division: int

        :param time_of_advantageous_mutation: The number of generations after which an advantageous mutation can occur.
        :type  time_of_advantageous_mutation: int

        :param number_of_clonal: Scale up each mutation by this factor.
        :type  number_of_clonal: int

        :param tumour_multiplicity: Run in single or double tumour mode. Possible values: "single", "double".
        :type  tumour_multiplicity: str
        """

        # Store parameters on the object.
        self.matrix_size = matrix_size
        self.number_of_generations = number_of_generations
        self.division_probability = division_probability
        self.advantageous_division_probability = advantageous_division_probability
        self.death_probability = death_probability
        self.advantageous_death_probability = advantageous_death_probability
        self.mutation_rate = mutation_rate
        self.advantageous_mutation_probability = advantageous_mutation_probability
        self.mutations_per_division = mutations_per_division
        self.time_of_advantageous_mutation = time_of_advantageous_mutation
        self.number_of_clonal = number_of_clonal
        self.tumour_multiplicity = tumour_multiplicity

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
    def advantageous_division_probability(self):
        return self.__advantageous_division_probability
    @advantageous_division_probability.setter
    def advantageous_division_probability(self, val):
        self.__advantageous_division_probability = check_set_number(val, float, 1, 0.0, 1.0)

    @property
    def death_probability(self):
        return self.__death_probability
    @death_probability.setter
    def death_probability(self, val):
        self.__death_probability = check_set_number(val, float, 0, 0.0, 1.0)

    @property
    def advantageous_death_probability(self):
        return self.__advantageous_death_probability
    @advantageous_death_probability.setter
    def advantageous_death_probability(self, val):
        self.__advantageous_death_probability = check_set_number(val, float, 0.0, 0.0, 1.0)

    @property
    def mutation_rate(self):
        return self.__mutation_rate
    @mutation_rate.setter
    def mutation_rate(self, val):
        self.__mutation_rate = check_set_number(val, float, 0.8, 0.0, 1.0)

    @property
    def advantageous_mutation_probability(self):
        return self.__advantageous_mutation_probability
    @advantageous_mutation_probability.setter
    def advantageous_mutation_probability(self, val):
        self.__advantageous_mutation_probability = check_set_number(val, float, 1.0, 0.0, 1.0)

    @property
    def mutations_per_division(self):
        return self.__mutations_per_division
    @mutations_per_division.setter
    def mutations_per_division(self, val):
        self.__mutations_per_division = check_set_number(val, int, 1, 0)

    @property
    def time_of_advantageous_mutation(self):
        return self.__time_of_advantageous_mutation
    @time_of_advantageous_mutation.setter
    def time_of_advantageous_mutation(self, val):
        self.__time_of_advantageous_mutation = check_set_number(val, int, 50000, 0)

    @property
    def number_of_clonal(self):
        return self.__number_of_clonal
    @number_of_clonal.setter
    def number_of_clonal(self, val):
        self.__number_of_clonal = check_set_number(val, int, 1, 0)

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

        :param outdir: The directory where simulation data is saved.
        :type  outdir: (str || path-like object)
        """

        if parameters is None:
            raise ValueError("No parameters given, simulation cannot execute.")
        self.parameters = parameters

        # Setup internal variables.
        self.__prop_of_driver = None
        self.__mtx = lil_matrix((self.parameters.matrix_size, self.parameters.matrix_size), dtype=int)
        self.__mut_container = None
        self.__lq_bipsy = None
        self.__xaxis_histogram = None
        self.__death_list = None
        self.__biopsy_timing = None
        self.__beneficial_mutation = []
        self.__growth_plot_data = None
        self.__mutation_counter = None
        self.__s = self.parameters.mutations_per_division
        self.__export_tumour = True
        self.__export_tumour_growth = False
        self.__tumour_multiplicity = self.parameters.tumour_multiplicity

        # Handle direct parameters.
        self.seed = seed
        self.outdir = outdir
        

        self.__ploidy=2

        
        self.__mutMultiplier=[self.__s]*100000
        #self.__mutMultiplier=np.random.poisson(self.__s, 1000)
        
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

        self.__outdir = '../outdir'
        self.__seeddir = '../seeddir'
        self.__logdir  = '../logdir'
        self.__simdir = '../simdir'

        if outdir is None:
            return

        # Not None, so we want to store output. Set flag accordingly.
        self.__export_tumour = True

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


    def extendSample(self, sampleCenter, sampleSize):
        """ Takes a subset of cells from the tumour positioned around single input cell with specific coordinates. Output is list containing tuples of cells belonging to the sample.
        :param sampleCenter: coordinates of cell that will be center of the sample
        :type  sample: tuple
        """
        
        biopsy_size=math.ceil(sampleSize*len(self.__pool))
                
        if biopsy_size>1:
            for z in range(1,len(self.__pool)):  #look at z-tier neighbours around sampleCenter
                expandedSample=[]
                for i in range(-z, z+1):  #look at z-tier neighbours around sampleCenter
                    for j in range(-z, z+1):
                        nc=(sampleCenter[0]-i, sampleCenter[1]+j)
                        if nc in self.__pool: #if surrounding cell in the pool add it to sample list
                            expandedSample.append(nc)
                if len(expandedSample)>biopsy_size:  # if chunk is larger than wanted percentage of total tumour
                    while len(expandedSample)>biopsy_size: #remove last value until desired chunk size
                        expandedSample=expandedSample[:-1]
                    break
           
            return expandedSample
        else:
            return [sampleCenter]


    def dump(self):
        """ Serialize the object. The current simulation will be stored in a machine readable
        format to <OUTDIR>/cancer_<SEED>/cancer_sim.py.dill, where <OUTDIR> is the specified output
        directory or (if the latter was not defined) a temporary directory."""

        with open(self.dumpfile, 'wb') as fp:
            dill.dump(self, fp)

    def run(self):
        """ Run the simulation. """

        # Setup square matrix.
        matrix_size=self.parameters.matrix_size

        if self.__tumour_multiplicity == 'single':
            LOGGER.info('Running in single tumour mode.')
            initLoc=(int(matrix_size*0.5),int(matrix_size*0.5))

            LOGGER.info("First cell at %s.", str(initLoc))

            self.__mtx[initLoc]=1
            self.__mut_container=[(0, 0), (0, 1)]
            self.__pool=[initLoc]

        #start the pool of cancer cells by adding the initial cancer cell into it
        if self.__tumour_multiplicity == 'double':
            LOGGER.info('Running in sdsa mode.')

            ### COMMENT: Should these be given as parameters?
            distance_between_tumours=0.05
            initLoc=(int(matrix_size*0.45),int(matrix_size*0.5))
            secondinitLoc=(int(matrix_size*0.65),int(matrix_size*0.51))

            self.__mtx[initLoc]=1
            self.__mtx[secondinitLoc]=2
            self.__mut_container=[(0, 0), (0, 1), (0,2)]
            self.__pool=[initLoc, secondinitLoc]

        # create lists used in loops
        lq_bipsy=[]
        self.__growth_plot_data=[]

        LOGGER.info('Tumour growth in progress.')

        death_list=[]
        prop_of_driver=[]

        start = timer()

        seed=self.__seed
        prng.seed(seed)

        #run growth function
        #output variable (true_vaf) is list of tuples with mutation id and frequency of mutation in the tumour [(mut_id, frequency),...]
        true_vaf=self.tumourGrowth()
        
        
        ##Exports
        if self.__export_tumour is True:   #export tumour as pickled matrix together with accompanying data
            self.export_tumour_matrix(true_vaf)

       
        if self.__export_tumour_growth is True:   #export a graph containing change in tumour size over time
            self.growth_plot()
        

#############here sampling
        
        #samplesCoordinatesList=[self.__pool[5],self.__pool[3],self.__pool[10]] #list of coordinates that serve as center of sampling [(x,y),(x,y),(x,y)]
        samplesCoordinatesList=[self.__pool[5]] #list of coordinates that serve as center of sampling [(x,y),(x,y),(x,y)]

        for centerCellCoordinates in samplesCoordinatesList:          
            extendedSample=self.extendSample(centerCellCoordinates, sampleSize=0.1)    #get sample of certain size
            
            dna_from_sample=self.mutation_reconstruction(extendedSample)

            countedSample=self.countMutations(dna_from_sample, getFrequencies=True)
            self.export_sample(countedSample, centerCellCoordinates)
            

            ##optional, increase mutation number 
            increasedMutNumberSample=self.increase_mut_number(countedSample)
            

            #introduce sequencing noise, works only with increased number of mutations
            noisyData=self.simulate_seq_depth(increasedMutNumberSample)
            self.export_histogram(noisyData, centerCellCoordinates)
                       
            #self.export_sample(countedSample, centerCellCoordinates)

            

#############


        end=timer()
        LOGGER.info("Consumed Wall time of this run: %f s.", end - start)


    def export_histogram(self, sampleData, sampleCoordinates):    
        xaxis_histogram=np.arange(0.0,1,0.01)   
        detectionLimit=0.05
              

        plt.hist([s[1] for s in sampleData if s[1]>detectionLimit], bins=xaxis_histogram)
        
        plt.xlabel('Mutation frequency')
        plt.ylabel('Number of mutations')
        
        figure_path = os.path.join(self.__outdir,'sampleHistogram_'+str(sampleCoordinates[0])+'_'+str(sampleCoordinates[1])+'.pdf')
        plt.savefig(figure_path)
        
        #plt.show()

    def export_sample(self, sampleData, sampleCoordinates):
        #  """ Export (write to disk) frequencies of samples.

        # :param sampleData: List of mutations and their frequencies
        # :type  sampleData: list
        
        # :param sampleCoordinates: coordinates of central sample cell
        # :type  sampleCoordinates: tuple (i,j) of cell indices 
        # """

        with open(os.path.join(self.__outdir, 'sample_out_'+str(sampleCoordinates[0])+'_'+str(sampleCoordinates[1])+'.txt'),'w') as sampleVafEx:
            sampleVafEx.write('mutation_id'+'\t'+'frequency'+'\n')
            for i in sampleData:
                sampleVafEx.write(str(i[0])+'\t'+str(i[1])+'\n')
                

    def export_tumour_matrix(self, tumour_matrix):
        """ Export (write to disk) the matrix of tumour cells.

        :param tumour_matrix: The tumour matrix to export
        :type  tumour_matrix: array like

        """

        LOGGER.info('Exporting simulation data')
        
        # save VAF to text file
        with open(os.path.join(self.__simdir, 'mtx_VAF.txt'),'w') as vafEx:
            vafEx.write('mutation_id'+'\t'+'frequency'+'\n')
            for i in tumour_matrix:
                vafEx.write(str(i[0])+'\t'+str(i[1])+'\n')
                #vafEx.write(str(i)+'\n')
        
        # Pickle the data.
        with open(os.path.join(self.__simdir, 'mtx.p'),'wb') as fp:
            pickle.dump(self.__mtx, fp)
        with open(os.path.join(self.__simdir, 'mut_container.p'),'wb') as fp:
            pickle.dump(self.__mut_container, fp)
        with open(os.path.join(self.__simdir, 'death_list.p'),'wb') as fp:
            pickle.dump(self.__death_list, fp)

    def growth_plot(self):
    #Plots number of cancer cells over time and outputs it in .pdf

        if self.outdir is None:
            return

        plt.plot([x/2 for x in range(len(self.__growth_plot_data))], self.__growth_plot_data)

        plt.xlabel('Division cycle')
        plt.ylabel('Number of tumour cells')
        figure_path = os.path.join(self.__simdir,'growthCurve.pdf')
        plt.savefig(figure_path)

        LOGGER.info("Growth curve graph written to %s.", figure_path)

        plt.clf()


    def countMutations(self,mutationList, getFrequencies):
        """ Count number each time mutation is detected in the sample 

        :param mutationList: mutation profiles of each cell in the sample
        :type  mutationList: list of lists 

        """
        mutCount=[]
                

        reduced=list(itertools.chain(*[j for j in mutationList]))  #flatten the list of mutations
                
        for i in set(reduced): #count number of unique mutations in whole tumour at time step
            mutCount.append((i, float(reduced.count(i))))
        

        mutCount=sorted(mutCount,key=itemgetter(0)) #sort list of mutations based on the mutation id just in case they are not sorted
        mutFreq=[]
        if getFrequencies==True:
            for mutation in mutCount:
                mutFreq.append((mutation[0],(mutation[1]/mutCount[0][1])/self.__ploidy)) #getting grequency of each mutation by dividing absolute number of detected mutations for each mutation with number of mutations "1" as this one is a proxy for sample size as every cancer cell has mutaiton "1" 
            
            return mutFreq

    
        return mutCount

 
   
    def simulate_seq_depth(self, extended_vaf):
        """ Ads a beta binomial noise to sampled mutation frequencies

        :param extended_vaf: The list of cells to take a sample from.
        :type  extended_vaf: list
        """
        
        #depth=np.random.poisson(params.read_depth, len(input_mut))
        depth=np.random.poisson(100, len(extended_vaf))   
        
        AF=np.array([i[1] for i in extended_vaf])
        
           
        samp_alleles=np.random.binomial(depth, AF)
        
        VAF = samp_alleles/depth
       
        return [(extended_vaf[i][0], VAF[i], extended_vaf[i][2]) for i in range(len(extended_vaf))]



    def increase_mut_number(self, originalMutList):
        """ Scale up the number of mutations according to the 'number_of_clonal' 'and mut_per_division' parameter.

        :param solid_pre_vaf: The list of mutations to scale.
        :type  solid_pre_vaf: list

        """
        
        extemdedMutList=[]
        
        target_mut_solid=[]

        for i in originalMutList:
            #first mutation duplicate N number of times
            # adding additional clonal mutations
          
            if i[0]==1:
                for j in range(self.parameters.number_of_clonal):
                    extemdedMutList.append((i[0] , float(i[1]),j))
               
            else:
                # for all subsequent mutations duplicate number
                # of them based on poisson distribution in variable self.__mutMultiplier
         
                
                for j in range(self.__mutMultiplier[i[0]]):
                    extemdedMutList.append((i[0] , float(i[1]),j))
                    # Log.
                    #LOGGER.debug("Will append %d additional mutations." % (number_of_additional_mutations))

                    # Append.
                    #target_mut_solid += [(i[0], float(i[1]))]*number_of_additional_mutations

        # Return the multiplied mutations.
        
        return extemdedMutList

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

    def deathStep(self, step):

        """ Takes a group of random cells and kills them

        :param int step: The time step in the simulation

        """
        for i in prng.sample(self.__pool, math.floor(self.parameters.death_probability*len(self.__pool))):
            self.terminate_cell(i, step)
    

    def mutation_reconstruction(self,cellsToReconstruct):
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
        for i in cellsToReconstruct:

            # Get cell.
            cell = self.__mtx[i]
            LOGGER.debug("Untangling cell %d.", cell)

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

    def tumourGrowth(self):
        """ Run the tumour growth simulation.  """

        # setup a counter to keep track of number of mutations that occur in this run.
        if self.__tumour_multiplicity == 'single':
            mutation_counter=1
        if self.__tumour_multiplicity == 'double':
            mutation_counter=2

        # Loop over time steps.
        for step in range(self.parameters.number_of_generations):
            LOGGER.debug("Cell matrix: \n%s", str(self.__mtx.todense()))
            LOGGER.debug('%d/%d generation started', step, self.parameters.number_of_generations)
            
            # setup a temporary list to store the mutated cells in this iteration.
            temp_pool=[]

            # reshuffle the order of pool to avoid that cells with low number divide always first.
            shuffle(self.__pool)

            LOGGER.debug('list of cancer cells %s', str(self.__pool))

            # Loop over all cells in the pool.
            for cell in self.__pool:
                LOGGER.debug('cell to divide %s', str(cell))

                # Get the existing neighboring cells.
                neigh=self.neighbours(cell)

                # first condition: if available neighbors
                if neigh:
                    # if cell has beneficial mutation.
                    if self.__mtx[cell] in self.__beneficial_mutation:
                        # cell divides with greater probability.
                        if prng.random()<self.parameters.advantageous_division_probability:
                            mutation_counter = self.division(cell, True, neigh, step, mutation_counter, temp_pool)

                    # cell does not have beneficial mutation -> normal division.
                    else:
                        if prng.random()<self.parameters.division_probability:
                            mutation_counter = self.division(cell, False, neigh, step, mutation_counter, temp_pool)

            # add new cancer cells to a pool of cells available for division next round
            [self.__pool.append(v) for v in temp_pool]
            self.__growth_plot_data.append(len(self.__pool))

            self.deathStep(step)
            self.__growth_plot_data.append(len(self.__pool))

            # at the end reconstruct mutational frequencies from the whole tumour
            
            if step == self.parameters.number_of_generations-1:

                LOGGER.info("All generations finished. Starting tumour reconstruction.")
                reconstructed = self.mutation_reconstruction(self.__pool)
                
                LOGGER.info("Reconstruction done,  get statistics.")
               
                mutationCounts=self.countMutations(reconstructed, getFrequencies=True)
                
                
                #LOGGER.info("Statistics done,  amplify mutations.")
                                                
                LOGGER.debug('Head of bulk_vaf: %s', str(mutationCounts[0:10]))

                return mutationCounts

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
        if beneficial:

            # Draw a free neighbor.
            place_to_divide=prng.choice(neighbors)
            pool.append(place_to_divide)
            mutation_counter = self.mutation(cell, neighbors, step, mutation_counter, pool,place_to_divide, True)

            return mutation_counter

        else:
            # Draw a free neighbor.
            place_to_divide=prng.choice(neighbors)
            pool.append(place_to_divide)

            LOGGER.debug('index of the mother cell: %s', str(self.__mtx[cell]))
            LOGGER.debug('random neighbor to divide: %s', str(place_to_divide))

            mutation_counter = self.mutation(cell, neighbors, step, mutation_counter, pool, place_to_divide, False)

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
        if prng.random()<self.parameters.mutation_rate:

            # Increment mutation counter.
            mutation_counter=mutation_counter+1

            # New cell gets the index number of largest number of mutation
            self.__mtx[place_to_divide]=len(self.__mut_container)
            self.__mut_container.append((self.__mut_container[self.__mtx[cell]][1], mutation_counter))

            # Log
            LOGGER.debug('Neighbor cell has new index %d', self.__mtx[place_to_divide])
            LOGGER.debug("%d, %d", self.__mut_container[self.__mtx[cell]][1], mutation_counter)
            LOGGER.debug('mut container updated: %s', str(self.__mut_container))

            if beneficial:
                self.__beneficial_mutation.append(int(self.__mtx[place_to_divide]))

            else:
                # Decide whether an advantageous mutation occurs.
                ### FIXME: if there is no "normal mutation" then cell will never get advantageous mutation. if adv_mut_prob is 1 at turn x,  mutation it should happen
                if prng.random()<self.parameters.advantageous_mutation_probability \
                        and len(self.__beneficial_mutation)==0 \
                        and step==self.parameters.time_of_advantageous_mutation:
                    LOGGER.info('new beneficial mutation: %d', int(self.__mtx[place_to_divide]))
                    self.__beneficial_mutation.append(int(self.__mtx[place_to_divide]))

            # Mother cell mutates
            mutation_counter=mutation_counter+1
            self.__mut_container.append((self.__mut_container[self.__mtx[cell]][1], mutation_counter))

            # Update mutation list.
            self.__mtx[cell]=len(self.__mut_container)-1

        # No new mutation.
        else:
            LOGGER.info('No new mutation in normal division, inheriting from parent')
            self.__mtx[place_to_divide]=self.__mtx[cell]

        return mutation_counter

def main(arguments):
    """ The entry point for the command line interface.

    :param arguments: The command line arguments for the cancer simulation tool.
    :type  arguments: Namespace

    """

    parameters = CancerSimulatorParameters()

    if "params.py" in os.listdir(os.getcwd()):

        sys.path.insert(0, os.getcwd())
        import params

        # Catch legacy issue
        ms = None
        if hasattr(params, "matrixSize"):
            ms = params.matrixSize
        else:
            ms = params.matrix_size

        parameters = CancerSimulatorParameters(matrix_size = ms,
                number_of_generations = params.num_of_generations,
                division_probability = params.div_probability,
                advantageous_division_probability = params.fittnes_advantage_div_prob,
                death_probability = params.dying_fraction,
                advantageous_death_probability = params.fitness_advantage_death_prob,
                mutation_rate = params.mut_rate,
                advantageous_mutation_probability = params.advantageous_mut_prob,
                mutations_per_division = params.mut_per_division,
                time_of_advantageous_mutation = params.time_of_adv_mut,
                number_of_clonal = params.num_of_clonal,
                tumour_multiplicity = params.tumour_multiplicity,
                )

    # Set loglevel.
    loglevel = {0 : logging.WARNING,
                1 : logging.INFO,
                2 : logging.DEBUG,
                }

    if not arguments.loglevel in loglevel.keys():
        arguments.loglevel = 0

    LOGGER.setLevel(loglevel[arguments.loglevel])
    HANDLER.setLevel(loglevel[arguments.loglevel])

    casim = CancerSimulator(parameters, seed=arguments.seed, outdir=arguments.outdir)

    casim.run()


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
    parser.add_argument("seed",
                        help="The prng seed.",
                        type=int,
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

    main(arguments)

