# -*- coding: utf-8 -*-
#!/usr/bin/env python3
__author__ = 'Luka Opasic, MD'
__email__ = 'opasic@evolbio.mpg.de'
__version__ = '0.0.1'

from argparse import ArgumentParser
from fractions import gcd
from matplotlib.lines import Line2D
from random import shuffle
from scipy.sparse import lil_matrix
from scipy.spatial import distance
from time import sleep
from timeit import default_timer as timer
import gc
import itertools
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import params
import pickle
import random as prng
import scipy.sparse as sp
import seaborn as sns
import subprocess
import sys


class CancerSimulatorParameters():
    """
        :class CancerSimulatorParameters: Represents the parameters for a cancer simulation.
    """

    def __init__(self,
                 matrix_size                         = None,
                 number_of_generations               = None,
                 division_probability                = None,
                 advantageous_division_probability   = None,
                 death_probability                   = None,
                 fitness_advantage_death_probability = None,
                 mutation_rate                       = None,
                 advantageous_mutation_probability   = None,
                 mutations_per_division              = None,
                 time_of_advantageous_mutation       = None,
                 number_of_clonal                    = None,
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

        :param fitness_advantage_death_probability: The probability for a cell with advantageous mutation to die during one generation.
        :type  fitness_advantage_death_probability: float (0.0 <= division_probability <= 1.0)

        :param mutation_rate: The rate of mutation (probability per generation).
        :type  mutation_rate: float (0.0 <= division_probability <= 1.0)

        :param advantageous_mutation_probability: The rate for an advantageous mutation to occur during one generation.
        :type  advantageous_mutation_probability: float (0.0 <= division_probability <= 1.0)

        :param mutations_per_division: The number of mutations per division
        :type  mutations_per_division: int

        :param time_of_advantageous_mutation: The number of generations after which an advantageous mutation occurs.
        :type  time_of_advantageous_mutation: int

        """
        # Store parameters on the object.
        self.matrix_size = matrix_size
        self.number_of_generations = number_of_generations
        self.division_probability = division_probability
        self.advantageous_division_probability = advantageous_division_probability
        self.death_probability = death_probability
        self.fitness_advantage_death_probability = fitness_advantage_death_probability
        self.mutation_rate = mutation_rate
        self.advantageous_mutation_probability = advantageous_mutation_probability
        self.mutations_per_division = mutations_per_division
        self.time_of_advantageous_mutation = time_of_advantageous_mutation
        self.number_of_clonal = number_of_clonal

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
    def fitness_advantage_death_probability(self):
        return self.__fitness_advantage_death_probability
    @fitness_advantage_death_probability.setter
    def fitness_advantage_death_probability(self, val):
        self.__fitness_advantage_death_probability = check_set_number(val, float, 0.0, 0.0, 1.0)

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

class CancerSimulator(object):
    """
        :class CancerSimulator: Represents the Monte-Carlo simulation of cancer tumour growth on a 2D(3D) grid.
    """

    def __init__(self, parameters=None):
        """
        Construct a new CancerSimulation.

        :param parameters: The cancer simulation parameters
        :type  parameters: CancerSimulationParameters

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
        self.__benefitial_mut = []
        self.__growth_plot_data = None
        self.__mutCounter = None
        self.__s = [self.parameters.mutations_per_division]*100000

    def run(self):
        # Setup square matrix.
        matrix_size=self.parameters.matrix_size

        #ctDNA_lifetime=self.parameters.ctDNA_lifetime

        #s = np.random.poisson(self.parameters.mutations_per_division, 100000)
        # introducing cancer cell
        initLoc=(int(matrix_size/2),int(matrix_size/2))
        print(initLoc)

        # value within matrix represents index of the mutation container
        # in this case number one point towards (0,1).
        self.__mtx[initLoc]=1
        self.__mut_container=[(0, 0), (0, 1)]

        # create lists used in loops
        lq_bipsy=[]
        self.__growth_plot_data=[]

        self.__pool=[initLoc]   #start the pool of cancer cells by adding the initial cancer cell into it

        print('Tumour growth in progress...')

        death_list=[]
        prop_of_driver=[]

        start = timer()

        seed=1
        prng.seed(seed)

        true_vaf=self.tumourGrowth()

        print(self.__mtx.todense())

        end=timer()
        print("time", end - start)


    def sampling(self, sample):
        """ TODO: Add a short documentation of this function.

        :param <+variable name+>: <+variable doc+>
        :type  <+variable name+>: <+type of variable+>

        """
        dna_from_sample=[self.__mutation_reconstruction(self.__mtx[i]) for i in sample]

        biopsy_raw_vaf=self.bulk_seq(dna_from_sample, self.parameters.number_of_generations, benefitial=False, sampling_or_fullTumour="Sample", )

        vaf=self.increase_mut_number(biopsy_raw_vaf)

        vaf=sorted(vaf, key=lambda x: x[0])

        return vaf

    def simulate_seq_depth(self, extended_vaf):
        """ TODO: Add a short documentation of this function.

        :param <+variable name+>: <+variable doc+>
        :type  <+variable name+>: <+type of variable+>

        """

        cellnum=extended_vaf[0]

        AF=np.array(extended_vaf)

        AF = AF/self.parameters.ploidy

        AF = AF * self.parameters.cellularity

        AF=np.extract(AF>self.parameters.detectionlimit*cellnum, AF)

        depth=np.random.poisson(self.parameters.read_depth, len(AF))

        samp_alleles=np.random.binomial(depth, AF/cellnum)

        VAF = samp_alleles/depth

        return VAF

    def increase_mut_number(self, solid_pre_vaf):
        """ TODO: Add a short documentation of this function.

        :param <+variable name+>: <+variable doc+>
        :type  <+variable name+>: <+type of variable+>

        """

        solid_extended_vaf=[]

        target_mut_solid=[]

        for i in solid_pre_vaf:

            if i[0]==1: #first mutation duplicate N number of times, adding additional clonal mutations

                for rep in range(self.parameters.number_of_clonal):

                    target_mut_solid.append((i[0] , float(i[1])))

            else:
                # for all subsequent mutations duplicate number
                # of them based on poisson distribution in variable s
                for rep in range(self.__s[i[0]]):
                    target_mut_solid.append((i[0], float(i[1])))

        return target_mut_solid

    def terminate_cell(self, cell, step):
        """ TODO: Add a short documentation of this function.

        :param <+variable name+>: <+variable doc+>
        :type  <+variable name+>: <+type of variable+>

        """

        self.__pool.remove(cell)
        death_list.append((self.__mtx[cell], step))
        self.__mtx[cell]=0

    def death_one_cell_chunk(self, step):
        """ TODO: Add a short documentation of this function.

        :param <+variable name+>: <+variable doc+>
        :type  <+variable name+>: <+type of variable+>

        """

        NtoDie=self.parameters.prop_of_expired_cells*len(self.__pool)
        NtoDie=math.ceil(NtoDie)
        toDie=prng.sample(self.__pool, math.ceil(NtoDie))

        for i in toDie:
            terminate_cell(i, self.__pool, step)

    def bulk_seq(self, DNA, step, benefitial, sampling_or_fullTumour):
        """ TODO: Add a short documentation of this function.

        :param <+variable name+>: <+variable doc+>
        :type  <+variable name+>: <+type of variable+>

        """

        #print(DNA)
        #print('dna', len(DNA))
        vaf_bulk=[]
        cellnum=[]

        #print('benefitial', benefitial)

        reduced=list(itertools.chain(*[j for j in DNA]))  #flatten the list of mutations

        # NOTE: np.hist?
        for i in set(reduced): #count number of unique mutations in whole tumour at time step
            vaf_bulk.append((i, float(reduced.count(i))))

        prop_of_driver=[]
        tum_size=len(DNA)
        # print('tum_size', tum_size)
        # print('vaf_bulk', vaf_bulk)

        if benefitial:

            for i in vaf_bulk:
                if i[0]==benefitial[0]:
                    prop_of_driver.append((i[0], i[1]/tum_size, step))

        return vaf_bulk

    def mutation_reconstruction(self, cellToUntangle):
        """ TODO: Add a short documentation of this function.

        :param <+variable name+>: <+variable doc+>
        :type  <+variable name+>: <+type of variable+>

        """

        mut_prof=[]

        m=self.__mut_container[cellToUntangle][0]

        mut_prof.append(self.__mut_container[cellToUntangle][1])

        for i in range(cellToUntangle, 0, -1):
            if self.__mut_container[i][1]==m:
                mut_prof.append(self.__mut_container[i][1])
                m=self.__mut_container[i][0]

        return mut_prof[::-1]

    def neighbours(self, rndNode):
        """ Returns the nearest-neighbor cells around the given node.

        :param rndNode: The node for which to calculate the neighbors.
        :type  rndNode: tuple (i,j) of cell indices.

        """
        # make list of all surrounding nodes
        neighboursList=[
                (rndNode[0]-1, rndNode[1]+1),
                (rndNode[0]  , rndNode[1]+1),
                (rndNode[0]+1, rndNode[1]+1),
                (rndNode[0]-1, rndNode[1]  ),
                (rndNode[0]+1, rndNode[1]  ),
                (rndNode[0]-1, rndNode[1]-1),
                (rndNode[0]  , rndNode[1]-1),
                (rndNode[0]+1, rndNode[1]-1)]

        # return nodes that are not cancerous, do not contain mutation index
        return [y for y in neighboursList if self.__mtx[y]==0]

    def place_to_divide(self):
        """ TODO: Add a short documentation of this function.

        :param <+variable name+>: <+variable doc+>
        :type  <+variable name+>: <+type of variable+>

        """

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

    def tumourGrowth(self):
        """ Run the tumour growth simulation.  """

        # setup a counter to keep track of number of mutations that occur in this run.
        mutCounter=1

        # setup a list of mutations with positive fitness effect. empty in the beginning.
        benefitial_mut=[]

        # Loop over time steps.
        for step in range(self.parameters.number_of_generations):
            print(self.__mtx.todense())
            print(' ')
            print('>>>>>>>>>>>>>>>>>>>>>>>>')
            print(' ')
            print('step in cancer growth', step)

            # bulk_vaf=self.bulk_seq([mutation_reconstruction(self.__mtx[i]) for i in pool], step, benefitial_mut)

            # setup a temporary list to store the mutated cells in this iteration.
            temp_pool=[]

            # reshuffle the order of pool to avoid that cells with low number divide always first.
            shuffle(self.__pool)

            print('list of cancer cells', self.__pool)
            #print('mut_container_', mut_container)

            # Loop over all cells in the pool.
            for cell in self.__pool:
                print(' ')
                print('________')

                print('cell to divide ', cell)
                print(' ')

                #renamed randomNode to fit with rest of the code, need to change that.
                rndNode=cell

                # Get the existing neighboring cells.
                neigh=self.neighbours(cell)

                # place_to_divide=place_to_divide()

                if neigh:
                # first condition: if available neighbors

                    # if cell has benefitial mutation.
                    if self.__mtx[cell] in benefitial_mut:

                        # cell divides with greater probability.
                        if prng.random()<self.parameters.fitness_advantageous_division_probability:
                            place_to_divide=prng.choice(neigh)
                            temp_pool.append(place_to_divide)

                            # daughter cells mutates
                            if prng.random()<self.parameters.mutation_rate:
                                # print('new mutation will occur', mutCounter+1)

                                mutCounter=mutCounter+1
                                # New cell gets the index number of largest number of mutation
                                self.__mtx[place_to_divide]=len(self.__mut_container)
                                self.__mut_container.append((self.__mut_container[self.__mtx[rndNode]][1], mutCounter))

                                # print((self.__mut_container[self.__mtx[rndNode]][1], mutCounter+1))
                                # add the index to the list of the benefitial ones
                                benefitial_mut.append(int(self.__mtx[place_to_divide]))

                                # mother cell mutates
                                mutCounter=mutCounter+1
                                self.__mut_container.append((self.__mut_container[self.__mtx[rndNode]][1], mutCounter))

                                self.__mtx[cell]=len(self.__mut_container)-1

                            else:
                            # if there is no new mutation just copy the index from mother cell
                                self.__mtx[place_to_divide]=self.__mtx[rndNode]

                    # second condition...if it is not on the list of cells with advantage, use normal division probability
                    else:
                        if prng.random()<self.parameters.division_probability:
                            place_to_divide=prng.choice(neigh)
                            print('index of the mother cell', self.__mtx[cell])
                            print('random neighbor to divide', place_to_divide)

                            # temp_pool will be updated to pool of cancer cells
                            # after all cells attempt to divide, this prevents
                            # that new cells divide in the same turn.
                            ### NOTE: But you're updating self.__mtx during the step.
                            temp_pool.append(place_to_divide)

                            # here is main thing, updating new cells and mutations
                            if prng.random()<self.parameters.mutation_rate:
                                # new cell gets the index number of largest number of mutation
                                mutCounter=mutCounter+1
                                self.__mtx[place_to_divide]=len(self.__mut_container)
                                print('neigh cell got new index', len(self.__mut_container))
                                print((self.__mut_container[self.__mtx[rndNode]][1], mutCounter))

                                self.__mut_container.append((self.__mut_container[self.__mtx[rndNode]][1], mutCounter))
                                print('mut container updated', self.__mut_container)

                                # NOTE: Why the second condition (len(benefitial_mut)==0) ???
                                if prng.random()<self.parameters.advantageous_mutation_probability and len(benefitial_mut)==0 and step==self.parameters.time_of_advantageous_mutation:
                                    print('new benefitial mutation!!!!!!', int(self.__mtx[place_to_divide]))
                                    benefitial_mut.append(int(self.__mtx[place_to_divide]))

                                # mother cell mutates
                                mutCounter=mutCounter+1
                                self.__mut_container.append((self.__mut_container[self.__mtx[rndNode]][1], mutCounter))

                                # print('mut container updated second time', self.__mut_container)

                                # update mutation list.
                                # print('mother cell gets new index', len(self.__mut_container)-1)
                                self.__mtx[cell]=len(self.__mut_container)-1

                            else:
                                print('no new mutation in normal division, inhereting from parent')
                                self.__mtx[place_to_divide]=self.__mtx[rndNode]
                                temp_pool.append(place_to_divide)

                    # temp_pool.append(place_to_divide)
                    # pool.remove(cell)

                # print('self.__mtx', self.__mtx.toarray())

            # add new cancer cells to a pool of cells available for division next round
            [self.__pool.append(v) for v in temp_pool]
            self.__growth_plot_data.append(len(self.__pool))


            # at the end reconstruct mutational frequencies from the whole tumour
            if step == self.parameters.number_of_generations-1:

                bulk_vaf=self.bulk_seq([self.mutation_reconstruction(self.__mtx[i]) for i in self.__pool], step, benefitial_mut, sampling_or_fullTumour="Full")

                bulk_vaf=self.increase_mut_number(bulk_vaf)
                print('1', bulk_vaf[0:10])

                return bulk_vaf

def main(arguments):
    """ TODO: Add a short documentation of this function.

    :param <+variable name+>: <+variable doc+>
    :type  <+variable name+>: <+type of variable+>

    """

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

#Avoid execution of main if the script is imported ad a module
if __name__ == "__main__":
    # Entry point
    # Setup the command line parser.
    parser = ArgumentParser()

    # Seed parameter.
    parser.add_argument("seed",
                        help="The prng seed.",
                        )

    # Parse the arguments.
    arguments = parser.parse_args()

    main(arguments)

