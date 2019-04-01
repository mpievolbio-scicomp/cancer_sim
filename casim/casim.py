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

def sampling(sample):
    """ TODO: Add a short documentation of this function.

    :param <+variable name+>: <+variable doc+>
    :type  <+variable name+>: <+type of variable+>

    """
    dna_from_sample=[mutation_reconstruction(mtx[i]) for i in sample]

    biopsy_raw_vaf=bulk_seq(dna_from_sample, params.num_of_generations, benefitial=False, sampling_or_fullTumour="Sample", )

    vaf=increase_mut_number(biopsy_raw_vaf)

    vaf=sorted(vaf, key=lambda x: x[0])

    return vaf

def simulate_seq_depth(extended_vaf):
    """ TODO: Add a short documentation of this function.

    :param <+variable name+>: <+variable doc+>
    :type  <+variable name+>: <+type of variable+>

    """

    cellnum=extended_vaf[0]

    AF=np.array(extended_vaf)

    AF = AF/params.ploidy

    AF = AF * params.cellularity

    AF=np.extract(AF>params.detectionlimit*cellnum, AF)

    depth=np.random.poisson(params.read_depth, len(AF))

    samp_alleles=np.random.binomial(depth, AF/cellnum)

    VAF = samp_alleles/depth

    return VAF

def increase_mut_number(solid_pre_vaf):
    """ TODO: Add a short documentation of this function.

    :param <+variable name+>: <+variable doc+>
    :type  <+variable name+>: <+type of variable+>

    """

    solid_extended_vaf=[]

    target_mut_solid=[]

    for i in solid_pre_vaf:

        if i[0]==1: #first mutation duplicate N number of times, adding additional clonal mutations

            for rep in range(params.num_of_clonal):

                target_mut_solid.append((i[0] , float(i[1])))

        else:
            # for all subsequent mutations duplicate number
            # of them based on poisson distribution in variable s
            for rep in range(s[i[0]]):
                target_mut_solid.append((i[0], float(i[1])))

    return target_mut_solid

def terminate_cell(cell, pool, step):
    """ TODO: Add a short documentation of this function.

    :param <+variable name+>: <+variable doc+>
    :type  <+variable name+>: <+type of variable+>

    """

    pool.remove(cell)
    death_list.append((mtx[cell], step))
    mtx[cell]=0

def death_one_cell_chunk(pool, step):
    """ TODO: Add a short documentation of this function.

    :param <+variable name+>: <+variable doc+>
    :type  <+variable name+>: <+type of variable+>

    """

    NtoDie=params.prop_of_expired_cells*len(pool)
    NtoDie=math.ceil(NtoDie)
    toDie=prng.sample(pool, math.ceil(NtoDie))

    for i in toDie:
        terminate_cell(i, pool, step)

def bulk_seq(DNA, step, benefitial, sampling_or_fullTumour):
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

def mutation_reconstruction(cellToUntangle):
    """ TODO: Add a short documentation of this function.

    :param <+variable name+>: <+variable doc+>
    :type  <+variable name+>: <+type of variable+>

    """

    mut_prof=[]

    m=mut_container[cellToUntangle][0]

    mut_prof.append(mut_container[cellToUntangle][1])

    for i in range(cellToUntangle, 0, -1):
        if mut_container[i][1]==m:
            mut_prof.append(mut_container[i][1])
            m=mut_container[i][0]

    return mut_prof[::-1]

def neighbours(rndNode):
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
    return [y for y in neighboursList if mtx[y]==0]

def place_to_divide():
    """ TODO: Add a short documentation of this function.

    :param <+variable name+>: <+variable doc+>
    :type  <+variable name+>: <+type of variable+>

    """

    a = prng.randint(0,params.matrixSize-1)
    b = prng.randint(0,params.matrixSize-1)
    random_place_to_divide=(a,b)

    if mtx[random_place_to_divide]==0:
        return a, b

    else:
        while mtx[random_place_to_divide]!=0:
            a = prng.randint(0,params.matrixSize-1)
            b = prng.randint(0,params.matrixSize-1)
            random_place_to_divide=(a,b)

        return a, b

def tumourGrowth(pool, n_time_steps, div_probability, death_probability):
    """ TODO: Add a short documentation of this function.

    :param pool: The pool of cancer cells.
    :type  pool: list

    :param n_time_steps: The number of time steps to simulate.
    :type  n_time_steps: int

    :param div_probability: The probability of cell division in each timestep.
    :type  div_probability: float (0 <= div_probability <= 1).

    :param death_probability: The cell death probability per timestep.
    :type  death_probability: float (0 <= death_probability <= 1).

    """

    # setup a counter to keep track of number of mutations that occur in this run.
    mutCounter=1

    # setup a list of mutations with positive fitness effect. empty in the beginning.
    benefitial_mut=[]

    # Loop over time steps.
    for step in range(n_time_steps):
        print(mtx.todense())
        print(' ')
        print('>>>>>>>>>>>>>>>>>>>>>>>>')
        print(' ')
        print('step in cancer growth', step)

        # bulk_vaf=bulk_seq([mutation_reconstruction(mtx[i]) for i in pool], step, benefitial_mut)

        # setup a temporary list to store the mutated cells in this iteration.
        temp_pool=[]

        # reshuffle the order of pool to avoid that cells with low number divide always first.
        shuffle(pool)

        print('list of cancer cells', pool)
        #print('mut_container_', mut_container)

        # Loop over all cells in the pool.
        for cell in pool:
            print(' ')
            print('________')

            print('cell to divide ', cell)
            print(' ')

            #renamed randomNode to fit with rest of the code, need to change that.
            rndNode=cell

            # Get the existing neighboring cells.
            neigh=neighbours(cell)

            # placeToDivide=place_to_divide()

            if neigh:
            # first condition: if available neighbors

                # if cell has benefitial mutation.
                if mtx[cell] in benefitial_mut:

                    # cell divides with greater probability.
                    if prng.random()<params.fittnes_advantage_div_prob:
                        placeToDivide=prng.choice(neigh)
                        temp_pool.append(placeToDivide)

                        # daughter cells mutates
                        if prng.random()<mut_rate:
                            # print('new mutation will occur', mutCounter+1)

                            mutCounter=mutCounter+1
                            # New cell gets the index number of largest number of mutation
                            mtx[placeToDivide]=len(mut_container)
                            mut_container.append((mut_container[mtx[rndNode]][1], mutCounter))

                            # print((mut_container[mtx[rndNode]][1], mutCounter+1))
                            # add the index to the list of the benefitial ones
                            benefitial_mut.append(int(mtx[placeToDivide]))

                            # mother cell mutates
                            mutCounter=mutCounter+1
                            mut_container.append((mut_container[mtx[rndNode]][1], mutCounter))

                            mtx[cell]=len(mut_container)-1

                        else:
                        # if there is no new mutation just copy the index from mother cell
                            mtx[placeToDivide]=mtx[rndNode]

                # second condition...if it is not on the list of cells with advantage, use normal division probability
                else:
                    if prng.random()<div_probability:
                        placeToDivide=prng.choice(neigh)
                        print('index of the mother cell', mtx[cell])
                        print('random neighbor to divide', placeToDivide)

                        # temp_pool will be updated to pool of cancer cells
                        # after all cells attempt to divide, this prevents
                        # that new cells divide in the same turn.
                        ### NOTE: But you're updating mtx during the step.
                        temp_pool.append(placeToDivide)

                        # here is main thing, updating new cells and mutations
                        if prng.random()<mut_rate:
                            # new cell gets the index number of largest number of mutation
                            mutCounter=mutCounter+1
                            mtx[placeToDivide]=len(mut_container)
                            print('neigh cell got new index', len(mut_container))
                            print((mut_container[mtx[rndNode]][1], mutCounter))

                            mut_container.append((mut_container[mtx[rndNode]][1], mutCounter))
                            print('mut container updated', mut_container)

                            # NOTE: Why the second condition (len(benefitial_mut)==0) ???
                            if prng.random()<params.advantageous_mut_prob and len(benefitial_mut)==0 and step==params.time_of_adv_mut:
                                print('new benefitial mutation!!!!!!', int(mtx[placeToDivide]))
                                benefitial_mut.append(int(mtx[placeToDivide]))

                            # mother cell mutates
                            mutCounter=mutCounter+1
                            mut_container.append((mut_container[mtx[rndNode]][1], mutCounter))

                            # print('mut container updated second time', mut_container)

                            # update mutation list.
                            # print('mother cell gets new index', len(mut_container)-1)
                            mtx[cell]=len(mut_container)-1

                        else:
                            print('no new mutation in normal division, inhereting from parent')
                            mtx[placeToDivide]=mtx[rndNode]
                            temp_pool.append(placeToDivide)

                # temp_pool.append(placeToDivide)
                # pool.remove(cell)

            # print('mtx', mtx.toarray())

        # add new cancer cells to a pool of cells available for division next round
        [pool.append(v) for v in temp_pool]
        growth_plot_data.append(len(pool))


        # at the end reconstruct mutational frequencies from the whole tumour
        if step == params.num_of_generations-1:

            bulk_vaf=bulk_seq([mutation_reconstruction(mtx[i]) for i in pool], step, benefitial_mut, sampling_or_fullTumour="Full")

            bulk_vaf=increase_mut_number(bulk_vaf)
            print('1', bulk_vaf[0:10])
            return bulk_vaf

def main(arguments):
    """ TODO: Add a short documentation of this function.

    :param <+variable name+>: <+variable doc+>
    :type  <+variable name+>: <+type of variable+>

    """

    global prop_of_driver
    global mtx
    global mut_container
    global mut_rate
    global lq_bipsy
    global xaxis_histogram
    global death_list
    global biopsy_timing
    global benefitial_mut
    global growth_plot_data
    global mutCounter
    global s

    start = timer()

    seed=arguments.seed
    prng.seed(seed)

    # Load single parameter values.
    # matrix size (nrows = ncols)
    matrixSize=params.matrixSize  # value presents the size of both x and y axis
    mut_rate=params.mut_rate
    num_of_generations=params.num_of_generations
    #ctDNA_lifetime=params.ctDNA_lifetime

    benefitial_mut=[]
    div_probability=params.div_probability
    death_probability=params.death_probability

    #s = np.random.poisson(params.mut_per_division, 100000)
    s = [params.mut_per_division]*100000

    # initiate tumour as lil sparce matrix with integer values
    mtx=lil_matrix((matrixSize, matrixSize), dtype=int)

    # introducing cancer cell
    initLoc=(int(matrixSize/2),int(matrixSize/2))
    print(initLoc)

    # value within matrix represents index of the mutation container
    # in this case number one point towards (0,1).
    mtx[initLoc]=1
    mut_container=[(0, 0), (0, 1)]

    # create lists used in loops
    lq_bipsy=[]
    growth_plot_data=[]

    pool=[initLoc]   #start the pool of cancer cells by adding the initial cancer cell into it

    print('Tumour growth in progress...')

    death_list=[]
    prop_of_driver=[]

    ##############MAIN FUNCTION
    true_vaf=tumourGrowth(pool, num_of_generations, div_probability, death_probability)
    ##################

    #print(mtx)
    print(mtx.todense())

    end=timer()
    print("time", end - start)

#Avoid execution of main if the script is imported ad a module
if __name__ == "__main__":
    # Entry point
    # Setup the command line parser.
    parser = ArgumentParser()

    # Parameters argument.
    #parser.add_argument("parameters",
                        #default='params',
                        #help="Python module that contains the parameters for this run."
                        #)

    # Seed parameter.
    parser.add_argument("seed",
                        help="The prng seed.",
                        )

    # Parse the arguments.
    arguments = parser.parse_args()

    main(arguments)

