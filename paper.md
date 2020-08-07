---
title: 'CancerSim: A Cancer Simulation Package for python3'
tags:
  - stochastic simulation
  - tumour growth
  - tumour sampling
  - cancer biology
authors:
  - name: Luka Opasic
    orcid: 0000-0001-7595-1722
    affiliation: 1
  - name: Jacob Scott
    orcid: 0000-0003-2971-7673
    affiliation: 2
  - name: Arne Traulsen
    orcid: 0000-0002-0669-5267
    affiliation: 1
  - name: Carsten Fortmann-Grote
    orcid: 0000-0002-2579-5546
    affiliation: 1
affiliations:
  - name: Max-Planck-Institute for Evolutionary Biology, Plön, Germany
    index: 1
  - name: Cleveland Clinic, Cleveland, OH, US
    index: 2
date: 17 March 2020
bibliography: references.bib
---


Summary
----------

Cancer is a group of complex diseases characterized by excessive cell
proliferation, invasion, and destruction of the surrounding tissue
[@Kumar2017]. Its high division and mutation rates
lead to excessive intratumour genetic heterogeneity which makes cancer
highly adaptable to environmental pressures such as therapy
[@Turajlic2019]. Throughout most of its existence
tumour is inaccessible to direct observation and experimental
evaluation. Therefore, computational modelling can be useful to study
many aspects of cancer. Some examples where theoretical models can be of
great use include early carcinogenesis, as lesions are clinically
observable when they already contain millions of cells, seeding of
metastases, and cancer cell dormancy
[@Altrock2015].

Here, we present CancerSim, a software that simulates somatic evolution of
tumours. The software produces virtual spatial tumours with variable extent of
intratumour genetic heterogeneity and realistic mutational profiles. Simulated
tumours can be subjected to multi-region sampling to obtain mutation profiles
that are realistic representation of the sequencing data. This makes the
software useful for studying various sampling strategies in clinical cancer
diagnostics. An early version of this cancer evolution model was used to
simulate tumours subjected to sampling for classification of mutations based on
their abundance [@Opasic2019]. Target users are scientists working in the field
of mathematical oncology and students with interest in studying somatic
evolution of cancer.  


Our model is abstract, not specific to any neoplasm type and does not
consider a variety of biological features commonly found in neoplasm
such as vasculature, immune contexture, availability of nutrients, and
architecture of the tumour surroundings. It resembles the most to
superficially spreading tumours like carcinoma in situ, skin cancers, or
gastric cancers, but it can be used to model any tumour on this abstract
level.

The tumour is simulated using a two-dimensional, on-lattice, agent-based
model. The tumour lattice structure is established by a sparse matrix
whose non-zero elements correspond to the individual cells. Each cell is
surrounded by eight neighbouring cells (Moore neighbourhood). The value
of the matrix element is an index pointing to the last mutation cell
acquired in the list of mutations which is updated in each simulation
step.

The simulation advances in discrete time-steps. In each simulation step,
every tumour cell in the tumour that has an unoccupied neighbour can
divide with a certain probability (params.div\_\_probability). The
daughter cell resulting from a cell division inherits all mutations from
the parent cell and acquires a new mutation with a given probability
(params.mut\_prob). Different division probabilities can be introduced
for some cells in order to simulate variability in fitness of cells that
acquired a beneficial or deleterious mutation. The simulation allows the
acquisition of more than one mutational event per cell
(params.mut\_per\_division). In that case, variable amounts of
sequencing noise [@Williams2016] can be added to make
the output data more biologically realistic.

Throughout the cancer growth phase, CancerSim stores information about
the parent cell and a designation of newly acquired mutations for every
cell. Complete mutational profiles of cells are reconstructed a
posteriori based on the stored lineage information.

The division rules which allow only cells with empty neighbouring nodes
to divide, cause exclusively peripheral growth and complete absence of
dynamics in the tumour centre. To allow for variable degree of growth
inside the tumour, we introduced a death process. At every time step,
after all cells attempt their division, a number of random cells die and
yield their position to host a new cancer cell in a subsequent time
step.

After the simulation, the tumour matrix, and the lists of lineages and
frequencies of each mutation in the tumour are exported to files.
Furthermore, the virtual tumour can be sampled and a histogram over the
frequency of mutations will be visualised. Alternatively, a saved tumour
can be loaded from file and then subjected to the sampling process.

Download and Installation
-------------------------

CancerSim is written in Python (version \>3.5). We recommend to install
it directly from the source code hosted at github <https://github.com/mpievolbio-scicomp/cancer_sim>.

 Detailed instructions including creation of a
`conda` environment are given in the online documentation at <https://cancer-sim.readthedocs.io/en/master/include/README.html#installation>.

Testing
-------

Although not strictly required, we recommend to run the test suite after
installation. Simply execute the `run_tests.sh` shell script:

    $> ./run_tests.sh

This will generate a test log named `casim_test@<timestamp>.log` with
`<timestamp>` being the date and time when the test was run.

The test suite is automatically run after each commit to the code base.
Results are published on
[travis-ci.org](https://travis-ci.org/mpievolbio-scicomp/cancer_sim).

High--level functionality
-------------------------

The parameters of the cancer simulation are given via a python module or
programmatically via the `CancerSimulationParameters` class. The file
`params.py` is a documented parameter module:

```    
################################################################################
#                                                                              #
# Commented casim parameter input file.                                        #
# Valid settings are indicated in parentheses at the end of each comment line. #
# [0,1] stands for the closed interval from 0 to 1, including the limits; ||   #
# means "or".                                                                  #
#                                                                              #
################################################################################

# Number of mesh points in each dimension (>0)
matrix_size = 1000

# Number of generations to simulate (>0).
number_of_generations = 20

# Probability of cell division per generation ([0,1]).
division_probability = 1

# Probability of division for cells with advantageous mutation ([0,1]).
adv_mutant_division_probability = 1

# Fraction of cells that die per generation ([0,1]).
death_probability = 0.1

# Fraction of cells with advantageous mutation that die per generation ([0,1]).
adv_mutant_death_probability = 0.0

# Probability of mutations ([0,1]).
mutation_probability = 1

# Mutation probability for the adv. cells ([0,1]).
adv_mutant_mutation_probability = 1

# Number of mutations per cell division (>=0).
number_of_mutations_per_division = 10

# Number of generation after which adv. mutation occurs (>=0).
adv_mutation_wait_time = 10

# Number of mutations present in first cancer cell (>=0).
number_of_initital_mutations = 150

# Tumour multiplicity ("single" || "double").
tumour_multiplicity = "double"

# Sequencing read depth (read length * number of reads / genome length).
read_depth = 100

# Fraction of cells to be sampled ([0,1]).
sampling_fraction = 0.9
    
# Plot the tumour growth curve (True || False).
plot_tumour_growth = True
    
# Export the tumour growth data to file (True || False).
export_tumour = True
```

The simulation is started from the command line. The syntax is

    $> python -m casim.casim [-h] [-o DIR] seed

The mandatory command line argument `seed` is the random seed. The optional
`DIR` specifies the output directory.

Reference Manual
----------------

The API reference manual is available at
<https://cancer-sim.readthedocs.io>.

Examples
--------

See our quickstart example in
`docs/source/include/notebooks/quickstart_example.ipynb` or use the following link to [launch it in Binder](https://mybinder.org/v2/gh/mpievolbio-scicomp/cancer_sim.git/master?filepath=docs%2Fsource%2Finclude%2Fnotebooks%2Fquickstart_example.ipynb).

References
----------

