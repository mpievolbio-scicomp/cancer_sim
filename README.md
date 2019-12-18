[![Build Status](https://travis-ci.org/mpievolbio-scicomp/cancer_sim.svg?branch=master)](https://travis-ci.org/mpievolbio-scicomp/cancer_sim)
[![Documentation Status](https://readthedocs.org/projects/cancer-sim/badge/?version=latest)](https://cancer-sim.readthedocs.io/en/latest/?badge=latest)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mpievolbio-scicomp/cancer_sim/master?filepath=docs%2Fsource%2Finclude%2Fnotebooks%2Fquickstart_example.ipynb))
Background
----------

Cancer is a group of complex diseases characterized by excessive cell proliferation, invasion, and destruction of the surrounding tissue [@kumar:book:2017].  Its high division and mutation rates lead to excessive intratumour genetic heterogeneity which makes cancer highly adaptable to environmental pressures such as therapy [@turajlic:NRG:2019]. Throughout most of its existence tumour is inaccessible to direct observation and experimental evaluation.  Therefore, computational modelling can be useful to study many aspects of cancer. Some examples where theoretical models can be of great use include early carcinogenesis as lesions are clinically observable when they already contain millions of cells, seeding of metastases, and cancer cell dormancy [@altrock:NatRevCancer:2015].


Here, we present CancerSim, a software that simulates somatic evolution of tumours. The software produces virtual spatial tumours with variable extent of intratumour genetic heterogeneity and realistic mutational profiles.  Simulated tumours can be subjected to multi-region sampling to obtain mutation profiles that are realistic representation of the sequencing data. This makes the software useful for studying various sampling strategies in clinical cancer diagnostics.  Early version of this cancer evolution model was used to simulate tumours subjected to sampling for classification of mutations based on their abundance [@opasic:BMCCancer:2019].  Target users are scientists working in the field of mathematical oncology and students with interest in studying somatic evolution of cancer.

Our model is abstract, not specific to any neoplasm type and does not consider a variety of biological features commonly found in neoplasm such as vasculature, immune contexture, availability of nutrients, and architecture of the tumour surroundings.  It resembles the most to superficially spreading tumours like carcinoma in situ, skin cancers, or gastric cancers but it can be used to model any tumour.

The tumour is simulated using a two-dimensional, on-lattice, agent-based model.  The tumour lattice structure is established by a sparse matrix whose non-zero elements correspond to an individual cell.  Each cell is surrounded by eight neighbouring cells (Moore neighbourhood).  The value of the matrix element is an index pointing to the last mutation cell acquired in the list of mutations which is updated in each simulation step.

The simulation advances in discrete time-steps. In each simulation step, every tumour cell in the tumour that has an unoccupied neighbour can divide with a certain probability (params.div__probability).  The daughter cell resulting from a cell division inherits all mutations from the parent cell and acquires a new mutation with a given probability (params.mut_prob).  Different division probabilities can be introduced for some cells in order to simulate variability in fitness of cells that acquired beneficial or deleterious mutation.  The simulation allows the acquisition of more than one mutational event per cell (params.mut_per_division). In that case variable amount of sequencing noise [@williams:NG:2016] will be added to make the output data more biologically realistic.

Throughout the cancer growth phase,  CancerSim stores information about the parent cell and a designation of newly acquired mutations for every cell, Complete mutational profiles of cells are reconstructed a posteriori based on the stored lineage information.

The division rules which allow only cells with empty neighbouring nodes to divide, cause exclusively peripheral growth and complete absence of dynamics in the tumour centre.  To make the system more realistic and allow for variable degree of growth inside the tumour, we introduced a death process.  At every time step, after all cells attempt their division, a number of random cells die and yield their position to host a new cancer cell in a subsequent time step.

After the simulation, the tumour matrix, and the  lists of lineages and frequencies of each mutation in the tumour are exported to files.
Furthermore, the virtual tumour can be sampled and a histogram over the
 frequency of mutations will be visualised. Alternatively, a saved tumour can be loaded from file and then subjected to the sampling process.


## Installation
CancerSim is written in Python (version >3.5). It can be installed in various
ways, see below. Software dependencies are listed in the `requirements.txt`
file.


### PIP
```
$> pip install casim
```

### Conda
```
$> conda install -c conda-forge casim
```

### From the source code repository:
```
$> pip install .
```

We implemented growth visualizer for two-dimensional tumour **Where is it?**

High--level functionality
-------------------------
The parameters of the cancer simulation are given via a python module or
programmatically via the ```CancerSimulationParameters``` class. A documented
example `params.py` is included in the source code (under `test/params.py`) and reproduced here:

```
$> cat test/params.py
# Number of mesh points in each dimension
matrix_size                      = 100

# Number of generations to simulate.
num_of_generations              = 20

# Number of divisions per generation.
div_probability                 = 1

# Number of division for cells with mutation.
fittnes_advantage_div_prob      = 1

# Fraction of cells that die per generation.
dying_fraction                   = 0.1

# Fraction of cells with mutation that die per generation.
fitness_advantage_death_prob    = 0.0

# Rate of mutations.
mut_prob                        = 1

# Mutation probability for the adv. cells.
advantageous_mut_prob           = 1

# Number of mutations per cell division.
mut_per_division                = 10

# Time after which adv. mutations occur.
time_of_adv_mut                 = 2

# Number of mutations present in first cancer cell.
num_of_clonal                   = 15

# Tumour multiplicity.
tumour_multiplicity             = None

# Read depth.
read_depth                      = 100

# Fraction of cells to be sampled.
# sampling_fraction             = 0.1
```

The simulation is started from the command line. The syntax is

    $> python -m casim.py [-h] [-o DIR] seed

 The mandatory command line argument ```seed``` is the
random seed. Using the same seed on two simulation runs with identical
parameters results in identical results, this may be used for testing and
debugging. The optional argument ```DIR``` specifies the directory where to
store the simulation log and output data. If not given, output will be stored
in the directory `casim_out` in the current directory. For each seed, a subdirectory
`cancer_SEED` will be created.  If that subdirectory already exists because an
earlier run used the same seed, the run will abort. This is a safety catch to avoid overwriting data from previous runs.

### Example 1

    $> python -m casim.casim 1

### Example 2

    $> mkdir sim_out
    $> python -m casim.casim.py -o sim_out

Results will be stored in the newly created directory ```sim_out/```.

Reference Manual
----------------
The API reference manual is available at [https://cancer-sim.readthedocs.io](https://cancer-sim.readthedocs.io).

Examples
--------
See our quickstart example in
`docs/source/include/notebooks/quickstart_example.ipynb`. Or [launch it in Binder]([![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mpievolbio-scicomp/cancer_sim/master?filepath=docs%2Fsource%2Finclude%2Fnotebooks%2Fquickstart_example.ipynb))

