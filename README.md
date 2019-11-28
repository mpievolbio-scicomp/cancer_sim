CancerSim: A Cancer Simulation Package for python3
==============================================


Background
----------
Cancer is a group of complex diseases characterized by excessive cell proliferation, invasion and destruction of the surrounding tissue.
Its high division and mutation rates lead to excessive intratumour genetic heterogeneity which makes cancer highly adaptable to environmental pressures such as therapy.
Throughout most of its existence tumour is inaccessible to direct observation.
Some examples where computational models can be of great use include early carcinogenesis as lesions are clinically observable when they already contain millions of cells, seeding of metastases, cancer cell dormancy.
Here, we present a software that simulates somatic evolution of cancer and produces virtual spatial tumours with variable extent of intratumour genetic heterogeneity.

Tumour is simulated using two-dimensional, on-lattice agent-based model.
Our model is abstract, not specific to any neoplasm type and does not consider variety of biological features commonly found in neoplasm such as vasculature, immune contexture, availability of nutrients and architecture of the tumour surroundings.
It is, however, most suitable for simulating mostly superficially spreading tumours like carcinoma in situ, skin cancers, gastric cancers etc.

Tumour lattice structure is established by a sparse matrix where each position on a matrix corresponds to an individual cell.
Each cell is surrounded by eight neighbouring cells (Moore neighbourhood).
Value in the matrix is an index pointing to the last mutation cell acquired in the stored list of mutations.

Simulation runs in discrete time-steps. Each turn every cell in the tumour that has unoccupied neighbour can divide with certain probability (params.div__probability).
With each division daughter cell inherits all mutations from parent cell and acquires new mutation with a given probability (params.mut_rate).
Different division probabilities can be introduced for some cells in order to simulate variability in fitness of cells who acquired beneficial or deleterious mutation.

Throughout the cancer growth phase, for every cell, CancerSim stores information about parent cell and a designation of newly acquired mutation.
Complete mutational profiles of cells are reconstructed a posteriori based on the stored lineage information.

Rule that allows only cells with empty neighbouring nodes to divide leads to predominantly peripheral growth of tumour and complete absence of dynamics in the centre.
To make the system more realistic and allow for some degree of growth inside the tumour, we introduced a death process.
Every time step, after all cells attempt their division, a number of random cells expires leaving their position free to host a new cancer cell.

After the simulation, frequency of each mutataion in the tumour is quantified and exported to a file.
We introduces variable amount sequencing noise (cite williams) which makes output data more biologically realistic.
Further, tumour can be sampled and number, or frequency of mutations quantified and visualized as histogram of mutation frequencies.

Simulation is written in Python and can be imported as Anaconda package.

We implemented growth visualizer for two-dimensional tumour.

High--level functionality
-------------------------
The parameter values of the cancer simulation are given via a python module or
programmatically via the ```CancerSimulationParameters``` class. A documented
example is included in the released distribution and reproduced here:

    $> cat/params.py
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
    mut_rate                        = 0.8

    # Probability that a mutation is advantagous (??).
    advantageous_mut_prob           = 1

    # Number of mutations per cell division.
    mut_per_division                = 1

    # Time after which adv. mutations occur.
    time_of_adv_mut                 = 2

    # Factor to scale the simulation to actual cell count numbers.
    num_of_clonal                   = 1

    # Tumour multiplicity.
    tumour_multiplicity             = None


The simulation is started from the command line. The syntax is

    $> python -m casim.py [-h] [-o DIR] seed

 The mandatory command line argument ```seed``` is the
random seed. Using the same seed on two simulation runs with identical
parameters results in identical results, this may be used for testing and
debugging. The optional argument ```DIR``` specifies the directory where to
store the simulation log and output data. That directory must exist.

### Example 1

    $> python -m casim.casim 1

### Example 2

    $> mkdir sim_out
    $> python -m casim.casim.py -o sim_out

Results will be stored in the newly created directory ```sim_out/```.

Reference Manual
----------------
The API reference manual is available at [our pages site](https://c.fortmanngrote.pages.gwdg.de/cancer_sim).

Examples
--------
See our quickstart example in
`docs/source/include/notebooks/quickstart_example.ipynb`. Or [launch it in Binder](https://mybinder.org/v2/git/https%3A%2F%2Fgitlab.gwdg.de%2Fc.fortmanngrote%2Fcancer_sim/develop?filepath=https%3A%2F%2Fgitlab.gwdg.de%2Fc.fortmanngrote%2Fcancer_sim%2Fblob%2Fdevelop%2Fdocs%2Fsource%2Finclude%2Fnotebooks%2Fquickstart_example.ipynb)



