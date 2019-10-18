CaSim: A Cancer Simulation Package for python3
==============================================

Background
----------
Cancer is a group of complex diseases characterized by excessive cell proliferation, invasion and destruction of the surrounding tissue.
Its high division and mutation rates lead to excessive intratumour heterogeneity which makes cancer highly adaptable to environmental pressures such as therapy.
We have built a cancer model that produces tumours with variable levels of intratumour heterogeneity.
It is able to simulate cancer in multiple spatial dimensions:

(i) Well-mixed cancer which corresponds to haematological neoplasms,
(ii) Two dimensional on-lattice simulations suitable for modeling superficially spreading tumours like carcinoma in situ and finally
(iii) three dimensional on-lattice simulations for spatial solid tumours.

It contains two specific division mechanism both ran in discrete time-steps. In the first scenario each time-step only one cancer cell is chosen for division and in the second scenario every cancer cell has a certain probability to divide during one time-step.
In former case tumour growth is linear, in latter it is exponential if every cell divides each time-step.
Different division probabilities can be introduced for some cells in order to simulate variability in fitness of cells who acquired beneficial or deleterious mutation.

During each division mother cell can give birth to one daughter cell and remain unaltered.
It can also give birth to two daughter cells, one placed in available neighbour location and one daughter cell which replaces the mother cell at the same location.


Our model is abstract and does not consider variety of biological features commonly found in neoplasm such as vasculature, immune contexture, availability of nutrients and architecture of the tumour surroundings.


Simulated tumours can be pickled via dill package and further subjected to virtual biopsy sampling with frequencies of mutations present within the each sample as an output.
To make output results more biologically sound, recovered frequencies can be passed through function designed to simulate variable sequencing coverage depth and to introduce sequencing noise into the data.

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
The API reference manual is available at https://cancer-sim.gitlab.io.

Examples
--------
See our collection of jupyter notebooks.



