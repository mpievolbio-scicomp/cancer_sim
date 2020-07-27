# Number of mesh points in each dimension
matrix_size = 1000

# Number of generations to simulate.
number_of_generations = 20

# Probability of cell division per generation.
division_probability = 1

# Probability for cell division for cells with adv. mutation.
adv_mutant_division_probability = 1

# Fraction of cells that die per generation.
death_probability = 0.1

# Fraction of cells with advantageous mutation that die per generation.
adv_mutant_death_probability = 0.0

# Probability that a single mutation event will occur at the time of division
mutation_probability = 1

# Probability that advantageous mutation will occur in time step params.adv_mutation_interval
adv_mutant_mutation_probability = 1

# Time after which advantageous mutations occur.
adv_mutation_interval = 10

# Number of mutations present in the first cancer cell
number_of_initital_mutations = 150

# Number of mutations per cell division.
number_of_mutations_per_division = 50

# Tumour multiplicity.
tumour_multiplicity = None

# Read depth
read_depth = 100

# Fraction of cells to sample (from interval [0,1))
sampling_fraction = 0.9

# Plot the growth curve.
plot_tumour_growth = True

# Write the tumour data to file.
export_tumour = False

