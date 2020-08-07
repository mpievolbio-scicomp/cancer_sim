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

# Number of generations after which adv. mutation occurs (>=0).
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


