import cProfile

from casim.casim import CancerSimulatorParameters, CancerSimulator

parameters = CancerSimulatorParameters(
                                       matrix_size =                           20  ,
                                       number_of_generations =                 10  , # vary
                                       division_probability =                   0.5, # 1 => exp. growth
                                       advantageous_division_probability =      0.3,
                                       death_probability =                      0.0,
                                       fitness_advantage_death_probability =    0.0,
                                       mutation_rate =                          0.2, # 1 =>
                                       advantageous_mutation_probability =      0.8,
                                       mutations_per_division =                10  , # if mutation event occurs, have this number of mutation
                                       time_of_advantageous_mutation =      30000  , # large to not kick in
                                       number_of_clonal =                       2  , # initial number of mutations in first cancer cell
                                       tumour_multiplicity =                'single',
                                      )


cProfile.run(commands)
