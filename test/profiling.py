import cProfile

from casim.casim import CancerSimulatorParameters, CancerSimulator

import logging
logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.DEBUG)

parameters = CancerSimulatorParameters(
                                       matrix_size =                           100,
                                       number_of_generations =                 20  , # vary
                                       division_probability =                  1.0, # 1 => exp. growth
                                       advantageous_division_probability =      0.3,
                                       death_probability =                      0.0,
                                       fitness_advantage_death_probability =    0.0,
                                       mutation_rate =                          1.0, # 1 =>
                                       advantageous_mutation_probability =      0.8,
                                       mutations_per_division =                10  , # if mutation event occurs, have this number of mutation
                                       time_of_advantageous_mutation =      30000  , # large to not kick in
                                       number_of_clonal =                       2  , # initial number of mutations in first cancer cell
                                       tumour_multiplicity =                'single',
                                      )


simulator = CancerSimulator(parameters, seed=None, outdir="profiling_out_new")
commands = "simulator.run()"
cProfile.run(commands)

#seed  = 1
#age = "new"
#simulator = CancerSimulator(parameters, seed=1, outdir="profiling_out_%s" % age)
#simulator.run()

