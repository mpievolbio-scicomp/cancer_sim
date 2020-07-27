#import cProfile
from time import sleep, time
import os
import sys
from io import StringIO

from casim.casim import CancerSimulatorParameters, CancerSimulator, HANDLER, LOGGER
import logging
level = logging.WARNING
HANDLER.setLevel(level)
LOGGER.setLevel(level)

##generations = [10,15,23,45,67,100,150,230,450,670,1000]
generations = [10,15,23,45,67,100,150, 230, 450, 670, 1000]
divprobs = [1.0]
#divprobs = [0.33]

n_run = sys.argv[1]
for dp in divprobs:
    for g in generations:
        parameters = CancerSimulatorParameters(
            matrix_size =                          4096,
            number_of_generations =                 g  , # vary
            division_probability =                    dp, # 1 => exp. growth
            adv_mutant_division_probability =       0.3,
            death_probability =                       0.0,
            adv_mutant_death_probabilityability =     0.0,
            mutation_rate =                           1.0, # 1 =>
            adv_mutant_mutation_probability =       0.8,
            mutations_per_division =                 10  , # if mutation event occurs, have this number of mutation
            adv_mutation_interval =       30000  , # large to not kick in
            number_of_initital_mutations =                        2  , # initial number of mutations in first cancer cell
            tumour_multiplicity =                'single',
                                              )

        simulator = CancerSimulator(parameters, seed=None, outdir="%s" % (n_run))

        t0 = time()
        simulator.run()
        T = time() - t0

        with open(os.path.join(simulator.outdir,"log.txt"), 'a') as fp:
            fp.write("%s\t%f\t%d\t%e\n" % (simulator.seed,dp,g,T))

        simulator.dump()
        sleep(1)
