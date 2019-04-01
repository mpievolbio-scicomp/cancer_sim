import cProfile

from casim.casim import main
from collections import namedtuple

commands = "arguments = namedtuple('arguments', ['seed']); arguments.seed=1; main(arguments=arguments)"

cProfile.run(commands)
