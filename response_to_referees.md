
# Reply to software paper review issue #16
> This issue is part of the JOSS review at [openjournals/joss-reviews#2436](https://github.com/openjournals/joss-reviews/issues/2436)
> 
> Please note that this represents my initial review and does not take into account recent efforts by the submitting author to address issues.
> # Software paper review
> 
> The work of Opasic et al. introduces the CancerSim software, which is a Python package used to simulate the somatic evolution of cancers. CancerSim employs a 2D lattice model, which allows the generational mutational profile and spatial heterogeneity of tumors to be modeled in an abstract manner as governed by key parameters related to cell division and mutation probabilities.
> ## Major Comments
> 
>     * [ ]  A detailed description of the simulation outputs is missing in both the paper and the software documentation (Documentation: Functional documentation). This may make it more difficult for users to fully evaluate if the software fits their needs (before using it), but perhaps more importantly, it is not initially clear how users can proceed with further analysis of the simulation outputs.
We have added documentation of the simulation outputs to the paper, the README
and the API documentation (`CancerSimulator.run()` method). Related commits:
3760a0, 4fb51c, 1e1ee0, 61220b .
> 
>     * [ ]   Although short descriptions of the key tumor models are provided in the paper¸ a detailed description of the key tumor model parameters and how they relate to each other is missing from both the paper and documentation (Documentation: Functional documentation).
All parameters are now documented in the paper, in the README, the online
documentation, and in the template
`params.py` file (as comments). In addition, the API reference manual for the
`CancerSimulationParameters` class was amended. Relevant commits: d6de50,
f95e93, d6a38d, 4a39b4d, f046b6, c96f26, 
> 
>     * [ ]  Although the summary of CancerSim itself is fairly detailed, there is no discussion of CancerSim in a broader context of cancer modeling methods or how CancerSim compares to other simulators. As such, I didn't get any sense of how CancerSim fits into the broader landscape of tumor modeling and cancer simulation methodologies. (Software paper: State of the field)
We added a short discussion of a recently published cancer simulation code
[waclaw2015] to the introduction. In the interest of overall paper length, we kept
this discussion to a minimum. The main point is that our model abstracts away
many of the fine grained model parameters of advanced tumour growth models. It
thereby offers an entry point for newcomers to the field to perform their first
steps in cancer modeling and to study the effect of choice of sampling positions
and sample size on the observed tumour profile in comparison to the whole tumour
profile. Relevant commits: 4a39b4
> 
> 
> ## Minor Comments
> 
>     * [ ]  This comment is related to the second major comment. There is no discussion or guidance on how users should determine key model parameters in a practical setting; e.g., should parameters be extracted directly from experimental measurements of cell division, or extracted from experimentally determined sequencing data (summary mentions reproducing sequencing data), or directly calibrated to some experimental data. I think some discussion and guidance in this area is needed to help users to evaluate the practical applicability of CancerSim to modeling their specific tumor systems, but will also be helpful for students and other non-specialists in further connecting the model to the tumor biology and available experimental methods.
A short recommendation for users how to apply their own parameters has been added to
the paper and the README. Further clarification about the role of each parameter
and how to choose its value and permissive values is given by the extended documentation of the
template parameter file `params.py` and the table explaining all parameters in the paper. 
Relevant commits: df2fade
> 
>     * [ ]  I think you may need to clarify that CancerSim installs as the `casim` package, both in the paper and the documentation.
A short note in this regard was added to the README, the online documentation,
and the paper. Relevant commits: e4161ff
> 
>     * [ ]  Summary, second paragraph, third sentence: This sentence is a little unclear to me. Does it mean that CancerSim can reproduce multi-region sequencing/profiling data? If so, I'd recommend re-wording a bit to make this more clear.
We clarified that sentence (4581588).
Indeed, our code allows to generate mutation
profiles from multiple samples taken at different positions of the tumour. This
feature was previously hidden from the user, i.e. the user would have to
hardcode the coordinates into the code `casim.py`. We have now exposed the
sampling coordinates in the public API through the parameter
`sampling_positions` in `CancerSimulationParameters`.
Relevant commits: 839c23a3, 7e32183
> 
>     * [ ]  Summary, second paragraph, fourth sentence: I'd recommend listing a few examples to help make the point more clearly.
We added needle biopsy and liquid biopsy as examples.
> 
>     * [ ]   Summary, first paragraph, third sentence: please change "existence tumour" to "existence a tumour".
done
> 
>     * [ ]  Summary, third paragraph, first sentence: I think a comma after "neoplasm type" is a little easier to read.
done.
> 
>     * [ ]  Summary, third paragraph, second sentence: I'd suggest rewriting "It resembles the most to superficially" as "It most closely resembles superficially"
done
> 
>     * [ ]  Summary, fourth paragraph, sentence three: Please replace "mutation cell" with "mutation the cell".
done
> 
>     * [ ]  Summary, paragraph five, sentence three: Can the division probabilities change during the course of a simulation in response to a cell acquiring a specific mutation? Or are the division probabilities just fixed at the beginning in a way that models beneficial or deleterious effects?
Both is possible: In a single run, the parameters `division_probability`,
`death_probability`, `mutation_probability` and their counterparts for carriers
of beneficial/deleterious mutations, 
`adv_mutant_division_probability`, `adv_mutant_mutation_probability`, and
`adv_mutant_death_probability` determine the growth rate of tumour cells
without and with beneficial/deleterious mutations.
These parameters are constant over the course of a run. It is however
possible to save a run with all its internal state variables to disk, reload the
run and continue the reloaded run after modifying selected parameters. This
opens up the possibility to modify growth and mutation parameters to model some
of the mentioned effects like dormancy or cancer therapy. This feature is now
documented in more detail in the paper, README, online documentation and in a
new example notebook. 
Relevant commit: 93e47a04



# Reply to software paper review issue #11


> This issue is part of the JOSS review at [openjournals/joss-reviews#2436](https://github.com/openjournals/joss-reviews/issues/2436)
> ## Software paper review
> 
> The work of Opasic and colleagues introduces CancerSim, a python package to study the evolution of somatic mutations in cancer cells, by introducing a simple simulation procedure on a 2D lattice and along discrete time steps. The package features key parameters, like cell division and mutation probabilities, tracks parent-daugther cell relationships and allows for adjustment of cell fitness of mutated cells. It therefore resembles an abstract model of superficially spreading tumors.
> 
> Although the paper consists of a well-written summary and an in-depth algorithm outline, description of key parameters remains superficial and applicability of the package in theoretical modelling is only vaguely delineated.
> ### Major comments
> 
>     * [ ]  The second paragraph of the **Summary** dives straight into specialist concepts like "cancer evolution", "somatic mutations", "genetic heterogeneity" and "sequencing". To enhance understanding for a diverse, non-specialist audience, please enhance the short conceptual introduction in the first paragraph, picking up those concepts. (Software paper: Summary)
We explained or replaced the technical/scientific terms by more general language
to enhance the accessibility.
Relevant commits: 9e2be3e
> 
>     * [ ]  There exists some discrepancy between parameter description and parameter naming (e.g. the number of divisions per generation vs. `div_probability`, number of division for cells with mutation vs. `fittnes_advantage_div_prob`). Please make sure to revise parameter names to indicate e.g. if a probability or absolute number is required, and include descriptions of which and how parameters interplay or depend on each other. (Documentation: Functional documentation)
All parameter names have been made consistent across the source code, the API
reference manual, README, paper, and online documentation. We explain the role
of all parameters in the template parameter module `params.py`. This file is
reproduced in the README, followed by additional explanations of each parameter
and their interplay. The `CancerSimulationParameters` API reference manual also
explains all parameters. A table with all parameter names, their function and
permitted values is added to the paper, replacing the code listing of the
`params.py` template parameter module. Relevant commits:
d6de50, b45278, 4a39b4d.
> 
>     * [ ]  A description of simulation output is missing in the paper and the online help. Having a detailed description of the output files lets readers quickly judge wheter or not a software will be suitable for their needs. (Documentation: Functional documentation)
We have added documentation of the simulation outputs to the paper, the README
and the API documentation (`CancerSimulator.run()` method). Also the quickstart
example now explains the output and renders the produced figures (mutation
profiles and growth curve.) Related commits:
3760a0, 4fb51c, 1e1ee0, 61220b .
> 
>     * [ ]  Having a detailed description of simulation output will also help readers to envision how the software can be applied in theoretical modelling and, importantly, which related research questions can be addressed e.g. in follow-up studies. Please delineate the softwares potential in much more detail, by suggestions or elaboration on possible downstream analysis. (Software paper: Statement of need). 
The word limit in JOSS allows us only to give a short discussion on potential
applications of our software. This has been added to the paper (Summary, last
paragraph.
Relevant commit: 5187dee, 5753fb5
> 
> 
> ### Minor comments
> 
>     * [ ]  The first paragraph mentions cancer cell dormancy as an example for beneficial use of theoretical models to understand carcinogenesis. Unfortunately, CancerSim does not allow cells to be simulated dormant, which renders the example detrimental to the papers statement.
Simulation of dormant cancer cells requires changing the growth parameters, in
particular `adv_mutant_division_probability` and `adv_mutant_death_probability`
during the course of a simulation. This can be achieved by snapshooting the
simulation after a certain number of generations corresponding to the dormancy
period, reloading the snapshot into memory, change the relevant parameters and
continue the simulation. While in principle possible in the unrevised
version, we have now documented this "dump-load-modify-restart" workflow
properly in a notebook which is referenced in the README and mentioned in the
paper. Some bugs had to be fixed in order to fully exploit this capability.
Relevant commits: 93e47a04 
> 
>     * [ ]  State of the field: Please include a short statement of which other tools are available for somatic mutation simulation and briefly compare features/drawbacks.
We added a short discussion of a recently published cancer simulation code
[waclaw2015] to the introduction. In the interest of overall paper length, we kept
this discussion to a minimum. The main point is that our model abstracts away
many of the fine grained model parameters of advanced tumour growth models. It
thereby offers an entry point for newcomers to the field to perform their first
steps in cancer modeling and to study the effect of choice of sampling positions
and sample size on the observed tumour profile in comparison to the whole tumour
profile. Relevant commits: 4a39b4
>> 

>     * [ ]  The summary states that simulated tumors can be subjected for multi-region sampling. Along my last major comment, please elaborate how this can be done in light of algorithm output.
Indeed, our code allows to generate mutation
profiles from multiple samples taken at different positions of the tumour. This
feature was previously hidden from the user, i.e. the user would have to
hardcode the coordinates into the code `casim.py`. We have now exposed the
sampling coordinates in the public API through the parameter
`sampling_positions` in `CancerSimulationParameters`. The default behaviour
corresponds to the unrevised version of the code, where one sampling
position would be chosen randomly.
Relevant commits: 839c23a3, 7e32183

> 
>     * [ ]  In the fifth paragraph of the summary, please change "In each simulation step, every tumour cell in the tumour that has an unoccupied neighbour" to "In each simulation step, every tumour cell **on the lattice** that has an unoccupied neighbour"
done.
> 
>     * [ ]  The fifth paragraph of the summary also states that variability in fitness can be introduced for some cells. Please be more specific: How and for which cells can varaibility in fitness be introduced? How can beneficial and deleterious mutation rates be encoded?
All model parameters for tumour growth, `division_probability`,
`mutation_probability`, and `death_probability` have a counterpart for cells
that carry a beneficial or detrimental mutation,
`adv_mutant_division_probability`, `adv_mutant_mutation_probability`, and
`adv_mutant_death_probability`. The role of each of these parameters is
explained in the README, the table in the paper, in the online documentation, in
the template parameter module `params.py` and in the API reference manual for
the `CancerSimulationParameters` class. Additionally and as explained above, the
user can change selected parameters by snapshooting and reloading the run. This
is now documented by a new example notebook which is referenced in the paper and
in the README.
> 
>     * [ ]  Can other sources of sequencing noise be introduced (e.g. sequencing errors at non-mutated sites)?
While adding other sources of sequencing noise is certainly possible within our
model, it would go beyond the scope of our objective for this version. Here,
we consider sequencing noise at mutated sites as the major contribution and
neglect all non-mutated sites as minorly important in all cases except single
cell sampling. We will consider more
sequencing noise in subsequent versions of the code. A note in this regard was
added to the paper (Summary, last paragraph).



