# Compagnion code and script

This folder contains the script and all the data used in the manuscript "Essential ingredients in Joint Species Distribution Models: how to optimise inference and prediction in species rich communities?"

- Scripts prefixed with 00 correspond to the pre-processing of the data before giving it to the models.

- Scripts prefixed with 01 correspond to the HMSC models. Beware, each of these scripts requires a lot of computing power. They took about 3 weeks to run on an HPC. Thus, we have saved the results of the models so that they can be more easily analysed.

- The R script prefixed with 03 contain the code to run the diagnostics on the MCMC chains.

- Finally, the script prefixed with 04 is used to generate all the figures in the manuscript. To do so, it uses the output saves of the models, as specified above.