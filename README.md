# Jeffreys_etal
covariates on trends for Rufous Hummingbird

This repo holds the code and data necessary to replicate the BBS component of Jeffreys et al. Detecting Rufous Hummingbird habitat change.

There are two key R scripts:
Species_data_prep_mean_habitat_and_slope.R - this script prepares the BBS data for Rufous Hummingbird, selecting the routes on which the species has been observed, setting up the spatial relationships among routes, and joining the BBS observations at each route to the annual estimates of habitat suitability. All data are then combined into an R list object that will form the data object for the Stan model fitting functions in the following script.
Fig_final_2006_2021_model.R - this script loads the data object prepared in the previous script, fits the Stan model, tests for convergence, and saves the stored model output.

The final results are explained in the quarto document that defines the model, summarises the results and includes maps of the predicted trends and abundances.


