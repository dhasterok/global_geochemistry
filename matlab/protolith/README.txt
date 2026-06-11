This directory contains data, codes and analysis for the classification of metamorphic protoliths using machine learning methods.  Only the preferred classifier is included in this directory.  The results included in this directory are associated with the paper by Hasterok et al., (Computers & Geosciences, submitted).

This README.txt file and the associated files can be found at the Github repository (http://github.com/dhasterok/global_geochemistry/tree/master/protolith/).

List of files:
--------------------------------------
classification_results_20190215.xlsx - a list of models tested the results on a test dataset.
misclassify_plot.m - creates a plot of results for the preferred classifier (RUSBoost, no PCA, with 30 learners and 1000 splits).
prep_for_cluster.m - A Matlab function that extracts data from the global geochemical dataset by Gard et al. (Earth System Sci. Data, submitted) that is used for training and validation of the protolithClassifier.
proto_xls2csv.m - A Matlab function used to format an Excel spreadsheet for protolith classification.
protolith_plots.mlx - A Matlab live script that contains the code to recreate figures from Hasterok et al. (Computers & Geosciences, submitted)
protolith_predictor.m - A Matlab function for predicting protolith estimates from a geochemical file (Excel spreadsheet). note: for this function to work, the protolith classification file must be downloaded from https://dx.doi.org/10.5281/zenodo.2586461
protolith_prep.m - A Matlab function that calls prep_for_cluster.m with selected options (e.g., data transformation, PCA, reserve dataset size).
protolith_template.xlsx - An example geochemistry input file.
train_RUSBoost_Classifier_30l_1000s_20190222.m - A Matlab function that can be used to train the a new protolithClassifier using the same options as presented in Hasterok et al. (Computers & Geosciences, revised).
trained_RUSBoost_Model_30l_1000s_20190226.mat - A Matlab datafile that contains the data tables used to train the protolithClassifier and post-training validation.
training_analysis.m - A Matlab function that creates figures illustrating model performance from the tested machine learning algorithms and options.
tree_analysis.m - A Matlab function that examines in detail the preferred classifier performance.

Associated files:
--------------------------------------
(symbolic links included in this directory)
fefix.m - in ../processing folder (converts all FeO and Fe2O3 to FeO total)
cat2ox.m - in ../processing folder (converts cations to oxides if necessary)
molecularwt.m - in ../toolbox folder (computes molecular weights given a chemical formula)
oxide_norm.m - in ../processing folder (normalizes oxides)

Associated geochemical database files:
--------------------------------------
The analyses used to train the metamorphic classifiers and the references to the geochemical data are available at Zenodo repository (https://dx.doi.org/10.5281/zenodo.2586461).

Reference:
--------------------------------------
Hasterok, D., M. Gard, C.M.B. Bishop, and D. Kelsey, (revised), Chemical identification of metamorphic protoliths using machine learning methods, Computers & Geosciences.
Gard, M., Hasterok, D., Halpin, J., 2019. Global whole-rock geochemical database compilation. Earth System Science Data Discussions, Earth System Science Data Discussions 1–23. doi:10.5194/essd-2019-50.
