This directory contains data, codes and analysis for the classification of metamorphic protoliths using machine learning methods.  Only the preferred classifier is included in this directory.  The results included in this directory are associated with the paper by Hasterok et al., (Computers & Geosciences, submitted).

This README.txt file and the associated files can be found at the Github repository (http://github.com/dhasterok/global_geochemistry/tree/master/protolith/).

List of files:
--------------------------------------
classification_results_20190215.xlsx - a list of models tested the results on a test dataset.
misclassify_plot.m - creates a plot of results for the preferred classifier (RUSBoost, no PCA, with 30 learners and 1000 splits).
prep_for_cluster.m - A Matlab function that extracts data from the global geochemical dataset by Gard et al. (Earth System Sci. Data, submitted) that is used for training and validation of the protolithClassifier.
protolith_plots.mlx - A Matlab live script that contains the code to recreate figures from Hasterok et al. (Computers & Geosciences, submitted)
protolith_prep.m - A Matlab function that calls prep_for_cluster.m with selected options (e.g., data transformation, PCA, reserve dataset size).
train_RUSBoost_Classifier_30l_1000s_20190222.m - A Matlab function that can be used to train the a new protolithClassifier using the same options as presented in Hasterok et al. (Computers & Geosciences, submitted).
trained_RUSBoost_Model_30l_1000s_20190226.mat - A Matlab datafile that contains the data tables used to train the protolithClassifier and post-training validation.
training_analysis.m - A Matlab function that creates figures illustrating model performance from the tested machine learning algorithms and options.
tree_analysis.m - A Matlab function that examines in detail the preferred classifier performance.

Associated geochemical database files:
--------------------------------------
The analyses used to train the metamorphic classifiers and the references to the geochemical data are available at Zenodo repository (https://dx.doi.org/10.5281/zenodo.2586461).

Reference:
--------------------------------------
Hasterok, D., M. Gard, C.M.B. Bishop, and D. Kelsey, (submitted to Computers & Geosciences), Chemical identification of metamorphic protoliths using machine learning methods.
