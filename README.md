# global_geochemistry
A global geochemistry dataset for rocks

Description:
This repository contains whole rock geochemical data for a global set of samples from a variety of governmental and academic studies.


Requirements:
MATLAB is required to run .m codes and load any .mat files


File structure:

global_geochemistry/classification
  Contains codes that classify whole rock samples on the basis of chemistry
  

global_geochemistry/protolith
  Contains codes associated with protolith classification, including the preparation of the dataset, training, validation, and analysis
  
  protolith_prep.m - sets the parameters of the dataset to be extracted
  prep_for_cluster.m - extracts the data, recenters if necessary and makes some plots of the data distribution
  
  train_RUSBoost_Classifier_30l_1000s_20190222.m - trains a classifier for protolith classes (igneous/sedimentary)
  trained_RUSBoost_Model_30l_1000s_20190226.mat - trained classifier

  classification_results_20190215.xlsx - a table of classification results from a number of machine learning algorithms
  training_analysis.m - plots the results contained within classification_results_(file date).xlsx
  
  tree_analysis.m - produces a plot of classification results for a training and validation dataset using the classifier
    trained_RUSBoost_Model_(classifier parameters).mat
  misclassify_plot.m - produces a plot of classification results by rock type
  

The repository is currently under construction.

For issues or questions, contact:

Derrick Hasterok
University of Adelaide
+61 4 8313 4540
derrick.hasterok@adelaide.edu.au
