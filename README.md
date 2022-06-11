# global_geochemistry
Codes for analyzing global geochemistry

Description:
This repository contains whole rock geochemical data for a global set of samples from a variety of governmental and academic studies.


Requirements:
MATLAB is required to run .m codes and load any .mat files

Some of the codes utilize geological province boundaries to assign metadata that can be used to examine chemical patterns in different tectonic environments.  These models are contained within a separate GitHub repository dhasterok/global_tectonics.

File structure:


global_geochemistry/processing


global_geochemistry/classification
  Contains codes that classify whole rock samples on the basis of chemistry
  
  
global_geochemistry/physprop
  Codes that estimate physical properties from geochemistry.
  
  
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


global_geochemistry/plotting
  Codes for visualizing the dataset.


global_geochemistry/ternary
  Make ternary scatter, gridded, or contour plots of data on ternary diagrams.


global_geochemistry/maptools
  Codes for making maps.


global_geochemistry/temporal
  Codes that deal specifically with temporal processing of the dataset.


global_geochemistry/stats
  Codes for statistical analysis of the data, generally non-routine methods including appropriately dealing with censored (below detection) data.
  
  
global_geochemistry/ref_models
  Contains geochemical reference models.
  
  
The repository is currently under construction.

For issues or questions, contact:

Derrick Hasterok
University of Adelaide
+61 4 8313 4540
derrick.hasterok@adelaide.edu.au
