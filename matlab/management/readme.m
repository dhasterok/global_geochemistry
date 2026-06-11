%**********************************************************
% Whole-rock geochemical database MATLAB manipulation guide
%               Matthew Gard, 01/08/2019
%**********************************************************
% V1.0.0 - Initial version
%**********************************************************

% Basic MATLAB knowledge is assumed

% Loading the data
%==================
% MATLAB's table format is an excellent variable format for manipulating
% large databases
% To load the table, simply use the readtable function
% e.g. data = readtable('filepath/filename.csv');

% Symbology used in this document
%=================================
% |  means OR
% &  means AND
% >  means greater than
% >= means greater than or equal to
% <  means less than
% <= means less than or equal to
% == means equals

% Filtering the data
%====================
% Filtering the data using indices for string type columns, or numeric
% columns
% Indices store a vector or 0 or 1's corresponding to True or False if the
% row satisfies some conditions. As we have the data loaded in table
% format, this makes filtering even easier
%
% For example:
%     Create an indices vector for sio2 >= 65 and sio2 <= 80
%     ind = data.sio2 >= 65 & data.sio2 <= 80;
%
%     Find all TAS classified basalts or granites
%     ind = strcmpi(data.rock_type,'granite') | strcmpi(data.rock_type,'basalt');
%
% You can then save a new table containing only these samples by using the
% indices vector you created
% e.g. 
%     ind = strcmpi(data.rock_type,'basalt');
%     basaltdata = data(ind,:);
