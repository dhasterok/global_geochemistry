function data = ecread3(filename)
% ECREAD - Reads Earthchem xls output
%
%   DATA = ECREAD(FILENAME) where FILENAME is an excel spreadsheet containing
%   earthchem.org data.  The output, DATA, is a structure array containing the
%   earthchem data with each field being a column from the original.

[~,~,data] = xlsread(filename);

return
