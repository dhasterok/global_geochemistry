function protolith_predictor
% PROTOLITH_PREDICTOR - predicts a protolith from major element oxides.
%
%   protolith_predictor will open a dialog to select an infile '*.xls'
%   or '*.xlsx to process and predict a protoliths.  The function will
%   output the name and path to the outfile.
%
%   protolith_predictor(infile) can be used to specify the infile directly
%   instead of using the dialog to open the infile.
%
%   protolith_predictor(path,infile) can be used to open the infile if it
%   is in a different path than the working directory.
%
%   Examples:
%       % to open a dialog to select the file of geochemical data
%       protolith_predictor 
%
%       % or to specify the geochemical data file directly
%       protolith_predictor('xlsx/','protolith_template.xlsx');

% Original: 2019 June 18 by D. Hasterok dhasterok@gmail.com

% file with protolithClassifier function;
% note: for this function to work, the protolith classification file must
% be downloaded from https://dx.doi.org/10.5281/zenodo.2586461
pclass = 'trained_RUSBoost_Model_30l_1000s_20190226.mat';

% convert chemistry file to format needed for protolith classification
outfile = proto_xls2csv;

% reading outfile
rawdata = readtable(outfile);

% calculate total iron and iron ratio
fprintf('Convert all Fe to FeO and calculate Fe2+/Fe_total ratio...\n');
%sum(data{:,'fe2o3_tot'} > 0 | data{:,'fe2o3'} > 0 | data{:,'feo'} > 0 | data{:,'feo_tot'} > 0)
processdata = fefix(rawdata);
%sum(data{:,'feo_tot'} > 0)

% convert cation data to oxides
fprintf('Convert cations to oxide data when missing...\n');
processdata = cat2ox(processdata);

% normalize oxide weights
fprintf('Normalizing oxides...\n');
oxlist = {'SiO2'; 'TiO2'; 'Al2O3'; 'FeO'; ...
    'MgO'; 'CaO'; 'Na2O'; 'K2O'; 'P2O5'};
processdata = oxide_norm(processdata,oxlist);

% classify protoliths
fprintf('Classifying protolith using %s...\n',pclass);
processdata.protolithClass = cell([height(processdata) 1]);
processdata.protolithClass(:) = {''};

load(pclass);
[processdata.protolithClass,protolithScore] = ...
    predict(protolithClassifier.ClassificationEnsemble,processdata);

% normalizing factor (This is output from tree_analysis.m as processed on
% the global dataset
x = 0.926879;
processdata.protolithScore = protolithScore(:,1)/x-1;

t = processdata(:,{'sample_name','protolithClass','protolithScore'});

% write out protolith estimates
pfile = [outfile(1:end-4),'_classified.csv'];
fprintf('Writing file of protolith predictions... %s\n',pfile);
writetable(t,pfile);

return