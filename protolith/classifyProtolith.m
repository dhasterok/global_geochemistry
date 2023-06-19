function [protolith_est,labels,protolithScore] = classifyProtolith(data,trainFcn)

% copy true rock group to temporary variable, otherwise it will use the
% true rock group if it exists as one of the classification 
%true_rg = data.rock_group;
%data.rock_group(:) = {''};
%import classreg.learning.classif.CompactClassificationEnsemble/predict

load(['protolith/',trainFcn]);

oxlist = {'SiO2'; 'TiO2'; 'Al2O3'; 'FeO'; 'MgO'; 'CaO'; 'Na2O'; 'K2O'; 'P2O5'};

protolith_est = protolithClassifier.predictFcn(data);
[labels,protolithScore] = predict(protolithClassifier.ClassificationEnsemble,data);
%data.rock_group = true_rg;


% Make a table with results for protolith classifier
fprintf('                                   predicted protolith\n');
fprintf('                                 ------------------------\n');
fprintf('                       total      igneous    sedimentary\n');

ind = rockgroup(data,'igneous');
fprintf('igneous                  %7i   %7i   %7i \n', ...
    sum(ind), ...
    sum(ind & strcmp('igneous',protolith_est)), ...
    sum(ind & strcmp('sedimentary',protolith_est)));

ind = rockgroup(data,'sedimentary');
fprintf('sedimentary              %7i   %7i   %7i \n', ...
    sum(ind), ...
    sum(ind & strcmp('igneous',protolith_est)), ...
    sum(ind & strcmp('sedimentary',protolith_est)));

ind = rockgroup(data,'metaigneous');
fprintf('metaigneous              %7i   %7i   %7i \n', ...
    sum(ind), ...
    sum(ind & strcmp('igneous',protolith_est)), ...
    sum(ind & strcmp('sedimentary',protolith_est)));

ind = rockgroup(data,'metasedimentary');
fprintf('metasedimentary          %7i   %7i   %7i \n', ...
    sum(ind), ...
    sum(ind & strcmp('igneous',protolith_est)), ...
    sum(ind & strcmp('sedimentary',protolith_est)));

ind = rockgroup(data,'metamorphic');
fprintf('metamorphic (unreported) %7i   %7i   %7i \n', ...
    sum(ind), ...
    sum(ind & strcmp('igneous',protolith_est)), ...
    sum(ind & strcmp('sedimentary',protolith_est)));

return