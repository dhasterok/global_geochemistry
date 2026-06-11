function protolith_prep(data);

ind = strfind(pwd,'protolith');
if ~isempty(ind)
    filename = 'protolith/';
else
    filename = '';
end
filename = [filename,'classifier_data_',datestr(date,29)];

% unequal data
%[traindata,testdata,scores] = prep_for_cluster(data, ...
%    'subset',3,'reserve',0.1,'fields',{'major'});
%save([filename,'_raw.mat'],'traindata','testdata');

%[traindata,testdata,scores] = prep_for_cluster(data, ...
%    'subset',3,'reserve',0.1,'fields',{'major'},'recenter','clr');
%save([filename,'_clr.mat'],'traindata','testdata');

[traindata,testdata,scores] = prep_for_cluster(data, ...
    'subset',3,'reserve',0.1,'fields',{'major'},'recenter','ilr');
save([filename,'_ilr2.mat'],'traindata','testdata');


% pca filtered unequal data (testdata are not pca filtered)
%[traindata,testdata,scores] = prep_for_cluster(data, ...
%    'subset',3,'reserve',0.1,'fields',{'major'},'pca');
%save([filename,'_raw_pca.mat'],'traindata','testdata');

%[traindata,testdata,scores] = prep_for_cluster(data, ...
%    'subset',3,'reserve',0.1,'fields',{'major'},'recenter','clr','pca');
%save([filename,'_clr_pca.mat'],'traindata','testdata');

%[traindata,testdata,scores] = prep_for_cluster(data, ...
%    'subset',3,'reserve',0.1,'fields',{'major'},'recenter','ilr','pca');
%save([filename,'_ilr_pca.mat'],'traindata','testdata');


% equal data
%[traindata,testdata,scores] = prep_for_cluster(data, ...
%    'subset',3,'reserve',0.1,'fields',{'major'},'equal');
%save([filename,'_raw_equal.mat'],'traindata','testdata');

%[traindata,testdata,scores] = prep_for_cluster(data, ...
%    'subset',3,'reserve',0.1,'fields',{'major'},'recenter','clr','equal');
%save([filename,'_clr_equal.mat'],'traindata','testdata');

%[traindata,testdata,scores] = prep_for_cluster(data, ...
%    'subset',3,'reserve',0.1,'fields',{'major'},'recenter','ilr','equal');
%save([filename,'_ilr_equal.mat'],'traindata','testdata');

return