function tree_analysis(classifier,varargin)

if nargin == 3
    traindata = varargin{1};
    testdata = varargin{2};
elseif nargin == 2
    load([varargin{1},'.mat']);
end

% determine rock types from chemistry
%traindata = rock_class(traindata);
%testdata = rock_class(testdata);

% produce tables and figures
analysis(traindata,testdata,classifier);

return



function analysis(traindata,testdata,vrsn,method)

% load classifier
load([vrsn,'.mat']);

% predict protolith type
traindata.protolith = protolithClassifier.predictFcn(traindata);
testdata.protolith = protolithClassifier.predictFcn(testdata);

% scores
if isfield(protolithClassifier,'ClassificationEnsemble')
    [traindata.label,traindata.score] = predict(protolithClassifier.ClassificationEnsemble,traindata);
    [testdata.label,testdata.score] = predict(protolithClassifier.ClassificationEnsemble,testdata);
elseif isfield(protolithClassifier,'ClassificationKNN')
    [traindata.label,traindata.score] = predict(protolithClassifier.ClassificationKNN,traindata);
    [testdata.label,testdata.score] = predict(protolithClassifier.ClassificationKNN,testdata);
else
    error('Unknown classification scheme.  protolithClassifier.ClassificationXXX');
end

% write table
fprintf('       & & predicted protolith\n');
fprintf('       & total & igneous & sedimentary\n');

fprintf('{\\it training dataset}\n');
ind = rockgroup(traindata,'igneous');
fprintf('igneous & %7i & %7i & %7i \n', ...
    sum(ind), ...
    sum(ind & strcmp('igneous',traindata.label)), ...
    sum(ind & strcmp('sedimentary',traindata.label)));

ind = rockgroup(traindata,'sedimentary');
fprintf('sedimentary & %7i & %7i & %7i \n', ...
    sum(ind), ...
    sum(ind & strcmp('igneous',traindata.label)), ...
    sum(ind & strcmp('sedimentary',traindata.label)));

fprintf('{\\it validation dataset}\n');
ind = rockgroup(testdata,'igneous');
fprintf('igneous & %7i & %7i & %7i \n', ...
    sum(ind), ...
    sum(ind & strcmp('igneous',testdata.label)), ...
    sum(ind & strcmp('sedimentary',testdata.label)));

ind = rockgroup(testdata,'sedimentary');
fprintf('sedimentary & %7i & %7i & %7i \n', ...
    sum(ind), ...
    sum(ind & strcmp('igneous',testdata.label)), ...
    sum(ind & strcmp('sedimentary',testdata.label)));

figure;
% histograms of training dataset scores (normalized to 1)
x = min(traindata.score(strcmp(traindata.label,'igneous'),1));
subplot(221);
histogram(traindata.score(rockgroup(traindata,'igneous'),1)/x-1,'BinWidth',0.05,'DisplayStyle','stairs');
hold on;
plot([0 0],get(gca,'YLim'));
xlabel('Normalized Score');
title('Igneous (train)');
xlim([-1 1]);

subplot(222);
histogram(traindata.score(rockgroup(traindata,'sedimentary'),2)/x-1,'BinWidth',0.05,'DisplayStyle','stairs');
hold on;
plot([0 0],get(gca,'YLim'));
xlabel('Normalized Score');
title('Sedimentary (train)');
xlim([-1 1]);

% histograms of validation dataset scores (normalized to 1)
subplot(223);
histogram(testdata.score(rockgroup(testdata,'igneous'),1)/x-1,'BinWidth',0.05,'DisplayStyle','stairs');
hold on;
plot([0 0],get(gca,'YLim'));
xlabel('Normalized Score');
title('Igneous (validation)');
xlim([-1 1]);

subplot(224);
histogram(testdata.score(rockgroup(testdata,'sedimentary'),2)/x-1,'BinWidth',0.05,'DisplayStyle','stairs');
hold on;
plot([0 0],get(gca,'YLim'));
xlabel('Normalized Score');
title('Sedimentary (validation)');
xlim([-1 1]);

% figure of incorrect values for training dataset
misclassify_plot(traindata,traindata.label);

% figure of incorrect values for testing dataset
misclassify_plot(testdata,testdata.label);

return