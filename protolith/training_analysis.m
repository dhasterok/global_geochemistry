function training_analysis

restable = readtable('classification_results_20190215.xlsx');

figure;

ind = restable.pca == 0 & strcmp(restable.transform,'raw') & ...
    strcmp(restable.setsize,'equal');
scatter(restable.igneous_tp_accuracy(ind),restable.sedimentary_tp_accuracy(ind),16,'filled','o');
hold on;

ind = restable.pca == 0 & strcmp(restable.transform,'raw') & ...
    strcmp(restable.setsize,'proportional') & strcmp(restable.algorithm,'RUSBoosted');
scatter(restable.igneous_tp_accuracy(ind),restable.sedimentary_tp_accuracy(ind),16,'filled','o');

ind = restable.pca == 0 & strcmp(restable.transform,'raw') & ...
    strcmp(restable.setsize,'proportional') & ...
    ~(strcmp(restable.algorithm,'weighted') | strcmp(restable.algorithm,'RUSBoosted'));
scatter(restable.igneous_tp_accuracy(ind),restable.sedimentary_tp_accuracy(ind),16,'filled','o');

ind = restable.pca == 0 & strcmp(restable.transform,'raw') & ...
    strcmp(restable.setsize,'proportional') & strcmp(restable.algorithm,'weighted');
scatter(restable.igneous_tp_accuracy(ind),restable.sedimentary_tp_accuracy(ind),16,'filled','o');

axis equal;
axis tight
xlabel('igneous TP/(true igneous)');
ylabel('sedimentary TP/(true sedimentary)');
legend('equal','RUSBoosted proportional','proportional','KNN weighted proportional');
ylim([0 1]);
xlim([0.7 1]);
set(gca,'Box','on','XTick',[0:0.1:1],'YTick',[0:0.1:1]);


indx = [];
indy = [];
for i = 1:height(restable)
    if restable.pca(i) == 0 & strcmp(restable.setsize(i),'proportional')
        ix = i;
        
        iy = find(restable.pca == 1 & ...
            strcmp(restable.setsize,'proportional') & ...
            strcmp(restable.method,restable.method(i)) & ...
            strcmp(restable.algorithm,restable.algorithm(i)) & ...
            strcmp(restable.transform,restable.transform(i)) & ...
            restable.splits == restable.splits(i) & ...
            restable.learners == restable.learners(i));
        iy = iy(iy ~= i);
        
        if length(iy) > 1
            restable(ix,:)
            restable(iy,:)
            pause
        end
        if ~isempty(iy)
            indx = [indx; ix];
            indy = [indy; iy];
        end
    end
end

figure;
subplot(121);
hold on;
plot([0 1],[0 0],'k-');
scatter(restable.normalized_accuracy(indx), ...
    (restable.normalized_accuracy(indy)-restable.normalized_accuracy(indx)) ...
    ./restable.normalized_accuracy(indx)*100,20,'filled','ko');
xlabel('Normalized accuracy (unfiltered)');
ylabel('% Difference');
text(0.45,7,'(a) PCA: (filtered - unfiltered)/unfiltered');
axis square;
axis([0.4 1 -10 10]);
set(gca,'Box','on','XTick',[0.4:0.1:1]);

indx = [];
indy = [];
indz = [];
for i = 1:height(restable)
    if restable.pca(i) == 0 & strcmp(restable.setsize(i),'proportional') ...
            & strcmp(restable.transform(i),'raw')
        ix = i;
        
        iy = find(restable.pca == 0 & ...
            strcmp(restable.setsize,'proportional') & ...
            strcmp(restable.method,restable.method(i)) & ...
            strcmp(restable.algorithm,restable.algorithm(i)) & ...
            strcmp(restable.transform,'clr') & ...
            restable.splits == restable.splits(i) & ...
            restable.learners == restable.learners(i));
        iy = iy(iy ~= i)';
        
        iz = find(restable.pca == 0 & ...
            strcmp(restable.setsize,'proportional') & ...
            strcmp(restable.method,restable.method(i)) & ...
            strcmp(restable.algorithm,restable.algorithm(i)) & ...
            strcmp(restable.transform,'ilr') & ...
            restable.splits == restable.splits(i) & ...
            restable.learners == restable.learners(i));
        iz = iz(iz ~= i)';
        
        if ~isempty(iy)
            indx = [indx; ix];
            indy = [indy; iy];
            indz = [indz; iz];
        end
        if ~isempty(iz)
            %indxx = [indx; ix];
            %indz = [indz; iz];
        end
    end
end

subplot(122);
hold on;
plot([0 1],[0 0],'k-');
p(1) = scatter(restable.normalized_accuracy(indx), ...
    (restable.normalized_accuracy(indy)-restable.normalized_accuracy(indx)) ...
    ./restable.normalized_accuracy(indx)*100,20,'filled','ko');
p(2) = scatter(restable.normalized_accuracy(indx), ...
    (restable.normalized_accuracy(indz)-restable.normalized_accuracy(indx)) ...
    ./restable.normalized_accuracy(indx)*100,20,'ko');
legend(p,'clr','ilr');
xlabel('Normalized accuracy (Aitchison geometry)');
ylabel('% Difference');
text(0.45,7,'(b) lr-transformed: (lr - Aitchison)/Aitchison');
axis square;
axis([0.4 1 -10 10]);
set(gca,'Box','on','XTick',[0.4:0.1:1]);


figure;
hold on;

ind = restable.pca == 0 & strcmp(restable.transform,'raw') & ...
    strcmp(restable.setsize,'proportional') & strcmp(restable.algorithm,'Bagging');

ensembleresults(restable,ind,'normalized_accuracy');
ensembleresults(restable,ind,'igneous_tp_accuracy');
ensembleresults(restable,ind,'sedimentary_tp_accuracy');

title('Bagging');
hpax([1 4],'x');
axis square;
xlabel('Branches/Trees');
legend('Branches, Trees = 30','Trees, Branches = 20','Branches = Trees','Trees, Branches = 1000','Location','NorthWest');

figure;
hold on;
ind = restable.pca == 0 & strcmp(restable.transform,'raw') & ...
    strcmp(restable.setsize,'proportional') & strcmp(restable.algorithm,'AdaBoost');

ensembleresults(restable,ind,'normalized_accuracy');
ensembleresults(restable,ind,'igneous_tp_accuracy');
ensembleresults(restable,ind,'sedimentary_tp_accuracy');

title('AdaBoost');
hpax([1 4],'x');
axis square;
xlabel('Branches/Trees');
legend('Branches, Trees = 30','Trees, Branches = 20','Branches = Trees','Trees, Branches = 1000','Location','NorthWest');

figure;
hold on;
ind = restable.pca == 0 & strcmp(restable.transform,'raw') & ...
    strcmp(restable.filename,'classifier_data_2019-02-11_ilr.mat') & ...
    restable.learning_rate == 0.1 & ...
    strcmp(restable.setsize,'proportional') & strcmp(restable.algorithm,'RUSBoosted');

ensembleresults(restable,ind,'normalized_accuracy');
ensembleresults(restable,ind,'igneous_tp_accuracy');
ensembleresults(restable,ind,'sedimentary_tp_accuracy');

title('RUSBoost (2019-02-11)');
hpax([1 4],'x');
axis square;
xlabel('Branches/Trees');
legend('Branches, Trees = 30','Trees, Branches = 20','Branches = Trees','Trees, Branches = 1000','Location','NorthWest');

figure;
hold on;
ind = restable.pca == 0 & strcmp(restable.transform,'raw') & ...
    strcmp(restable.filename,'classifier_data_2019-02-26_ilr.mat') & ...
    strcmp(restable.setsize,'proportional') & strcmp(restable.algorithm,'RUSBoosted');

ensembleresults(restable,ind,'normalized_accuracy');
ensembleresults(restable,ind,'igneous_tp_accuracy');
ensembleresults(restable,ind,'sedimentary_tp_accuracy');
title('RUSBoost (2019-02-26)');

hpax([1 4],'x');
axis square;
xlabel('Branches/Trees');
legend('Branches, Trees = 30','Trees, Branches = 20','Branches = Trees','Trees, Branches = 1000','Location','NorthWest');

return


function ensembleresults(restable,ind,field)

switch field
    case 'normalized_accuracy'
        m = 'o';
    case 'igneous_tp_accuracy'
        m = '^';
    case 'sedimentary_tp_accuracy'
        m = 'v';
end
C = [0.000   0.447   0.741;
    0.850   0.325   0.098;
    0.929   0.694   0.125;
    0.494   0.184   0.556];
[x,ix] = sort(log10(restable.splits(ind & restable.learners == 30)));
y = restable{ind & restable.learners == 30,field};
y = y(ix);
plot(x,y,'-','Color',C(1,:));
scatter(x,y,20,C(1,:),'filled',m);

[x,ix] = sort(log10(restable.learners(ind & restable.splits == 20)));
y = restable{ind & restable.splits == 20,field};
y = y(ix);
plot(x,y,'-','Color',C(2,:));
scatter(x,y,20,C(2,:),'filled',m);

[x,ix] = sort(log10(restable.splits(ind & restable.splits == restable.learners)));
y = restable{ind & restable.splits == restable.learners,field};
y = y(ix);
plot(x,y,'-','Color',C(3,:));
scatter(x,y,20,C(3,:),'filled',m);

[x,ix] = sort(log10(restable.learners(ind & restable.splits == 1000)));
y = restable{ind & restable.splits == 1000,field};
y = y(ix);
plot(x,y,'-','Color',C(4,:));
scatter(x,y,20,C(4,:),'filled',m);

return