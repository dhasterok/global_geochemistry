function [p,ndata] = spider(refdata,sdata,normref,varargin);

addpath ref_models

% load reference chemistry table
eref = readtable('earthref.xlsx');
ind = (strcmpi(eref.model,normref) | strcmpi(eref.reference,normref)) & ...
    eref.sigma == 0;

opt = 1;
if ischar(normref)
    if sum(ind) > 1
        if nargin < 3
            error(['A layer is required to use ',normref,'.']);
        else
            layer = varargin{opt};
            ind = ind & strcmpi(eref.layer,layer);
            opt = opt + 1;
        end
    else
        layer = eref.layer{ind};
    end
    refchem = eref(ind,:);
else
    refchem = normref;
    normref = '';
    layer = 'user defined';
end

ellist = '';
Q = [0.05 0.25 0.5 0.95 1];
if nargin >= opt + 2
    while opt + 2 < nargin
        switch lower(varargin{opt})
            case 'elements'
                ellist = varargin{opt+1};
                opt = opt + 2;
            case 'quantiles'
                Q = varargin{opt+1};
                opt = opt + 2;
            otherwise
                error('Unknown option.');
        end
    end
end

[p,ndata] = spiderplot(sdata,refchem,ellist,Q);

str = sprintf('Normalized to %s',layer);
title(['Normalized to ',layer]);

str = sprintf('abundance/%s (%s)',layer,normref);
ylabel(str);

return


function [p,ndata] = spiderplot(data,refchem,list,Q)

if isempty(list)
    list = 'all';
end

if ischar(list)
    switch list
        case 'all'
            ellist = {'Cs', 'Rb', 'Ba', 'Th', ...
                'U', 'Nb', 'Ta', 'K', ...
                'La', 'Ce', 'Pb', 'Mo', ...
                'Pr', 'Sr', 'Ga', 'Zr', 'Hf', ...
                'Nd', 'Sm', 'Eu', 'Li', ...
                'Ti', 'Gd', 'Tb', 'Dy', ...
                'Y', 'Ho', 'Er', ...
                'Tm', 'Yb', 'Lu', 'Zn', ...
                'Mn', 'V', 'Sc', 'Co', ...
                'Cu', 'Ni', 'Cr'};
        case 'ree'
            ellist = {'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', ... 
                'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu'};

    end
else
    ellist = list;
end
       
for i = 1:length(ellist)
    field = [lower(ellist{i}),'_ppm'];
    
    if strcmp(field,'k_ppm')
        data.k_ppm = 2*molecularwt('K')/molecularwt('K2O')*data.k2o*10000;
    elseif strcmp(field,'ti_ppm')
        data.ti_ppm = molecularwt('Ti')/molecularwt('TiO2')*data.tio2*10000;
    elseif strcmp(field,'ba_ppm')
        data.ba_ppm = molecularwt('Ba')/molecularwt('BaO')*data.bao*10000;    
    elseif strcmp(field,'cr_ppm')
        data.cr_ppm = 2/3*molecularwt('Cr')/molecularwt('Cr2O3')*data.cr2o3*10000;
    elseif strcmp(field,'p_ppm')
        data.p_ppm = 2*molecularwt('P')/molecularwt('P2O5')*data.p2o5*10000;   
    end
    
    ind = data{:,field} > 0;
    el_norm = log10(data{ind,field}/refchem{1,field});
    
    el_stat(i,:) = quantile(el_norm,Q);
end
x = [1:length(ellist)]';
if length(Q) == 5
    p = plot(x,el_stat(:,[1 5]),'-','Color',[0.7 0.7 0.7]);
    hold on;
    p = [p; fill([x; flipud(x)],[el_stat(:,2); flipud(el_stat(:,4))],[0.8 0.8 0.8])];
    p = [p; plot(x,el_stat(:,3),'-')];
elseif length(Q) == 3
    p = fill([x; flipud(x)],[el_stat(:,1); flipud(el_stat(:,3))],[0.8 0.8 0.8]);
    hold on;
    p = [p; plot(x,el_stat(:,2),'-')];
else
    p = plot(x,el_stat,'-');
end
set(gca,'XTick',x,'XTickLabel',ellist);
xlim([x(1)-1 x(end)+1]);

ndata = array2table(10.^el_stat','VariableNames',ellist);

return
