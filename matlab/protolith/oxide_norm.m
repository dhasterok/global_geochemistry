function data = oxide_norm(data,varargin)
% OXIDE_NORM - Computes the oxide norm for a set of given oxides
%
%   data = oxide_norm(data,oxides) computes the oxide norm for the oxides in the
%   given cell array 'oxides'.  The normalization of iron is set by the list 'oxides'.
%
%   Normalization is performed to total iron converted to FeO.

% Oxide lists
% ---------------------------------
% all avaliable oxides
oxides = {'sio2','tio2','al2o3','feo_tot', ...
    'mgo','cao','na2o','k2o','p2o5', ... % common oxides
    'cr2o3','mno','nio','bao'};%,'sro'}; % less common oxides

% all volatiles
volatiles = {'h2o_tot','co2','so3','f_ppm','cl_ppm'};
% conversion factor for volatiles
cf = [1 1 1 1e-4 1e-4];


normlist = {};
if nargin > 1
    for i = 1:nargin-1
        if iscell(varargin{i})
            normlist = lower(varargin{i});
        else
            error('Unknown option.');
        end
    end
end

for i = 1:length(oxides)
    if ~any(strcmp(data.Properties.VariableNames,oxides{i}))
        data{:,oxides{i}} = nan([height(data) 1]);
    end
end
if ~any(strcmp(data.Properties.VariableNames,'loi'))
    data.loi = nan([height(data) 1]);
end


% fixes norm oxide list to lowercase and feo to feo_tot
normlist = lower(normlist);
if any(strcmp('feo',normlist))
    normlist(strcmp('feo',normlist)) = {'feo_tot'};
end


% Prep data
% ---------------------------------
% convert caco3 to cao + co2

% if any(strcmp(data.Properties.VariableNames,'caco3'))
%     ind = data.caco3 > 0;
%     data.cao(ind) = nansum([data.cao(ind),data.caco3(ind)*molecularwt('CaO')/molecularwt('CaCO3')],2);
%     data.co2(ind) = nansum([data.co2(ind),data.caco3(ind)*molecularwt('CO2')/molecularwt('CaCO3')],2);
%     data.caco3 = [];
% end

% MATT EDIT:
% If it has caco3, co2, and cao -> co2 and cao are already converted in
% these publications, we dont need to add it on
% If it has caco3, and co2, but no cao -> theres only 1, ignore it
% If it has caco3 and cao, convert to co2 -> these are inconsistent, theres
% 1187 of them (out of 30000ish with caco3). They vary. Just ignore them?
% If it has caco3 and neither cao and co2 (30000ish) -> fill cao and co2
if any(strcmp(data.Properties.VariableNames,'caco3'))
    % Only fill co2 and cao where theyre empty and caco3 exists
    ind = data.caco3 > 0 & data.co2 < 0 & data.cao < 0;
    data.cao(ind) = data.caco3(ind)*molecularwt('CaO')/molecularwt('CaCO3');
    data.co2(ind) = data.caco3(ind)*molecularwt('CO2')/molecularwt('CaCO3');
    data.caco3 = [];
end

% convert mgco3 to mgo + co2
% if any(strcmp(data.Properties.VariableNames,'mgco3'))
%     ind = data.mgco3 > 0;
%     data.mgo(ind) = nansum([data.mgo(ind),data.mgco3(ind)*molecularwt('MgO')/molecularwt('MgCO3')],2);
%     data.co2(ind) = nansum([data.co2(ind),data.mgco3(ind)*molecularwt('CO2')/molecularwt('MgCO3')],2);
% 
%     data.mgco3 = [];
% end
% MATT EDIT: Same as above
if any(strcmp(data.Properties.VariableNames,'mgco3'))
    ind = data.mgco3 > 0 & data.co2 < 0 & data.mgo < 0;
    data.mgo(ind) = data.mgco3(ind)*molecularwt('MgO')/molecularwt('MgCO3');
    data.co2(ind) = data.mgco3(ind)*molecularwt('CO2')/molecularwt('MgCO3');
    data.mgco3 = [];
end

% add h2o+ and h2o- to yield h2o_tot where h2o_tot does not exist
if all(~strcmp('h2o_tot',data.Properties.VariableNames))
    data{:,{'h2o_tot'}} = nan([height(data) 1]);
end
if all(~strcmp('h2o_plus',data.Properties.VariableNames))
    data{:,{'h2o_plus'}} = nan([height(data) 1]);
end
if all(~strcmp('h2o_minus',data.Properties.VariableNames))
    data{:,{'h2o_minus'}} = nan([height(data) 1]);
end
if all(~strcmp('co2',data.Properties.VariableNames))
    data{:,{'co2'}} = nan([height(data) 1]);
end
if all(~strcmp('so3',data.Properties.VariableNames))
    data{:,{'so3'}} = nan([height(data) 1]);
end
if all(~strcmp('cl_ppm',data.Properties.VariableNames))
    data{:,{'cl_ppm'}} = nan([height(data) 1]);
end
if all(~strcmp('f_ppm',data.Properties.VariableNames))
    data{:,{'f_ppm'}} = nan([height(data) 1]);
end
ind = ~isnan(data.h2o_tot) | data.h2o_tot < 0;
data.h2o_tot(ind) = nansum([data.h2o_plus(ind),data.h2o_minus(ind)],2);

% compute loi from volatiles if loi is not given, probably a minimum loi
loi = nansum(repmat(cf,height(data),1).*data{:,volatiles},2);

% CO2 may be excluded from LOI if the rock is a carbonate.
ind = 1.2 * data.loi < data.co2; %20% variance between CO2 and total
data.loi(ind) = data.loi(ind) + data.co2(ind);

loi(loi == 0) = NaN;
ind = isnan(data.loi) & ~isnan(loi);

data.loi(ind) = loi(ind);

data.computed_loi = logical(zeros([height(data) 1]));
data.computed_loi(ind) = 1;

%fprintf('LOI <10 wt.%%: %i,   >10 wt.%%: %i\n',sum(data.loi < 10),sum(data.loi > 10));

subplot(221);
histogram(data.loi,'BinEdges',[-10:2:110]);
xlabel('LOI (wt.%)');
axis square;

% determine normalization factor
nf = ones([height(data) 1]);
if length(normlist) > 0
    fprintf('Nomalizing to oxides (normlist): \n');
    for i = 1:length(normlist)
        fprintf(' %s',normlist{i});
    end
    fprintf('...\n');
    
    % extract oxide data and set negative values to zero for normalization
    oxdata = data{:,normlist};
    oxdata(oxdata < 0) = 0;
    
    % compute total of oxides
    oxtotal = nansum(oxdata,2);
    
    ind = oxtotal > 0;
    
    % compute normalization factor
    %if data.cao + data.mgo > data.sio2 & data.cao + data.mgo > 20
    nf(ind) = 100 ./ oxtotal(ind);
else
    disp('oxides')
    % find non-nan loi to normalize oxides
    fprintf('Nomalizing to oxides (default): \n');
    for i = 1:length(oxides)
        fprintf(' %s',oxides{i});
    end
    fprintf('...\n');
    
    % extract oxide data and set negative values to zero for normalization
    oxdata = data{:,oxides};
    oxdata(oxdata < 0) = 0;
    
    % compute total of oxides
    oxtotal = nansum(oxdata,2);
    
    %oxtotal
    fprintf('Oxide Total <40 wt.%%: %i, <90 wt.%%: %i,   >110 wt.%%: %i\n', ...
        sum(oxtotal < 40 & oxtotal ~= 0), ...
        sum(oxtotal < 90 & oxtotal ~= 0), ...
        sum(oxtotal > 110));
    
    %histogram(oxtotal)
    
    ind = ~isnan(data.loi) & oxtotal > 0 & oxtotal < 100;
    
    % compute normalization factor
    nf(ind) = (oxtotal(ind) + data.loi(ind)) ./ oxtotal(ind);
end

subplot(223);
histogram(oxtotal,'BinEdges',[-10:2:110]);
xlabel('Oxide total (wt.%)');
axis square;

subplot(222);
histogram(nf,'BinEdges',[-0.2:0.05:1.2]);
xlabel('Normalization Factor');
axis square;

% compute normalized oxide weights
fprintf('Note: below detection values (BDL) are exluded from oxide total\n');
fprintf('      and normalization factor, but BDL values are normalized.\n');
data.oxtotal = nf.*oxtotal;

data{:,oxides} = repmat(nf,1,length(oxides)).*data{:,oxides};

subplot(224);
histogram(nansum(data{:,oxides},2),'BinEdges',[-10:2:110]);
xlabel('Oxide total (wt.%)');
axis square;

%data.oxtotal = nansum(data{:,oxides},2);

%ind = find(50 < nansum(data{:,oxides},2) & nansum(data{:,oxides},2) < 60);

%data(ind(1:50),:)
%unique(data.filename(ind))

% compute normalized loi
data.loi = nf.*data.loi;

% compute normalized volatiles
data{:,volatiles} = repmat(nf,1,length(volatiles)).*data{:,volatiles};
%[data.loi(1:10) data.oxtotal(1:10)]
% histogram(nf,'BinEdges',[-1:0.02:1])

return

