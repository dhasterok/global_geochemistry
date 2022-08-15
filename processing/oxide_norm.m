function data = oxide_norm(data,varargin)
% OXIDE_NORM - Computes the oxide norm for a set of given oxides
%
%   data = oxide_norm(data,oxides) computes the oxide norm for the oxides in the
%   given cell array 'oxides'.  The normalization of iron is set by the list 'oxides'.
%
%   Normalization is performed to total iron converted to FeO.
%
%   The following option value pairs may be given:
%
%       'Normalization'     'anhydrous' (default) or 'hydrous',
%                           will normalize using a set of major oxides, 
%                           oxides can be set using 'Oxides'.  Note
%                           selecting hydrous will significantly reduce the
%                           size of the database.
%
%       'Oxides'            entered as a cell array of oxides.  The default
%                           list is:
%                           {'SiO2', 'TiO2', 'Al2O3', 'Cr2O3',
%                           'FeO', 'MnO', 'MgO', 'NiO', 'CaO', 
%                           'Na2O', 'K2O', 'P2O5'}
%
%       'TotalTol'          tolerance for sum of oxides, i.e., sum of
%                           majors + loi/(h2o + co2 + so2 + f + cl) must be
%                           within 100 +/- tol in wt%.  Default is 10.


% ---------------------------------
% Oxide lists
% ---------------------------------
% all avaliable oxides
oxides = {'sio2','tio2','al2o3','feo_tot', ...
    'mgo','cao','na2o','k2o','p2o5', ... % common oxides
    'cr2o3','mno','nio'}; % less common oxides

% all volatiles
volatiles = {'h2o_tot','co2','so3','f_ppm','cl_ppm'};

% conversion factor for volatiles
cf = [1 1 1 1e-4 1e-4];

tau = 100;

% ---------------------------------
% parse inputs
% ---------------------------------
p = inputParser;

addRequired(p,'data');
addParameter(p,'Oxides',oxides,@iscell);
addParameter(p,'Normalization','anhydrous',@ischar);
addParameter(p,'TotalTol',tau,@isnumeric);

parse(p,data,varargin{:});

oxides = lower(p.Results.Oxides);
tau = p.Results.TotalTol;

% remove volatiles and LOI from list of oxides
oxides = setdiff(oxides,volatiles);

for i = 1:length(oxides)
    if ~any(strcmp(data.Properties.VariableNames,oxides{i}))
        data{:,oxides{i}} = nan([height(data) 1]);
    end
end
if ~any(strcmp(data.Properties.VariableNames,'loi'))
    data.loi = nan([height(data) 1]);
end

% fixes norm oxide list to lowercase and feo to feo_tot
if any(strcmp('feo',oxides))
    oxides(strcmp('feo',oxides)) = {'feo_tot'};
end


% ---------------------------------
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

data.total_loi = false([height(data) 1]);
data.total_loi(ind) = 1;

%fprintf('LOI <10 wt.%%: %i,   >10 wt.%%: %i\n',sum(data.loi < 10),sum(data.loi > 10));

%subplot(221);
%histogram(data.loi,'BinEdges',[-10:2:110]);
%xlabel('LOI (wt.%)');
%axis square;

fprintf('Normalizing to %s conditions...\n',p.Results.Normalization);
fprintf('Nomalizing to oxides: \n');
for i = 1:length(oxides)
    fprintf(' %s',oxides{i});
end

% extract oxide data and set negative values to zero for normalization
oxdata = data{:,oxides};
oxdata(oxdata < 0) = 0;

% compute total of oxides
data.total_ox = nansum(oxdata,2);

% remove totals outsize acceptable range of tolerance
% note will keep all samples without majors.
if tau ~= 100
    ind = (100-tau < data.total_ox & data.total_ox < 100+tau) | ... % anhydrous total within bounds
        (100-tau < data.total_ox + data.total_loi & ...             % hydrous total within bounds
        data.total_ox + data.total_loi < 100+tau) | ...
        isnan(data.total_ox);                                       % no oxide data
    data = data(ind,:);
    
    fprintf('Tolerance filter of %i applied, %i/%i samples removed...\n',tau,size(oxdata,1),height(data));
else
    fprintf('No tolerance filter...\n');
end

% determine normalization factor
nf = ones([height(data) 1]);
switch p.Results.Normalization
    case 'anhydrous'
        fprintf('\n');
        % compute normalization factor
        %if data.cao + data.mgo > data.sio2 & data.cao + data.mgo > 20
        ind = data.total_ox > 0;
        nf(ind) = 100 ./ data.total_ox(ind);
    case 'hydrous'
        fprintf(' + LOI/volatiles\n');
        
        ind = data.total_ox > 0 & data.loi > 0;
        nf(ind) = 100 ./ (data.total_ox(ind) + data.loi(ind));
    otherwise
        error('Unknown normalization type.');
end

%figure;
%histogram(nf,'BinWidth',0.01)

% compute normalized oxide weights
fprintf('Note: below detection values (BDL) are exluded from oxide total\n');
fprintf('      and normalization factor, but BDL values are normalized.\n');
data{:,oxides} = repmat(nf,1,length(oxides)).*data{:,oxides};

data.total_ox = nf.*data.total_ox;
data.total_loi = nf.*data.total_loi;

% compute normalized loi
data.loi = nf.*data.loi;

% compute normalized volatiles
data{:,volatiles} = repmat(nf,1,length(volatiles)).*data{:,volatiles};

fprintf('Oxide Total <40 wt.%%: %i, <90 wt.%%: %i,   >110 wt.%%: %i\n', ...
    sum(data.total_ox < 40 & data.total_ox ~= 0), ...
    sum(data.total_ox < 90 & data.total_ox ~= 0), ...
    sum(data.total_ox > 110));
 

%subplot(223);
%histogram(oxtotal,'BinEdges',[-10:2:110]);
%xlabel('Oxide total (wt.%)');
%axis square;

%subplot(222);
%histogram(nf,'BinEdges',[-0.2:0.05:1.2]);
%xlabel('Normalization Factor');
%axis square;

%subplot(224);
%histogram(nansum(data{:,oxides},2),'BinEdges',[-10:2:110]);
%xlabel('Oxide total (wt.%)');
%axis square;

%data.oxtotal = nansum(data{:,oxides},2);

%ind = find(50 < nansum(data{:,oxides},2) & nansum(data{:,oxides},2) < 60);

%data(ind(1:50),:)
%unique(data.filename(ind))


%[data.loi(1:10) data.oxtotal(1:10)]
% histogram(nf,'BinEdges',[-1:0.02:1])

return

