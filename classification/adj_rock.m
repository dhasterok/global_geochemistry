function data = adj_rock(data, volc, AdjTAS, Cancrinite, Calcite)

if any(strcmp(data.Properties.VariableNames,'feo_tot'))
    data.feo = data.feo_tot;
end
if ~any(strcmp(data.Properties.VariableNames,'fe2o3'))
    data.fe2o3 = zeros([height(data),1]);
end

%histogram(data.feo,[0:0.1:20]);
%hold on;
%histogram(data.fe2o3,[0:0.1:20]);

oxides = {'sio2','tio2','al2o3','fe2o3','feo','mno','mgo','cao','na2o','k2o','p2o5','co2','loi'};

% Change NaN values 0
for i = 1:length(oxides)
    ind = data{:,oxides{i}} > 0;
    data{~ind,oxides{i}} = 0;
end
 
% If using CO2 to compute minerals, remove CO2 from LOI
if ~(Cancrinite | Calcite)
    data.co2 = zeros([height(data),1]);
end
ind = data.loi > 0;
data.loi(ind) = data.loi(ind) - data.co2(ind);
data.loi(data.loi < 0) = 0;
    
% normalize oxides on anhydrous basis
total = sum(data{:,oxides(1:end-1)},2);

data{:,oxides(1:end-1)} = 100*data{:,oxides(1:end-1)} ./ ...
    repmat(total,1,length(oxides)-1);

% Adjustment of Fe-oxidation (Middlemost, 1989)
%---------------------------------------------------------------------

S = data.sio2;
NK = data.k2o + data.na2o;

% tasfield
%   1: SiO2 verticies
%   2: Na2O + K2O verticies
%   2: Fe ratio
%   3: volcanic name
%   4: plutonic name
tasfield = {[77,100,100,69,69],[0,0,25,25,8],0.5,'rhyolite','granite';
    [57.6,69,30],[11.7,17.73,24.15],0.5,'phonolite','nephelinic syenite';
    [63,69,69,57.6],[7,8,17.73,11.7],0.5,'trachyte/trachydacite','quartz monzonite';
    [63,77,69,63],[0,0,8,7],0.4,'dacite','granodiorite';
    [53,57.6,52.5,48.4],[9.3,11.7,14,11.5],0.4,'tephriphonolite','foid syenite';
    [57,63,57.6,53],[5.9,7,11.7,9.3], 0.4, 'trachyandesite', 'monzonite';
    [57,63,63,57],[0,0,7,5.9], 0.35, 'andesite','diorite';
    [49.4,53,48.4,45],[7.3,9.3,11.5,9.4], 0.35,'phonotephrite','foide syenite';
    [52,57,53,49.4],[5,5.9,9.3,7.3], 0.35,'basaltictrachyandesite','monzodiorite';
    [52,57,57,52],[0,0,5.9,5], 0.3,'basalticandesite','gabrodiorite';
    [41,41,45,49.4,47],[6,7,9.4,7.3,6], 0.3,'tephrite/basanite','foid monzogabbro/diorite';
    [45,52,49.4],[5,5,7.3], 0.3,'trachybasalt','monzogabbro';
    [45,45,52,52],[0,5,5,0], 0.2, 'basalt','gabbro',;
    [41,41,47,45,45],[3,6,6,5,3], 0.2,'tephrite/basanite','foid monzogabbro/diorite';
    [41,45,45,41],[0,0,3,3], 0.15,'picrobasalt','gabbro';
    [46,52.5,30,0,0],[10,14,24.15,24.15,10], 0.4,'foidite','foidite';
    [0,46,41,0],[10,10,7,7], 0.3,'foidite','foidite';
    [0,41,41,0],[7,7,3,3], 0.2,'foidite','foidite';
    [0,41,41,0],[0,0,3,3], 0.15, 'foidite','ijolite'};

[nr,~] = size(tasfield);

c = 1;
data.rock = cell([height(data) 1]);
data.rock(:) = {'unknown'};

data.Fe_ratio = zeros([height(data) 1]);
for i = 1:nr
    ind = inpolygon(S,NK,tasfield{i,1},tasfield{i,2});
    
    data.Fe_ratio(ind) = tasfield{i,3};
    data.rock(volc(ind)) = tasfield(i,4);
    data.rock(~volc(ind)) = tasfield(i,5);
end
% Define type of magma
%---------------------------------------------------------------------

% magmatype
%   1: SiO2 verticies
%   2: Na2O + K2O verticies
%   2: type name
magmatype = {[66,100,100,66],[0,0,24,24], 'Hipersilicic';
    [52,66,66,52],[0,0,24,24], 'Intermediate';
    [45,52,52,45],[0,0,24,24], 'Basic';
    [0,45,45,0],[0,0,24,24], 'Ultrabasic'};

[nm,~] = size(magmatype);

data.magma_type = cell([height(data) 1]);
data.magma_type(:) = {''};
for i = 1:nm
    ind = inpolygon(S,NK,magmatype{i,1},magmatype{i,2});
    
    data.magma_type(ind) = {magmatype{i,3}};
end



% Define FeOx output
%---------------------------------------------------------------------
if ~AdjTAS
    data.Fe_ratio = oxide_fe(volc,data.sio2,data.na2o,data.k2o);
end

ind = data.fe2o3 == 0 | data.feo == 0;
data.Fe_ratio(ind) = data.Fe_ratio(ind);

data.Fe_ratio(~ind) = data.fe2o3(~ind)./data.feo(~ind);


% Define Fe.output (Pruseth, 2009)
%---------------------------------------------------------------------
Rfe = 2*molecularwt('FeO')/molecularwt('Fe2O3');
feot = data.feo + Rfe * data.fe2o3;
%histogram(data.Fe_ratio,'BinWidth',0.001);
data.feo =  feot./(1 + Rfe * data.Fe_ratio);
data.fe2o3 = (feot - data.feo)/Rfe;


% Recalculate major elements data on an anhydrous basis
%---------------------------------------------------------------------
total = sum(data{:,oxides(1:end-1)},2);

data{:,oxides(1:end-1)} = 100*data{:,oxides(1:end-1)} ./ ...
    repmat(total,1,length(oxides)-1);

return


function Fe_ratio = oxide_fe(volc,sio2,na2o,k2o)
% Adjustment of Fe-oxidation (LeMaitre, 1976)
%---------------------------------------------------------------------

% volcanics
Fe_ratio(volc) = 0.93 - 0.0042 * sio2(volc) - ...
    0.022*(na2o(volc) + k2o(volc));
% plutonics
Fe_ratio(~volc) = 0.88 - 0.0016 * sio2(~volc) - ...
    0.027*(na2o(~volc) + k2o(~volc));

Fe_ratio = Fe_ratio(:);

return
