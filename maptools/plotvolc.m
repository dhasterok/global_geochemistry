function varargout = plotvolc(varargin)

size = 12;
if nargin == 1
    size = varargin{1};
end

% PLOTVOLC - plots locations of holocene and quaternary volcnoes
volc.holo = readtable('../GIS/volcanoes/GVP_Volcano_List_Holocene.csv', ...
    'Format','%f%q%q%q%q%q%q%q%f%f%f%q%q%q');
volc.pleist = readtable('../GIS/volcanoes/GVP_Volcano_List_Pleistocene.csv', ...
    'Format','%f%q%q%q%q%q%q%q%f%f%f%q%q%q');
volc.possible = readtable('../GIS/volcanoes/possible_volcanic_centers.csv', ...
    'Format','%f%q%q%q%q%q%q%q%f%f%f%q%q%q');

if nargout == 1
    volc.holo.age = cell([height(volc.holo) 1]);
    volc.holo.age(:) = {'Holocene'};
    volc.pleist.age = cell([height(volc.pleist) 1]);
    volc.pleist.age(:) = {'Pleistocene'};
    volc.possible.age = cell([height(volc.possible) 1]);
    volc.possible.age(:) = {'Unknown'};
    
    varargout{1} = [volc.holo; volc.pleist; volc.possible];
end

if isempty(size)
    return
end

C = [0    0.4470    0.7410;
   0.8500    0.3250    0.0980;
   0.9290    0.6940    0.1250;
   0.4940    0.1840    0.5560;
   0.4660    0.6740    0.1880;
   0.6350    0.0780    0.1840;
   0.3010    0.7450    0.9330];

setting = unique(cat(1,cat(1,volc.holo.TectonicSetting,volc.possible.TectonicSetting),volc.pleist.TectonicSetting));
for i = 1:length(setting)
    ind = strcmp(volc.pleist.TectonicSetting,setting{i});
    
    s(i) = scatter(volc.pleist.Longitude(ind),volc.pleist.Latitude(ind),'^');
    set(s(i),'SizeData',size,'CData',C(i,:));
    hold on;
    
    ind = strcmp(volc.holo.TectonicSetting,setting{i});
    s2(i) = scatter(volc.holo.Longitude(ind),volc.holo.Latitude(ind),'^','filled');
    set(s2(i),'SizeData',size,'CData',C(i,:));
    
    ind = strcmp(volc.possible.TectonicSetting,setting{i});
    s3(i) = scatter(volc.possible.Longitude(ind),volc.possible.Latitude(ind),'^');
    set(s3(i),'SizeData',size,'CData',C(i,:));
end

legend(s,setting,'Location','north','Orientation','horizontal');

return