function data = metamorphic_class(data);

%msed = {{'metasedimentary','metased'},{'arenite'},{'slate'},{'pelite','pelitic'}, ...
%    {'greywacke','graywacke'},{'quartzite','chert'}, ...
%    {'calcsilicate','calc-silicate','marble','dolomite','limestone','carbonate'}, ...
%    {'paraschist'},{'paragneiss'}};

%mig = {{'granite','granitic'},{'basalt','basaltic','basanite'}, ...
%    {'monzonite'},{'diorite','dioritic','tonolite'},{'gabbro','gabbroic'}, ...
%    {'syenite'},{'andesite','andesitic'},{'dacitic'}, ...
%    {'igneous','volcanic','plutonic','mafics','felsics'}, ...
%    {'charonockite','enderbite','norite','mangerite','enderbite', ...
%    'jotunite','farsundite','opdalite'}, ...
%    {'orthogneiss'}, {'dike','dyke','sill'}};

% metamorphic facies list
facies = {{'zeolite'}, ...
    {'hornfels','hornfel','granofel'}, ...
    {'sanidinite'}, ...
    {'phrenite','pumpellyite','phyllit'}, ...
    {'blueschist'}, ...
    {'greenschist'}, ... % will include some amphibolite facies w/ schistose textures
    {'amphibolite','amphibolit','hornblendite'}, ...
    {'granulite','granulit','charnock', ...
    'enderbite','mangerite', ...
    'jotunite','farsundite','opdalite'}, ... % will include some amphibolite facies w/ granulite textures
    {'eclogite','eclogit'}}; % not really a facies, but indicates melting

% metamorphic textures list
texture = {{'slate','slaty'}, ...
    {'schist'}, ...
    {'gneiss'}, ...
    {'migmatite','migmatit'}};

% determine metamorphic facies from names/descriptions
data.rock_name = lower(data.rock_name);
data.rock_facies = lower(data.rock_facies);
for j = 1:length(facies)
    ind = logical(zeros([height(data) 1]));
    for k = 1:length(facies{j})
        tmp1 = strfind(data.rock_name,facies{j}{k});
        tmp2 = strfind(data.rock_facies,facies{j}{k});
        for i = 1:length(tmp1)
            if ~isempty(tmp1{i}) | ~isempty(tmp2{i})
                ind(i) = 1;
            end
        end
    end
    data.facies(ind) = {facies{j}{1}};
end

% determine metamorphic texture from names/descriptions
for j = 1:length(texture)
    ind = logical(zeros([height(data) 1]));
    for k = 1:length(texture{j})
        tmp1 = strfind(lower(data.rock_name),texture{j}{k});
        tmp2 = strfind(lower(data.rock_facies),texture{j}{k});
        for i = 1:length(tmp1)
            if ~isempty(tmp1{i}) | ~isempty(tmp2{i})
                ind(i) = 1;
            end
        end
    end
    data.texture(ind) = {texture{j}{1}};
end

% summary of the facies determinations
fprintf('\n%7s %-7s %-7s %s\n','N','metaig','metased','facies');
fprintf(' ------ ------- ------- ----------------------\n');
cig = 1;
csed = 1;
if any(strcmp(data.Properties.VariableNames,'protolith_est'))
    igind = ( strcmpi('metaigneous',data.rock_origin) | ...
        (strcmpi('igneous',data.protolith_est) & ...
        ~strcmpi('metasedimentary',data.rock_origin)) );
    sedind = ( strcmpi('metasedimentary',data.rock_origin) | ...
        (strcmpi('sedimentary',data.protolith_est) & ...
        ~strcmpi('metaigneous',data.rock_origin)) );
else
    igind = ( strcmpi('metaigneous',data.rock_origin) );
    sedind = ( strcmpi('metasedimentary',data.rock_origin) );
end

for i = 1:length(facies)
    ind = strcmpi(facies{i}{1},data.facies);
    
    igi = ind & igind;
    sedi = ind & sedind;
    fprintf('%7i %7i %7i %s\n',sum(ind),sum(igi),sum(sedi),facies{i}{1});
end

% summary of the texture determinations
fprintf('\n%7s %-7s %-7s %s\n','N','metaig','metased','texture');
fprintf(' ------ ------- ------- ----------------------\n');
for i = 1:length(texture)
    ind = strcmpi(texture{i}{1},data.texture);
    
    igi = ind & igind;
    sedi = ind & sedind;
    fprintf('%7i %7i %7i %s\n',sum(ind),sum(igi),sum(sedi),texture{i}{1});
end

return