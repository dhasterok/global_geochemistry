function [data,provinces] = addprovinces(data)

%cont = {'africa', 'antarctica', 'asia', 'australia', 'europe', 'india', ...
%    'north_america', 'south_america', 'lips', 'plates'};
cont = {'global','lips','plates'};

% add polar coordinates to database (for Antarctica)
[y,x] = pol2cart(data.longitude*pi/180,(data.latitude+90)*111.2);
%data = addvars(data,x,y,'After','longitude','NewVariableNames',{'polarx','polary'});

ind = data.latitude > -40;
data.polarx(ind) = NaN;
data.polary(ind) = NaN;

%figure;
%plotcoast;
%hold on;
newfields = {'plate_major','poly_name','plate_subplate','ocean_domain', ...
    'ocean_name','plate_type','crust_type','prov_name','prov_type', ...
    'prov_group','last_orogen','lip_name'};

for i = 1:length(newfields)
    if any(strcmp(data.Properties.VariableNames,newfields{i}))
        data(:,newfields{i}) = [];
    end
    data(:,newfields{i}) = {''};
end
% data.plate_major(:) = {''};
% data.poly_name(:) = {''};
% data.plate_subplate(:) = {''};
% data.ocean_domain(:) = {''};
% data.ocean_name(:) = {''};
% data.plate_type(:) = {''};
% data.crust_type(:) = {''};
% data.prov_name(:) = {''};
% data.prov_type(:) = {''};
% data.prov_group(:) = {''};
% data.last_orogen(:) = {''};
% data.lip_name(:) = {''};

for i = 1:length(cont)
    fprintf('Working on %s...\n',cont{i});
    switch cont{i}
        case 'antarctica'
            filename = ['../GIS/global_tectonics/plates&provinces/',cont{i},'_gprv.shp'];
        case 'lips'
            filename = ['../GIS/global_tectonics/polygon_data/Johansson_etal_2018_EarthByte_LIPs_v2.shp'];
        case 'plates'
            filename = ['../GIS/global_tectonics/plates&provinces/plates.shp'];
        otherwise
            filename = ['../GIS/global_tectonics/plates&provinces/',cont{i},'_gprv.shp'];
    end
    
    s = struct2table(shaperead(filename));
    if i == 1 & ~strcmp(cont{i},'lips') & ~strcmp(cont{i},'plates')
        provinces = s;
    elseif ~strcmp(cont{i},'lips') & ~strcmp(cont{i},'plates')
        provinces = [provinces; s];
    end
    %s.Properties.VariableNames
    for j = 1:height(s)
        switch cont{i}
            case 'antarctica'
                in = inpolygon(data.polarx,data.polary,s.X{j}'*1e-3,s.Y{j}'*1e-3);
            case 'lips'
                in = inpolygon(data.longitude,data.latitude,s.X{j}',s.Y{j}');
                if sum(in) > 0
                    data(in,'lip_name') = {s.NAME{j}};
                end
                continue;
            case 'plates'
                in = inpolygon(data.longitude,data.latitude,s.X{j}',s.Y{j}');
                if sum(in) > 0
                    data(in,'plate_major') = {s.plate{j}};
                    data(in,'poly_name') = {s.poly_name{j}};
                    data(in,'plate_subplate') = {s.subplate{j}};
                    data(in,'ocean_domain') = {s.domain{j}};
                    data(in,'ocean_name') = {s.sea_name{j}};
                    data(in,'plate_type') = {s.plate_type{j}};
                    data(in,'crust_type') = {s.crust_type{j}};
                end
                continue;
            otherwise
                in = inpolygon(data.longitude,data.latitude,s.X{j}',s.Y{j}');
        end
        %fprintf('    %s, found %i data...\n',s(j).prov_name,sum(in));
        if sum(in) > 0
            data(in,'prov_name') = {s.prov_name{j}};
            data(in,'prov_type') = {s.prov_type{j}};
            data(in,'prov_group') = {s.prov_group{j}};
            data(in,'last_orogen') = {s.lastorogen{j}};
        end
        %scatter(data.longitude(in),data.latitude(in),5,'o','filled');
    end
    
    %s.cont = repmat({cont{i}},height(s),1);
    %if i == 1
    %    provinces = s;
    %    vn = provinces.Properties.VariableNames;
    %else
    %    provinces = [provinces; s(:,vn)];
    %end
end

return