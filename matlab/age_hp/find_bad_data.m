spreadsheet = {};
major_elements = {'sio2','tio2','al2o3','cr2o3','feo_tot','mgo','cao','nio',...
    'k2o','na2o','sro','p2o5'};


%Major element searches:
for i = 1:length(major_elements)
    ind = find(data.(major_elements{i}) > 100);
    spreadsheet.(major_elements{i}) = [data.filename(ind,:) data.sample_name(ind,:)];
end


%U, Th searches:
hp_elements = {'u_ppm','th_ppm'};
for i = 1:length(hp_elements)
    ind = find(data.(hp_elements{i}) > 1000);
    spreadsheet.(hp_elements{i}) = [data.filename(ind,:) data.sample_name(ind,:)];
end