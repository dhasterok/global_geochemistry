function data = adjust_origin(data);

mig = {'ingeous','orthogneiss','orthoamphibolite'};

mp = {'granit','grano','diorit','tonolit','TTG','gabbro', ...
   'monzonit','syenit', ...
   'plutonic', ...
   'charonockit','enderbite','norit','mangerite','enderbite', ...
   'jotunite','farsundite','opdalite', ...
   'dike','dyke','sill'};

mv = {'basalt','basanit','rhyolit','dacit','andesit','trachy', ...
   'volcanic','lava','tuff'};
   
ms = {'sediment','conglomerat','carbon','graphit','slate','shale', ...
    'mud','silt','sand','lime','marble','dolomit','para'};

mind = rockgroup(data,'metamorphic') & ~rockgroup(data,'metaigneous') & ~rockgroup(data,'metasedimentary');
ind = logical(zeros([height(data) 1]));
for j = 1:length(mig)
    tmp1 = strfind(data.rock_name,mig{j});
    for i = 1:length(tmp1)
        if ~isempty(tmp1{i})
            ind(i) = 1;
        end
    end
end
ind = ind & mind;
data.rock_origin(ind) = {'metaigneous'};
fprintf('Changed %i rock origins to metaigneous.\n',sum(ind));


mind = rockgroup(data,'metamorphic') & ~rockgroup(data,'metaplutonic') & ~rockgroup(data,'metavolcanic') & ~rockgroup(data,'metasedimentary');
ind = logical(zeros([height(data) 1]));
for j = 1:length(mp)
    tmp1 = strfind(data.rock_name,mp{j});
    for i = 1:length(tmp1)
        if ~isempty(tmp1{i})
            ind(i) = 1;
        end
    end
end
ind = ind & mind;
data.rock_origin(ind) = {'metaplutonic'};
fprintf('Changed %i rock origins to metaplutonic.\n',sum(ind));

mind = rockgroup(data,'metamorphic') & ~rockgroup(data,'metaplutonic') & ~rockgroup(data,'metavolcanic') & ~rockgroup(data,'metasedimentary');
for j = 1:length(mv)
    tmp1 = strfind(data.rock_name,mv{j});
    for i = 1:length(tmp1)
        if ~isempty(tmp1{i})
            ind(i) = 1;
        end
    end
end
ind = ind & mind;
data.rock_origin(ind) = {'metavolcanic'};
fprintf('Changed %i rock origins to metavolcanic.\n',sum(ind));

mind = rockgroup(data,'metamorphic') & ~rockgroup(data,'metaplutonic') & ~rockgroup(data,'metavolcanic') & ~rockgroup(data,'metasedimentary');
for j = 1:length(ms)
    tmp1 = strfind(data.rock_name,ms{j});
    for i = 1:length(tmp1)
        if ~isempty(tmp1{i})
            ind(i) = 1;
        end
    end
end
ind = ind & mind;
data.rock_origin(ind) = {'metasedimentary'};
fprintf('Changed %i rock origins to metasedimentary.\n',sum(ind));

return