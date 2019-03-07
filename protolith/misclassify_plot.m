function misclassify_plot(data,protolith_est,ind)

figure;
subplot(121);
%misclassify_sub(data,protolith_est,'igneous','sedimentary');
misclassify_sub2(data,protolith_est,'igneous','sedimentary');

subplot(122);
%misclassify_sub(data,protolith_est,'sedimentary','igneous');
misclassify_sub2(data,protolith_est,'sedimentary','igneous');


return

function misclassify_sub(data,protolith_est,truetype,prototype)

fn_ind = rockgroup(data,truetype) & strcmp(prototype,protolith_est);
tp_ind = rockgroup(data,truetype) & strcmp(truetype,protolith_est);

rt = unique(data.rock_type(fn_ind));
n = zeros(size(rt));
for i = 1:length(rt)
    for j = 1:length(rt{i})
        n(i) = n(i) + nansum(strcmp(data.rock_type(fn_ind),rt{i}{j}))
        nrt(i) = n(i) + nansum(strcmp(data.rock_type(tp_ind),rt{i}{j}));
    end
    [j n(i)]
    p(i) = n(i)/nrt(i);
end
bar(n,'stacked')
set(gca,'XTick',[0:length(rt)+1],'XTickLabel',{'',rt{:},''});
for i = 1:length(p);
    t = text(i,n(i)+0.05*max(n),num2str(100*p(i),2));
    set(t,'HorizontalAlign','right');
end

xlim([0 length(rt)+1]);
ylabel(['No. misclassified as ',prototype]);
view([90 90]);
return

function misclassify_sub2(data,protolith_est,truetype,prototype)

fn_ind = rockgroup(data,truetype) & strcmp(prototype,protolith_est);
tp_ind = rockgroup(data,truetype) & strcmp(truetype,protolith_est);

rt = getrtypes;
n = zeros(size(rt));
m = zeros(size(rt));
for i = 1:length(rt)
    str{i} = rt{i}{1};
    for j = 1:length(rt{i})
        %[j, n(i), sum(strcmp(data.rock_type(fn_ind),rt{i}{j}))]
        n(i) = n(i) + sum(strcmp(data.rock_type(fn_ind),rt{i}{j}));
        m(i) = m(i) + sum(strcmp(data.rock_type(tp_ind),rt{i}{j}));
        if j == 2
            str{i} = [str{i},'/',rt{i}{2}];
        end
    end
    p(i) = n(i)/(n(i) + m(i));
end
map = parula(101);
hold on;
for i = 1:length(n)
    cind = round(10*sqrt(100*p(i)));
    h = bar(i,n(i));
    if ~isnan(cind)
        h.FaceColor = map(cind+1,:);
    end
end
colormap(parula);
caxis([0 10]);
cbar;
set(gca,'Box','on','XTick',[0:length(str)+1],'XTickLabel',{'',str{:},''});
for i = 1:length(p);
    if ~isnan(n(i))
        t = text(i,n(i)+0.05*max(n),num2str(100*p(i),2));
    end
    set(t,'HorizontalAlign','right');
end

xlim([0 length(str)+1]);
ylabel(['No. misclassified as ',prototype]);
view([90 90]);
return

function rt = getrtypes

rt = {{'silexite', 'quartzolite'};
    {'rhyolite', 'granite'};
    {'dacite', 'granodiorite'};
    {'andesite', 'diorite'};
    {'basaltic andesite', 'gabbroic diorite'};
    {'subalkalic basalt', 'subalkalic gabbro'};
    {'picrobasalt', 'peridotgabbro'};
    {'cumulate peridotite','crustal peridotite','crustal pyroxenite','peridotite'};
    {'trachyte', 'syenite'};
    {'trachydacite', 'quartz monzonite'};
    {'trachyandesite', 'monzonite'};
    {'basaltic trachyandesite', 'monzodiorite'};
    {'trachybasalt', 'monzogabbro'};
    {'alkalic basalt', 'alkalic gabbro'};
    {'phonolite', 'foid syenite'};
    {'tephriphonolite', 'foid monzosyenite'};
    {'phonotephrite', 'foid monzodiorite'};
    {'tephrite', 'foid gabbro'};
    {'ultra-high alkali volcanic','ultra-high alkali plutonic'};
    {'foidite','foidolite','intermediate foidite', 'intermediate foidolite','mafic foidite', 'mafic foidolite','ultramafic foidite','ultramafic foidolite'};
    {'boninite', 'sanukitoid'};
    {'picrite','alkali picrite','ferropicrite','alkali ferropicrite'};
    {'komatiite', 'meimechite', 'basaltic komatiite','intrusive gabbroic komatiite'};
    {'mantle peridotite', 'mantle pyroxenite'};
    {'carbonatite','calciocarbonatite','ferrocarbonatite','magnesiocarbonatite'};
    {'silicocarbonatite','silico-calciocarbonatite','silico-ferrocarbonatite','silico-magnesiocarbonatite'};
    {'quartzite'};
    {'quartz arenite'};
    {'litharenite'};
    {'sublitharenite'};
    {'arkose'};
    {'subarkose'};
    {'wacke'};
    {'shale'};
    {'iron-rich shale'};
    {'iron-rich sand'};
    {'laterite','bauxite','oxide'};
    {'limestone'};
    {'dolomite'}};
    
return