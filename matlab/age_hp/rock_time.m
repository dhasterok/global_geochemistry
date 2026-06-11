function fields = rock_time(data,age_div)

Exist_Column = strcmp('avg_age',data.Properties.VariableNames);
val = Exist_Column(Exist_Column==1);
if isempty(val)
    data.avg_age = data.age;
end


fields = {
    '', 'peridotite',[0];
    'picrobasalt', 'peridotgabbro',[0];
    'alkalic basalt', 'alkalic gabbro',[0];
    'subalkalic basalt', 'subalkalic gabbro',[0];
    'basaltic andesite', 'gabbroic diorite',[0];
    'andesite', 'diorite',[0];
    'dacite', 'granodiorite',[0];
    'rhyolite', 'granite',[0];
    'trachybasalt', 'monzogabbro',[0];
    'basaltic trachyandesite', 'monzodiorite',[0];
    'trachyandesite', 'monzonite',[0];
    'trachydacite', 'quartz monzonite',[0];
    'trachyte', 'syenite',[0];
    'tephrite', 'foid gabbro',[0];
    'phonotephrite', 'foid monzodiorite',[0];
    'tephriphonolite', 'foid monzosyenite',[0];
    'phonolite', 'foid syenite',[0];
    'ultramafic foidite','ultramafic foidolite',[0];
    'mafic foidite', 'mafic foidolite',[0];
    'intermediate foidite', 'intermediate foidolite',[0];
    'silexite', 'quartzolite',[0];
    'ultra-high alkali volcanic','ultra-high alkali plutonic',[0];
    'picrite', '',[0];
    'alkali picrite', '',[0];
    'komatiite', 'meimechite',[0];
    'boninite', '',[0]};

for i = 1:size(fields,1)
    for j = 1:2
        if strcmpi('',fields{i,j})
            continue;
        else
            ind = find(strcmpi(data.rock_type,fields{i,j}));
            fields{i,3} = fields{i,3} + length(ind);
        end
    end
end

[trash idx] = sort([fields{:,3}], 'descend');
temp = fields(idx,:);

%Top 3 fields
fields2plot = temp(1:6,:)

figNum = 1;
figure(figNum)
textwidth = 16.99757;
set(1, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth/1.6])

temp = data;


colorsets = {[0 0 0],[1 0 1],[0 1 1],[1 0 0],[0 1 0],[0 0 1]};


for i = 1:3
    data = temp;
    ind2 = [];
    ind2 = find(strcmpi(data.rock_type,fields2plot{i,1}));
    ind2 = union(ind2,find(strcmpi(data.rock_type,fields2plot{i,2})));
    data = data(ind2,:);
    for j = 1:length(age_div)-1
        ind = (age_div(j) <= data.avg_age & data.avg_age < age_div(j+1));
        agebin{j} = ind;
        avg_age{j} = data.avg_age(ind);
        heat_production{j} = data.heat_production(ind);
    end
    
    subplot(3,1,i)
    [Qage,Qhp] = whisker(avg_age,heat_production,'Color',[0.5 0.5 0.5],'Scale','log');
    plot(0,0)
    sqwavefill(Qhp,Qage(:,3),age_div,colorsets{i})
    xlim([age_div(1) age_div(end)]);
    strtitle = strcat(fields2plot{i,1},'/\n',fields2plot{i,2});
    
    strtitle = [fields2plot{i,1} sprintf('/\n%s',fields2plot{i,2})];
    
    
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);%,...
        %'FontName', 'Times');
    
    title(strtitle);
    %set(gca,'Box','on');
    %hpax([floor(min(Qhp(:,1))) ceil(max(Qhp(:,5)))]);
    hpax([-2 2]);
    axis normal
    set(gca,'Box','off');
    ylabel('')

    if i== 4
        ylabel('A [\muW m^{-3}]','FontSize',10);
    end
    
    if i == 5
        xlabel('Age [Ma]','FontSize',10)
    end
    
    x1 = 2200;
    y1 = 1.7;
    txt1 = horzcat('n = ',int2str(fields2plot{i,3}));
    text(x1,y1,txt1)
    %set(gca,'ticklength',3*get(gca,'ticklength'))
    
    a = gca;
    a.XAxis.Visible = 'off';
    if i == 6
        a.XAxis.Visible = 'on';
    end
    
end


pos = get(figNum,'Position');
set(figNum,'PaperPositionMode','Auto','PaperUnits','centimeters',...
'PaperSize',[pos(3), pos(4)])
%export_fig '/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/rock_types.pdf' -q101

figNum = 2;
figure(figNum)
textwidth = 16.99757;
set(1, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth/1.6])
for i = 4:6
    data = temp;
    ind2 = [];
    ind2 = find(strcmpi(data.rock_type,fields2plot{i,1}));
    ind2 = union(ind2,find(strcmpi(data.rock_type,fields2plot{i,2})));
    data = data(ind2,:);
    for j = 1:length(age_div)-1
        ind = (age_div(j) <= data.avg_age & data.avg_age < age_div(j+1));
        agebin{j} = ind;
        avg_age{j} = data.avg_age(ind);
        heat_production{j} = data.heat_production(ind);
    end
    
    subplot(3,1,i-3)
    [Qage,Qhp] = whisker(avg_age,heat_production,'Color',[0.5 0.5 0.5],'Scale','log');
    plot(0,0)
    sqwavefill(Qhp,Qage(:,3),age_div,colorsets{i})
    xlim([age_div(1) age_div(end)]);
    strtitle = strcat(fields2plot{i,1},'/\n',fields2plot{i,2});
    
    strtitle = [fields2plot{i,1} sprintf('/\n%s',fields2plot{i,2})];
    
    
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);%,...
        %'FontName', 'Times');
    
    title(strtitle);
    %set(gca,'Box','on');
    %hpax([floor(min(Qhp(:,1))) ceil(max(Qhp(:,5)))]);
    hpax([-2 2]);
    axis normal
    set(gca,'Box','off');
    ylabel('')

    if i== 4
        ylabel('A [\muW m^{-3}]','FontSize',10);
    end
    
    if i == 5
        xlabel('Age [Ma]','FontSize',10)
    end
    
    x1 = 2200;
    y1 = 1.7;
    txt1 = horzcat('n = ',int2str(fields2plot{i,3}));
    text(x1,y1,txt1)
    %set(gca,'ticklength',3*get(gca,'ticklength'))
    
    a = gca;
    a.XAxis.Visible = 'off';
    if i == 6
        a.XAxis.Visible = 'on';
    end
    
end
pos = get(figNum,'Position');
set(figNum,'PaperPositionMode','Auto','PaperUnits','centimeters',...
'PaperSize',[pos(3), pos(4)])
%export_fig '/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/rock_types.pdf' -q101




return