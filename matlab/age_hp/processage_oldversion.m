function varargout = processage(data,age_div,varargin)
% processage does all age plots and calculations for paper
%   processage(data,age_div) does plots of data at specified age divs
%   example w preloaded data from hp_process: processage(data,[0:200:4000],'close all')
%
%   data: table format
%   age_div: vector format age divisions


%% 1. Preparation and corrections to raw data
set(0,'DefaultFigureColor',[1 1 1])
%Close all plots and clear the command window
    clc
    %close all

%age_var: maximum age variance of points (e.g. 200Ma range to remove poorly
%defined age rocks such as 'proterozoic' size ranges
    age_var = age_div(end)./(length(age_div)-1);

%age_correction: returns data with age corrections made
%Includes: Removal of ages with greater than age_var variance, swapping min
%and max age where applicable
    data.avg_age = age_correction(data,age_var);
    
%Major and HP element corrections - negatives set to 0 and >100 values 
%divided by 1,000
    data = element_correction(data);
    
%Folder for output of figures
if ismac
    folder = '/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/';
elseif ispc
    folder = 'D:/Users/Matt/Documents/University/Doctor of Philosophy/Background Research/Papers/HP vs Age/figures/';
else
    folder = '';
end
%Set textwidth of LaTeX document
textwidth = 9;

c = 1;
if nargin > 2
    while c <= nargin - 2
        switch lower(varargin{c})
            case 'close all' %Close all
                close all
                c = c+1;
        otherwise
            error(['An unknown option was entered in argument ,', ...
               num2str(6 + c),'.']);
        end
    end   
end


%% 2. Relevant data indices  

%Crude silica range of useful data: Must be interval of 2
    %sio2range = [44 80];
    sio2range = [38 80];
    %sio2range = [0 100];
    
%Specify indices with following conditions:
% - Is not NaN: HP, avg_age
% - Rock must be continental, and of igneous origin
% - Silica must be between the range specified
% - Potassium and Sodium must be greater than 0 for total alkali plots
%ind = data.heat_production > 0 & ~isnan(data.avg_age) ...
%    & ~strcmpi(data.country,'ocean') & ~strcmpi(data.country,'unknown') ...
%    & rockgroup(data,'igneous protolith') & data.sio2 >= sio2range(1) & data.sio2 <= sio2range(2) ...
%    & data.k2o >= 0 & data.na2o >= 0;
ind = data.heat_production > 0 & ~isnan(data.avg_age) ...
    & ~strcmpi(data.country,'ocean') & ~strcmpi(data.country,'unknown') ...
    & rockgroup(data,'all igneous')...
    & data.sio2 >= sio2range(1) & data.sio2 <= sio2range(2);

ind2 = data.heat_production > 0 & ~isnan(data.avg_age) ...
    & ~strcmpi(data.country,'ocean') & ~strcmpi(data.country,'unknown') ...
    & rockgroup(data,'all igneous');%rockgroup(data,'igneous protolith');

data_removed = sum(ind2)-sum(ind);
fprintf('SiO2 correction removes %d data, %f percent of the data\n',data_removed,100*(1-(sum(ind)/sum(ind2))))


figure;
subplot(121);
h = histogram(data.density_model(~isnan(data.heat_production)),'DisplayStyle','stairs');
hold on;
m = mean(data.density_model(~isnan(data.heat_production)));
s = std(data.density_model(~isnan(data.heat_production)));
fprintf('Density Stats:\n');
fprintf('%5i  %4.0f +/- %4.0f\n',sum(~isnan(data.heat_production)),m,s);
plot([m,m],get(gca,'YLim'),'-');
histogram(data.density_model(ind),'DisplayStyle','stairs','BinEdges',h.BinEdges);
m = mean(data.density_model(ind));
s = std(data.density_model(ind));
fprintf('%5i  %4.0f +/- %4.0f\n',sum(ind),m,s);
plot([m,m],get(gca,'YLim'),'--');
golden;
ylabel('No. data');
xlabel('Density [kg m^{-3}]');

subplot(122);
h = histogram(log10(data.heat_production(~isnan(data.heat_production))),'DisplayStyle','stairs');
hold on;
m = mean(log(data.heat_production(~isnan(data.heat_production))));
s = std(log(data.heat_production(~isnan(data.heat_production))));
fprintf('Heat Production Stats:\n');
fprintf('%5i  %4.2f +/- %4.2f\n',sum(~isnan(data.heat_production)),m,s);
plot(log10(exp(1))*[m,m],get(gca,'YLim'),'-');

histogram(log10(data.heat_production(ind)),'DisplayStyle','stairs','BinEdges',h.BinEdges);
m = mean(log(data.heat_production(ind)));
s = std(log(data.heat_production(ind)));
fprintf('%5i  %4.2f +/- %4.2f\n',sum(ind),m,s);
plot(log10(exp(1))*[m,m],get(gca,'YLim'),'--');
hpax([-1 2],'x');
golden;
    
%Number of granites removed:
ind_gran = find(~isnan(data.heat_production) & ~isnan(data.avg_age)...
        & ~strcmpi(data.country,'ocean') & ~strcmpi(data.country,'unknown')...
        & strcmpi(data.rock_type,'granite') ...
        & data.k2o >= 0 & data.na2o >= 0);
ind_gran_sio2 = find(~isnan(data.heat_production) & ~isnan(data.avg_age)...
        & ~strcmpi(data.country,'ocean') & ~strcmpi(data.country,'unknown')...
        & strcmpi(data.rock_type,'granite') & data.sio2 >= sio2range(1)  & data.sio2 <= sio2range(2) ...
        & data.k2o >= 0 & data.na2o >= 0);

fprintf('Percentage of granites removed:%.2f%%\n',100*(1-(size(ind_gran_sio2,1)/size(ind_gran,1))))
    

%Use a reduced temp table with above indices (allowing for further indices
%easier)
temp = data(ind,:);


if ~exist(['supplementary_data_','2017_10_13','.csv'])
    fields = {'sample_id','filename','ref_id','latitude','longitude', ...
        'rock_type','age_min','avg_age','age_max', ...
        'sio2','tio2','al2o3','feo_tot','mgo','cao','na2o','k2o','p2o5', ...
        'th_ppm','u_ppm', ...
        'density_model','heat_production'};
    writetable(temp(temp.heat_production > 0 & temp.avg_age >= 0,fields),['supplementary_data_','2017_10_13','.csv']);
end

%% 3. Figures and Tables

%Table 1: Database statistics & indices for spreadsheets
    numrows = 2;
    %spreadsheet_names: {name,no. data,indices}
    spreadsheet_names = database_table(temp,numrows);
    replace_OZCHEM = find(strcmpi(spreadsheet_names(:,1),'postgres_geochemistry_AUS_2016Aug08.csv'));
    spreadsheet_names{replace_OZCHEM,1} = 'OZCHEM';

%Figure 1: Map
    %Sort the data by age first - guarantees youngest are plotted
    %underneath oldest for cleaner plot
    [~, order] = sort(temp.avg_age);
    sorteddata = temp(order,:);
    
    %Derricks map
    figure;
    plotcoast;
    hold on;
    %fig = mapdata(sorteddata,'color',sorteddata.avg_age,'scale',4*ones(size(sorteddata.heat_production,1),1));
    %set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth/1.6])
    %set(gca,'Units','normalized',...
    %    'FontUnits','points',...
    %    'FontWeight','normal',...
    %    'FontSize',8);
    
    p = plot(data.longitude(~ind),data.latitude(~ind),'ko');
    set(p,'MarkerSize',3,'MarkerFaceColor',[0.15 0.15 0.15],'Color','none');
    scatter(sorteddata.longitude,sorteddata.latitude,3,sorteddata.avg_age);
    
    ylabel('Latitude','FontSize',10)
    xlabel('Longitude','FontSize',10)
    %c = colorbar('Ticks',[0,1000,2000,3000,4000],...
    %     'TickLabels',{'0','1000','2000','3000','4000'});
    %c.Label.String = 'Age [Ma]';
    gtspal('contrast.xlsx');
    %colormap(flipud(load('contrast.csv')));
    %cb = colorbar;
    %cb.YLabel.String = 'Age [Ma]';
    %cb.YTick = [0:400:4000];
    %cb.TickDirection = 'out';
    x1=get(gca,'position');
    %x=get(c,'Position');
    %x(3)=0.02;
    %set(c,'Position',x)
    set(gca,'position',x1)
    fig = gcf;
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])
    %set(gca,'LooseInset',get(gca,'TightInset'))
    
    %Matts map
    hold on;
    fig = mapdata(sorteddata,'color',sorteddata.avg_age,'scale',4*ones(size(sorteddata.heat_production,1),1));
    set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth/1.6])
    set(gca,'Units','normalized',...
       'FontUnits','points',...
       'FontWeight','normal',...
       'FontSize',8);

    ylabel('Latitude','FontSize',10)
    xlabel('Longitude','FontSize',10)
    %c = colorbar('Ticks',[0,1000,2000,3000,4000],...
    %     'TickLabels',{'0','1000','2000','3000','4000'});
    %c.Label.String = 'Age [Ma]';
    gtspal('contrast.xlsx');
    %colormap(flipud(load('contrast.csv')));
    %cb = colorbar;
    %cb.YLabel.String = 'Age [Ma]';
    %cb.YTick = [0:400:4000];
    %cb.TickDirection = 'out';
    x1=get(gca,'position');
    %x=get(c,'Position');
    %x(3)=0.02;
    %set(c,'Position',x)
    set(gca,'position',x1)
    
figure()
    % Sample numbers histogram
%---------------------------
fprintf('Age Range     N\n');
fprintf('---------   -----\n');
for i = 1:length(age_div)-1
    N(i) = sum(age_div(i) <= temp.avg_age & temp.avg_age < age_div(i+1));
    fprintf('%4i %4i  %5i\n',age_div(i),age_div(i+1),N(i));
end
fprintf('-----------------\n');
fprintf('Total:     %6i\n',sum(N));
h = bar(age_div(1:end-1),log10(N),'histc');
set(h,'FaceColor',[0.5 0.5 0.5])
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
hpax([0 4],'y');    
ylabel('N data','FontSize',10);
xlabel('Age [Ma]','FontSize',10);
set(gca,'Box','on')
    
 
    fig = gcf;
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])
    %set(gca,'LooseInset',get(gca,'TightInset'))
    
%    export_fig '/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/map.pdf' -q101


% %Figure 2: Data source distribution
%     %spreadsheet_names holds the names and indices
%     %plot the two largest databases and the rest as misc.
%     fig = figure();
%     database1 = spreadsheet_names(1,:);
%     database2 = spreadsheet_names(2,:);
%     databasemisc = {'Misc.',0,[]};
%     for i = 3:size(spreadsheet_names,1)
%         databasemisc{1,2} = databasemisc{1,2} + spreadsheet_names{i,2};
%         databasemisc{1,3} = union(databasemisc{1,3},spreadsheet_names{i,3});
%     end
%     
%     %age_divs for each set
%     for i = 1:length(age_div)-1
%         indx.d1{i} = find(age_div(i) <= temp.avg_age(database1{1,3},:) & temp.avg_age(database1{1,3},:) < age_div(i+1));
%         indx.d2{i} = find(age_div(i) <= temp.avg_age(database2{1,3},:) & temp.avg_age(database2{1,3},:) < age_div(i+1));
%         indx.dmisc{i} = find(age_div(i) <= temp.avg_age(databasemisc{1,3},:) & temp.avg_age(databasemisc{1,3},:) < age_div(i+1));
%         n.d1(i) = length(indx.d1{i});
%         n.d2(i) = length(indx.d2{i});
%         n.dmisc(i) = length(indx.dmisc{i});
%     end
%     n.d1(n.d1 == 0) = NaN;
%     lognd1 = log10([n.d1';NaN]);
%     
%     n.d2(n.d2 == 0) = NaN;
%     lognd2 = log10([n.d2';NaN]);
%     
%     n.dmisc(n.dmisc == 0) = NaN;
%     logndmisc = log10([n.dmisc';NaN]);
%     
%     h = bar(age_div+100,[lognd1 lognd2 logndmisc],1,'grouped');
%     colormap gray
%     lh=legend(database1{1,1},database2{1,1},databasemisc{1,1})
%     xlim([age_div(1) age_div(end)]);
%     hpax([0 5]);
%     lh.FontSize = 9;
%     set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth/1.6])
%     set(gca,'Units','normalized',...
%         'FontUnits','points',...
%         'FontWeight','normal',...
%         'FontSize',8);
%     ylabel('No. Samples','FontSize',10);
%     xlabel('Age [Ma]','FontSize',10)
%     pos = get(fig,'Position');
%     set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
%     'PaperSize',[pos(3), pos(4)])
%     golden
    %set(gca,'LooseInset',get(gca,'TightInset'))
%    export_fig '/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/datadist.pdf' -q101

% figure;
% z = readtable('Condie&Aster2013_zircon.xlsx');
% for i = 1:length(age_div)-1
%     ind = find(age_div(i) <= temp.avg_age & temp.avg_age < age_div(i+1));
%     rt = rtstats(2,2,data.longitude(ind),data.latitude(ind),log10(data.heat_production(ind)));
%     
%     N(i) = sum(rt.n > 0);
% end
% subplot(211);
% h1 = bar(age_div/1000,[N 0],'histc');
% ylabel('No. Cells');
% xlim([0 4]);
% golden;
% 
% subplot(212);
% h2 = histogram(z.zircon_age/1000,'BinEdges',age_div/1000);
% ylabel('No. Zircons');
% xlabel('Age [Ga]');
% xlim([0 4]);
% golden
% 
% h1.EdgeColor = 'none';
% h2.EdgeColor = 'none';


%Figure 3: Heat production vs Age raw
    [agebin,mafic_ind,felsic_ind] = hpvsage_raw(temp,age_div,folder);
    
%Figure 4: SiO2 correction
    temp.hp_corrected = nan([length(temp.sio2),1]);
    temp.hp_corrected = hp_sio2_2(temp,age_div,sio2range(1):2:sio2range(2),felsic_ind,mafic_ind,folder);
    
    
    %Alkalic vs subalkalic ratio
    alkalicratio(temp,age_div)
    
    
%return

%     temp2 = temp;
%     temp2.heat_production = temp2.hp_corrected;
%     figure;
%     histogram(log10(temp2.hp_corrected));
%     hold on;
%     histogram(log10(temp2.heat_production));
%     hpax([-2 2],'x');
%     close all;
%     
%     ind3 = ~strcmp(temp2.country,'Australia') | ~strcmp(temp2.country,'AU') | ~strcmp(temp2.country,'ocean');
%     agefigures(temp2(ind3,:),age_div)   
% return

%Figure 5: Decay correction

%Do both of these in the same loop - speed up
% 
figure()
    temp.hp_present = nan([length(temp.sio2),1]);
    temp.hp_origin = nan([length(temp.sio2),1]);

    %[temp.hp_present,temp.hp_origin] = hpinitial(temp,age_div,fig,'sio2corrected');
    [temp.hp_present,temp.hp_origin] = hpinitial3(temp,age_div);
    %axis square
%     h = findobj(gca,'Type','patch');
    %[lh,lhicons] = legend(h([2 1]),'Original','Present');
    %lh = legend(h([2 1]),'Original','Present');
    %patchinlegend = findobj(lhicons,'type','patch');
    %set(patchinlegend,'facea',0.3)
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    ylabel('A [\muW m^{-3}]','FontSize',10);
    xlabel('Age [Ma]','FontSize',10)
%     lh.FontSize = 9;
    % Get current position (which I used to keep overall size the same)
%     currentLegendPosition = lh.Position;
    % Define new position
%     newLegendPosition = [0.3235 0.71 currentLegendPosition([3 4])];
    % Set new position
%     lh.Position = newLegendPosition;
    
    textwidth = 16.99757;
    set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth/1.6])
    
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])
%     export_fig '/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/agecorrect.pdf' -q101




%Figure 6: Decay correction without aus
    indaus = (...
        (strcmpi(temp.country,'aus') | strcmpi(temp.country,'australia') | strcmpi(temp.country,'au'))...
        & temp.avg_age > 1200 & temp.avg_age < 2000);
figure()
    %[temp.hp_present(~indaus,:),temp.hp_origin(~indaus,:)] = hpinitial(temp(~indaus,:),age_div,fig,'sio2corrected');
    [temp.hp_present(~indaus,:),temp.hp_origin(~indaus,:)] = hpinitial3(temp(~indaus,:),age_div);
    %axis square
%     h = findobj(gca,'Type','patch');
    %[lh,lhicons] = legend(h([2 1]),'Original','Present');
%     lh = legend(h([2 1]),'Original','Present');
    %patchinlegend = findobj(lhicons,'type','patch');
    %set(patchinlegend,'facea',0.3)
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    ylabel('A [\muW m^{-3}]','FontSize',10);
    xlabel('Age [Ma]','FontSize',10)
%     lh.FontSize = 9;
    % Get current position (which I used to keep overall size the same)
%     currentLegendPosition = lh.Position;
    % Define new position
%     newLegendPosition = [0.774 0.708 currentLegendPosition([3 4])];
    % Set new position
%     lh.Position = newLegendPosition;
    
    textwidth = 16.99757;
    set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth/1.6])
    
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])

    %export_fig '/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/agecorrect_noaus.pdf' -q101
    
    
    
    nout = max(nargout,1);
    if nout == 1
        varargout{1} = temp;
    end
    
% Figure 11: Frequency plots
    fft_hp(temp)
    
    
    
return
%Figure 7: TA vs SiO2
    %perc_hp_TA(temp,age_div)
    
%Figure 8: Jaupart decay
%Talk about this some more before finish it
    jaupartplot

%Figure 9: Rock types
    fields = rock_time(temp,[0:400:4000]);
    
%Figure 10: Comparison other models
    erosionplot(temp)
    textwidth = 16.99757;
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    set(gcf, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth])
    
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])


    
   
return
