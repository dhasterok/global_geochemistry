%Data already loaded - data, plutonic, sedimentary, volcanic, metamorphic

function processage_draft(data,age_div)
clc
close all

%Last thing to change: Colours, scale, labels, titles


%processage(data,age_div,name)
%Most common: processage(data,[0:200:4000])
%age_div: vector containing age ranges e.g. [0:200:4000]
%data: table format

%Recalculating values:
%https://australianmuseum.net.au/classification-of-igneous-rocks
%TAS classification scheme renormalizes without H2o and co2
%Therefore - use that

%Age correction and maximum age variance between min and max
%Set age_var to the be the average bin size for the age_divs - point will
%be at most one bin in error then on average
    age_var = age_div(end)./(length(age_div)-1);
    data.avg_age = age_correction(data,age_var);
    
%Major and HP element corrections - negatives set to 0 and >100 values 
%divided by 1,000
    data = element_correction(data);
    
%Restrict data:
%Depending on function:
%   1. sio2 between 30 and 90
%   2. has HP estimate
%   3. has an age or time period
%   4. Continental/country location
% ind = find(data.sio2 > 30 & data.sio2 < 90 & ...
%         ~isnan(data.heat_production) & ~isnan(data.avg_age)...
%         & ~strcmpi(data.country,'ocean') & ~strcmpi(data.country,'unknown'));
%     data = data(ind,:);
    
%Remove Australia option: Between 1600 and 2000, anomalous population
%PLOT WORLD VS WORLD WITHOUT AUS TO PICK RANGE in 100Mya bins?
%     remove_aus = 1;
%     if remove_aus
%         ind = find((strcmpi(data.country,'au') | strcmpi(data.country,'aus') | ...
%             strcmpi(data.country,'australia')) & data.avg_age < 2000 & data.avg_age > 1400);
%         ind = (strcmpi(data.country,'au') | strcmpi(data.country,'aus') | ...
%             strcmpi(data.country,'australia')) & data.avg_age < 2000 & data.avg_age > 1400;
%         %fprintf('%i Australian Proterozoic data removed\n',length(find(ind==1)))
%         data = data(~ind,:);
%     end

%Indices that persist through all: Must be intervals of 2
    sio2range = [44 78];
    
    ind = find(~isnan(data.heat_production) & ~isnan(data.avg_age)...
        & ~strcmpi(data.country,'ocean') & ~strcmpi(data.country,'unknown')...
        & rockgroup(data,'igneous protolith') & data.sio2 >= sio2range(1)  & data.sio2 <= sio2range(2) ...
        & data.k2o >= 0 & data.na2o >= 0);
    
    remove_aus = 1;
    if remove_aus
        indaus = (strcmpi(data.country(ind,:),'au') | strcmpi(data.country(ind,:),'aus') | ...
            strcmpi(data.country(ind,:),'australia')) & data.avg_age(ind,:) < 2000 & data.avg_age(ind,:) > 1400;
        fprintf('%i Australian Proterozoic data removed\n',length(find((indaus)==1)))
        ind2 = ind(~indaus);
        ind = intersect(ind,ind2);
    end

    
%Figure 14: 6 most common rock types
    %Restricted age div so sqwavefill works (doesnt deal with gaps yet)
    fields = rock_time(data(ind,:),[0:400:4000]);

%Figure 1: Map
%     ind = find(~isnan(data.heat_production) & ~isnan(data.avg_age)...
%         & ~strcmpi(data.country,'ocean') & ~strcmpi(data.country,'unknown')...
%         & rockgroup(data,'igneous protolith'));
    fig = mapdata(data(ind,:),'color',log10(data.heat_production(ind,:)),'scale',3.*ones(size(data.heat_production(ind,:),1),1))
    %fig = mapdata(data(ind,:),'color',data.avg_age(ind,:),'scale',10.*ones(size(data.heat_production(ind,:),1),1));
    textwidth = 16.99757;
    %set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth/2])
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    ylabel('Latitude','FontSize',10)
    xlabel('Longitude','FontSize',10)
    
    c = colorbar('Ticks',[-2,-1,0,1,2],...
         'TickLabels',{'0.01','0.1','0','10','100'});
    c.Label.String = 'Heat Production [\muW m^{-3}]';
    x1=get(gca,'position');
    x=get(c,'Position');
    x(3)=0.02;
    set(c,'Position',x)
    set(gca,'position',x1)
    
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])
    %pdf
%    print(fig,['/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/' ...
%        'fig' int2str(2)],'-dpdf','-r0')
    
%Table 1: Database statistics
% ind = (~isnan(data.heat_production) & ~isnan(data.avg_age)...
%     & ~strcmpi(data.country,'ocean') & ~strcmpi(data.country,'unknown')...
%     & rockgroup(data,'igneous protolith'));
%Num rows: Specifies how many entries get listed in table before being
%counted only in 'Misc.'
numrows = 2;
spreadsheet_names = database_table(data(ind,:),numrows);


%Figures 2 & 3: Database Statistics
    fig = figure();
    subplot(1,2,1)
    histogram(data.sio2(ind,:),'BinWidth',1)
    xlim([30 90])
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    xlabel('SiO_2 [wt%]','FontSize',10)
    ylabel('No. Samples','FontSize',10)
    subplot(1,2,2)
    histogram(log10(data.heat_production(ind,:)),'BinWidth',0.2)
    xlim([-2 2])
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    xlabel('A [\muW m^{-3}]','FontSize',10)
    textwidth = 16.99757;
    set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth/2])
    
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])
    %pdf
%    print(fig,['/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/' ...
%        'fig' int2str(3)],'-dpdf','-r0')
        return
%Figure 4: Heat production vs Age
    %sio2
%     ind = (~isnan(data.heat_production) & ~isnan(data.avg_age)...
%         & ~strcmpi(data.country,'ocean') & ~strcmpi(data.country,'unknown')...
%         & rockgroup(data,'igneous protolith') & data.sio2 >= 30  & data.sio2 <= 90);
    agebin = age_hp_box(data(ind,:),age_div,'sio2');


%Figure 5: Mafic vs Felsic bias
    %BE CAREFUL: mafic and felsic ind here are for current data(ind,:) not
    %total data table. Do not change ind and then use mafic and felsic ind.
    [mafic_ind,felsic_ind] = ero_proxy(data(ind,:),age_div,60);
    close
    
%Figure 6: Plots of mafic vs felsic HP through time
    hpage_cor(data(ind,:),age_div,felsic_ind,mafic_ind)
    close

%Figure 7: Showing correction for sio2
    data.hp_corrected = nan([length(data.sio2),1]);
    %[data.hp_corrected(ind,:),~] = hp_sio2(data(ind,:),sio2range(1):2:sio2range(2));
    [~,data.hp_corrected(ind,:)] = hp_sio2(data(ind,:),sio2range(1):2:sio2range(2));

%Different plots for crustal contamination checks
    %figure()
    %contam(data)
    
%Geographic location vs age bin interactive plot
%     ind = (~isnan(data.heat_production) & ~isnan(data.avg_age)...
%         & ~strcmpi(data.country,'ocean') & ~strcmpi(data.country,'unknown')...
%         & rockgroup(data,'igneous protolith'));
     geo_age(data(ind,:),age_div)
    

%Figure 8: Mafic vs felsic HP through time corrected for Sio2
    hpage_cor_2(data(ind,:),age_div,felsic_ind,mafic_ind);
    
    %Plot data in timescale bins instead - expanded data set, but lower
    %resolution bins
    %figure()
    %output = ts_agehp(data,'era');
    
    %Th/U and K/U plot
%     figure()
%     hp_ratio_time(data,age_div)

    %Liquidus plot
    %figure()
    %liq_age(data,age_div)
    
    %Tect age plot
    
%Figure 9:Total REE and other misc. vs HP
%REWRITE THIS COMPLETELY - indices are flawed.
    %totalREE(data(ind,:))
    
%Figure 10: Correction for alkali
    %figure()
    %hp_TA(data,age_div)
%     ind = (~isnan(data.heat_production) & ~isnan(data.avg_age)...
%         & ~strcmpi(data.country,'ocean') & ~strcmpi(data.country,'unknown')...
%         & rockgroup(data,'igneous protolith') & data.sio2 >= 45  & data.sio2 <= 80 ...
%         & data.k2o >= 0 & data.na2o >= 0);
    %new_hp_TA(data(ind,:),age_div)
    
%Figure 11: Mantle potential temperature through time
    %mantle_temp(data(ind,:),age_div)


%Figure 12: HP decay plot correction
    data.hp_present = nan([length(data.sio2),1]);
    data.hp_origin = nan([length(data.sio2),1]);
    [data.hp_present(ind,:),data.hp_origin(ind,:)] = hpinitial(data(ind,:),age_div,'sio2corrected');
    
    %pc_hp type plot - showing Australian distribution vs rest of world
    %histogram for specified age ranges - eons or eras?

    
%Figure 13: Correction for alkali - modified for percentiles instead
%     ind = (~isnan(data.heat_production) & ~isnan(data.avg_age)...
%         & ~strcmpi(data.country,'ocean') & ~strcmpi(data.country,'unknown')...
%         & rockgroup(data,'igneous protolith') & data.sio2 >= 45  & data.sio2 <= 80 ...
%         & data.k2o >= 0 & data.na2o >= 0);
    perc_hp_TA(data(ind,:),age_div)
    %ADD PERCENTILES/WHISKER ON THESE? I have the sets, look where i put
    %the median and change this to whisker info?
    %Split mafic and felsic?
    
    jaupartplot
    
    
    %Percentiles
    

return