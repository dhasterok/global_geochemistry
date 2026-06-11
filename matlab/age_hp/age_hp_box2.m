function agebin = age_hp_box2(data,age_div,name,varargin)
%--------------------------------------------------------------------------
%
%                 Age vs. HP (w/ chemical distributions)
%
%--------------------------------------------------------------------------

%Files pre-loaded: data, plutonic, volcanic, sedimentary
%Outputs: Stats includes median, ranges, number of data points
%age_div = [0:200:4000];

%Plot only mafic volcanics for Nb/La plot?
%agebin = age_hp_box2(data(data.sio2<60 & strcmpi(data.rock_origin,{'volcanic'}),:),[0:200:4000],'All data','sio2')




%data in table format

flag = 0;
if nargin == 4
    flag = 1;
    field = varargin{1};
end

%Too large?
max_dage = 200;

%Needs more conditions on these - see age_hp_error, including neg age, and
%single column values in max or min age
data.avg_age = data.age;
ind = isnan(data.avg_age);

dage = (data.age_max - data.age_min);
dage(dage > max_dage) = NaN;
  
data.avg_age(ind) = data.age_min(ind) + dage(ind)/2;

figure;
for i = 1:length(age_div)-1
    ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) ...
        & ~isnan(data.heat_production) & data.sio2 > 25 & ~strcmp(data.country,'unknown'));
    n(i) = length(ind);
    agebin.ind{i} = ind;
    avg_age{i} = data.avg_age(ind);
    heat_production{i} = data.heat_production(ind);
    
    if flag
        field_val{i} = data{ind,lower(field)};
    end
end
if flag
    subplot(311);
else
    subplot(211);
end
[agebin.Qage,agebin.Qhp] = whisker(avg_age,heat_production,'Color',[0.5 0.5 0.5],'Scale','log');
xlim([age_div(1) age_div(end)]);
ylabel('Heat Production [\muW m^{-3}]');
title(name);
set(gca,'Box','on');

if flag
    subplot(312);
else
    subplot(212);
end

n(n == 0) = NaN;
logn = log10([n';NaN]);
h = bar(age_div,logn,'histc');
hpax([0 ceil(max(logn))]);
set(h,'FaceColor',[0.5 0.5 0.5]);
xlim([age_div(1) age_div(end)]);
ylabel('Number of Samples');
age_div'
n'

if ~flag
    xlabel('Age [Ma]');
    return
end


subplot(313); hold on;
for i = 1:length(field_val)
    
    if isempty(field_val{i})
        continue;
    end
    
    N = 2/(age_div(i+1) - age_div(i));
    mirrhist(field_val{i},N,agebin.Qage(i,3),'y');
end
xlim([age_div(1) age_div(end)]);
ylabel(field);
xlabel('Age [Ma]');
set(gca,'Box','on');

figure()
[mafic_ind,felsic_ind] = ero_proxy(data,age_div);

figure()
subplot(2,1,1)
histage(data,mafic_ind,age_div,'mafic')
subplot(2,1,2)
histage(data,felsic_ind,age_div,'felsic')

figure()
hpage_cor(data,age_div,agebin.ind,mafic_ind,felsic_ind)

%figure()
% alt_index(data,agebin.ind,age_div)

%figure()
%tect_age = tectage(data);

figure()
contam(data);

return


function mirrhist(val,N,shift,coord)

if isempty(val)
    warning('Value vector was empty.');
    return;
end

edges = linspace(min(val),max(val),round(sqrt(length(val))));

x = [edges; edges];
x = x(:);

x = [x; flipud(x)];

n = histc(val,edges)';

if N ~=0
    n = n/(N*max(n));
end


y = [0 n(1:end-1); n(1:end-1) 0];
y = y(:);

y = [y; -flipud(y)] + shift;

if strcmp(coord,'y');
    fill(y,x,[0.5 0.5 0.5]);
else
    fill(x,y,[0.5 0.5 0.5]);
end

hline = refline([0 65]);
hline.Color = 'r';

hline = refline([0 55]);
hline.Color = 'b';

return

function [mafic_ind,felsic_ind] = ero_proxy(data,age_div)
    for i = 1:length(age_div)-1
        mafic_ind{i} = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) ...
        & ~isnan(data.heat_production) & data.sio2 > 25 & ~strcmp(data.country,'unknown')...
        & data.sio2<=60 & data.sio2>0);
    
        felsic_ind{i} = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) ...
        & ~isnan(data.heat_production) & data.sio2 > 25 & ~strcmp(data.country,'unknown')...
        & data.sio2>60 & data.sio2<=100);

        felsic_num(i) = length(felsic_ind{i});
        mafic_num(i) = length(mafic_ind{i});
    end
    
    h = bar(age_div(1:end-1),(felsic_num-mafic_num)./(felsic_num+mafic_num),'histc');
    set(h,'FaceColor',[0.5 0.5 0.5])
    ylabel('Composition bias [+ve Felsic, -ve Mafic]');
    xlabel('Age [Ma]');
    title('Felsic vs. Mafic bias');
    set(gca,'Box','on')
return

function histage(data,rock_ind,age_div,name)
    for i = 1:length(age_div)-1
        n(i) = length(rock_ind{i});
    end

    n(n == 0) = NaN;
    logn = log10([n';NaN]);
    h = bar(age_div,logn,'histc');
    hpax([0 ceil(max(logn))]);
    set(h,'FaceColor',[0.5 0.5 0.5]);
    xlim([age_div(1) age_div(end)]);
    ylabel('Number of Samples');
    xlabel('Age [Ma]')
    title(name)
return

function hpage_cor(data,age_div,data_ind,mafic_ind,felsic_ind)
    for i = 1:length(age_div)-1
        avg_age{i} = data.avg_age(data_ind{i});
        heat_production{i} = data.heat_production(data_ind{i});
        
        avg_age_mafic{i} = data.avg_age(mafic_ind{i});
        heat_production_mafic{i} = data.heat_production(mafic_ind{i});
        
        avg_age_felsic{i} = data.avg_age(felsic_ind{i});
        heat_production_felsic{i} = data.heat_production(felsic_ind{i});
    end
    subplot(3,1,1)
    whisker(avg_age,heat_production,'Color',[0.5 0.5 0.5]);
    xlim([age_div(1) age_div(end)]);
    set(gca,'Box','on');
    title('HP vs. Age: All data')
    ylim([0 10])
    
    subplot(3,1,2)
    whisker(avg_age_felsic,heat_production_felsic,'Color',[0.5 0.5 0.5]);
    xlim([age_div(1) age_div(end)]);
    set(gca,'Box','on');
    title('HP vs. Age: Felsic')
    ylim([0 10])
    
    subplot(3,1,3)
    whisker(avg_age_mafic,heat_production_mafic,'Color',[0.5 0.5 0.5]);
    xlim([age_div(1) age_div(end)]);
    set(gca,'Box','on');
    xlabel('Age [Ma]')
    title('HP vs. Age: Mafic')
    ylim([0 10])
    
    figure()
    subplot(121);
    title('SiO2 vs HP slope through time')
    hold on
    for i = 1:((length(age_div)-1)/2)
        p(i,:) = polyfit(vertcat(data.sio2(data_ind{i}),data.sio2(data_ind{i+1})), ...
                log10(vertcat(data.heat_production(data_ind{i}),data.heat_production(data_ind{i+1}))),1);
        plot(age_div(2*i),p(i,1),'o')
    end
    
    hold off
    ylabel('Slope')
    xlabel('Age [Ma]')
    axis square;
    
    subplot(122); hold on;
    for i = 1:((length(age_div)-1)/2)
        plot([35 85],polyval(p(i,:),[35 85]),'-');
    end
    hpax([-2 2]);
    axis square;
    
    
    figure()
    title('SiO2 vs HP all data')
    data_ind_all = [];
    for i = 1:length(age_div)-1
        data_ind_all = vertcat(data_ind_all,data_ind{i});
    end
    p = polyfit(data.sio2(data_ind_all),log10(data.heat_production(data_ind_all)),1);
    y = polyval(p,data.sio2(data_ind_all));
    hold on
    plot(data.sio2(data_ind_all),log10(data.heat_production(data_ind_all)),'.')
    plot(data.sio2(data_ind_all),y,'-r')
    hold off
    hpax([-2 2]);
    xlabel('SiO2')
    
return

function alt_index(data,data_ind,age_div)
    
% granite classification system Frost et al. 2001
disp('Determining granite class based on Frost et al. (2001)...');
data = graniteclass(data);

gfield = getgfield(data);

%save OzChem_processed.mat data

gclass = {'magnesian>alkalic>metaluminous';
    'magnesian>alkali-calcic>metaluminous';
    'magnesian>calc-alkalic>metaluminous';
    'magnesian>calcic>metaluminous';
    'magnesian>alkalic>peraluminous';
    'magnesian>alkali-calcic>peraluminous';
    'magnesian>calc-alkalic>peraluminous';
    'magnesian>calcic>peraluminous';
    'ferroan>alkalic>metaluminous';
    'ferroan>alkali-calcic>metaluminous';
    'ferroan>calc-alkalic>metaluminous';
    'ferroan>calcic>metaluminous';
    'ferroan>alkalic>peraluminous';
    'ferroan>alkali-calcic>peraluminous';
    'ferroan>calc-alkalic>peraluminous';
    'ferroan>calcic>peraluminous'};

gbet = {'S','I','A'};

gcolor = [110  65  60;
    139  73  54;
    184  77  41;
    213  77  29;
    234 118  36;
    240 156  38;
    247 190  41;
    254 226  44;
    163  25  91;
    149  27 129;
    102  36 131;
     41  35  92;
     45  46 131;
     29 113 184;
     54 169 225;
     47 172 102]/255;

for j = 1:length(gclass)
    gdata(j).A = [];
    gdata(j).Vp = [];
end

for i = 1:length(data.sio2)
    if isnan(data.heat_production(i)) | isnan(data.p_velocity(i))
        continue;
    end

    c = [];
    for j = 1:length(gclass)
        if strcmp(data.gclass(i),gclass{j})
            gdata(j).A = [gdata(j).A; data.heat_production(i)];
            gdata(j).Vp = [gdata(j).Vp; data.p_velocity(i)];
            break;
        end
    end
end

eV = [5.8:0.2:8.2];
eA = [-2:0.1:2];

Vp = [6 8.2];

    chemhpplot(data,'CIA_molar',[0 2]);
    CIA = [0:0.1:2];
    for i = 1:length(CIA)-1
        ind = (CIA(i) <= data.CIA_molar ...
            & data.CIA_molar < CIA(i+1));

        Atemp{i} = data.heat_production(ind);
        CIAtemp{i} = data.CIA_molar(ind);
    end
    subplot(2,1,1)
    [X,Y] = whisker(CIAtemp,Atemp,'Color',[0.5 0.5 0.5],'Scale','log');
    ylabel('Heat Production');
    xlabel('CIA (molar)');
    xlim([0 2]);
    hpax([-2 2],'y');
    %golden;
    pbaspect([3 2.5 1]);
    
    
    for i = 1:length(age_div)-1
        age_temp{i} = data.avg_age(data_ind{i});
        CIA_temp{i} = data.CIA_molar(data_ind{i});
    end
    
    subplot(2,1,2)
    [X,Y] = whisker(age_temp,CIA_temp,'Color',[0.5 0.5 0.5]);
    ylabel('CIA (molar)');
    xlabel('Age [Ma]');
    xlim([0 4000]);

return

function gfield = getgfield(data);

for i = 1:length(data.gclass)
    tmp = textscan(data.gclass{i},'%s %s %s','Delimiter','>');
    gfield{1}{i} = tmp{1};
    gfield{2,1}{i} = tmp{2};
    gfield{3}{i} = tmp{3};
end

return

function tect_age = tectage(data)
%     % check locations
%     ind = find(data.latitude < -90 | data.latitude > 90 ...
%         | data.longitude < -180 | data.longitude > 180);
%     if ~isempty(ind)
%         warning('Some latitude(s) or longitude(s) are out of bounds. Ignoring data...');
%     end
%     data.latitude(ind) = NaN;
%     data.longitude(ind) = NaN;
% 
%     nonan = (~isnan(data.latitude) & ~isnan(data.longitude) & ~isnan(data.heat_production));
%     
%     H.lat = data.latitude(nonan);
%     H.lon = data.longitude(nonan);
%     H.A = data.heat_production(nonan);
% 
%     % Load tectonic age
%     fprintf('\nLoading tectonic age...\n');
%     tage = worldgrid('tage');
%     % Find index of tectonic age data corresponding to heat flow locations
%     latind = round((90 - H.lat)/dl) + 1;
%     lonind = round((180 + H.lon)/dl) + 1;
%     
%     [nrow min(latind) max(latind) ncol min(lonind) max(lonind)]
%     ind = sub2ind([nrow,ncol],latind,lonind);
%     tect_age.tage = tage(ind);
%     tage_key = readkey('tage.key');
%     
%     
%     for i = 1:length(tage_key.val)
%         tect_age.tmin(i) = tage_key.agemin(i);
%         tect_age.tmax(i) = tage_key.agemax(i);
% 
%         ind = find(tage == tage_key.val(i));
%         tect_age.hp{i} = H.A(ind);
%     end
tect_age = 1;
return

function contam(data)
    ind = find(~isnan(data.nb_ppm) & data.nb_ppm>=0 & ...
        ~isnan(data.la_ppm) & data.la_ppm>=0 & ... 
        ~isnan(data.heat_production) & ~isnan(data.u_ppm) & ~isnan(data.th_ppm)...
        & ~isnan(data.k_ppm));
    %plot(sort(data.nb_ppm(ind)./data.la_ppm(ind),'ascend'),'.')
    h = histogram(data.nb_ppm(ind)./data.la_ppm(ind),[0:0.1:3]);
    xlim([0 3])
    title('Nb/La ratio')
    hold on
    %hline = refline([0 1]);
    %hline.Color = 'r';
    line([1,1],ylim,'LineWidth',2,'Color','r')
    hold off
    set(h,'FaceColor',[0.5 0.5 0.5]);
    set(gca,'Box','on');
    %ylim([0 4])
    %xlim([0,length(data.heat_production(ind))])
    
    figure()
    %plot(sort(data.th_ppm(ind)./data.u_ppm(ind),'ascend'),'.')
    h = histogram(data.th_ppm(ind)./data.u_ppm(ind),[0:0.5:20]);
    xlim([0 20])
    title('Th/U ratio')
    hold on
    %hline = refline([0 4]);
    %hline.Color = 'r';
    line([4,4],ylim,'LineWidth',2,'Color','r')
    hold off
    set(h,'FaceColor',[0.5 0.5 0.5]);
    set(gca,'Box','on');
    %ylim([0 16])
    %xlim([0,length(data.heat_production(ind))])
    
    figure()
    subplot(1,2,1)
    %plot(sort(data.k_ppm(ind)./data.u_ppm(ind),'ascend'),'.')
    h = histogram(data.k_ppm(ind)./data.u_ppm(ind),[0:0.1*(10^4):5*(10^4)]);
    xlim([0 5*(10^4)])
    title('Kppm/U ratio')
    %ylim([0,0.4*(10^5)])
    %xlim([0,length(data.heat_production(ind))])
    set(h,'FaceColor',[0.5 0.5 0.5]);
    set(gca,'Box','on');
    
    subplot(1,2,2)
    %k2o in wt%
    %ans * 10000 is in ppm
    %94.1906 is molecular weight of k2o
    %k2 component is 78.1966
    %k2 in ppm = 78.1966./94.1906
    
    h = histogram((data.k2o(ind).*10000.*(78.1966./94.1906))./data.u_ppm(ind),[0:0.1*(10^4):5*(10^4)]);
    %xlim([0 5*(10^4)])
    title('K/U ratio')
    %ylim([0,0.4*(10^5)])
    %xlim([0,length(data.heat_production(ind))])
    set(h,'FaceColor',[0.5 0.5 0.5]);
    set(gca,'Box','on');
    xlim([0 5*(10^4)])
    
    figure()
    subplot(121);
    plot(log10(data.la_ppm(ind)./data.th_ppm(ind)),log10(data.th_ppm(ind)./data.u_ppm(ind)),'.')
    hold on;
    plot(log10([1 10 10 1 1]),log10([2 2 10 10 2]),'-');
    title('Th/La vs Th/U')
    hpax([-1 3],'x')
    xlabel('Th/La')
    hpax([-1 3],'y')
    ylabel('Th/U')
    axis square;
    
    subplot(122);
    plot(log10(data.u_ppm(ind)),log10(data.th_ppm(ind)),'.');
    hold on;
    plot([-2 2],[-2 2]+log10(2),'-');
    plot([-2 2],[-2 2]+log10(10),'-');
    hpax([-1 2],'x');
    xlabel('U (ppm)');
    hpax([-1 2],'y');
    ylabel('Th (ppm)');
    axis square;
return