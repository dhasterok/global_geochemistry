function [hppresent, hporigin] = hpinitial(data,age_div,fig,varargin)
%Calculates and plots HP at formation i.e. correction for current day
%values to age it was formed.

%change this to use U/Th/K ratios instead

   
Exist_Column = strcmp('avg_age',data.Properties.VariableNames);
val = Exist_Column(Exist_Column==1);
if isempty(val)
    data.avg_age = data.age;
end

%Convert k_ppm to k2o if not listed
ind = data.k_ppm >= 0 & isnan(data.k2o);
data.k2o(ind,:) = (10^4)*data.k_ppm(ind,:);


%Only plot for k,u, and th positive, avg_age value
ind = find(data.k2o >= 0 & data.u_ppm >= 0 & data.th_ppm >= 0 & data.avg_age >= 0 & data.density_model > 0);


data.hp_present = nan([length(data.sio2),1]);
data.hp_origin = nan([length(data.sio2),1]);

if nargin > 2
    %sio2 corrected have no U/Th corrected... use median curve
    data.hp_present = data.hp_corrected;
    %Use median values and exponential curve for them as correction
    U = nanmedian(data.u_ppm);
    Th = nanmedian(data.th_ppm);
    Rb = nanmedian(data.rb_ppm);
    Sm = nanmedian(data.sm_ppm);
    Density = nanmedian(data.density_model);
    K2O = nanmedian(data.k2o);
    [Density K2O U Th Rb]
    [present,~] = radtime(Density,K2O,0,0,Th,U,'K2O','Formula','r88');
    h = waitbar(0,'Please wait...');
    for i = 1:length(ind)
        [past,~] = radtime(Density,K2O,0,0,Th,U,'K2O','Age',data.avg_age(ind(i),:),'Formula','r88');
        ratioincrease = past/present;
        data.hp_origin(ind(i),:) = data.hp_present(ind(i),:) * ratioincrease;
        waitbar(i / length(ind))
    end
    close(h)
    
    %fig = figure()
    %subplot(2,2,3:4)
    for i = 1:length(age_div)-1
        ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1));
        n(i) = length(ind);
       agebin.ind{i} = ind;
       avg_age{i} = data.avg_age(ind);
       hp_origin{i} = data.hp_origin(ind);
       hp_present{i} = data.hp_present(ind);
    end

    [agebin.Qage,agebin.Qhp] = whisker(avg_age,hp_origin,'Color',[0.5 0.5 0.5],'Scale','log');
    [agebin1.Qage,agebin1.Qhp] = whisker(avg_age,hp_present,'Color',[0.5 0.5 0.5],'Scale','log');
    plot(0,0)
    sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,1],'Original')
    hold on
    sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,[1,0,0],'Present')
    hold off
    xlim([age_div(1) age_div(end)]);
    
    set(gca,'Box','on');
    hpax([-1 2]);

    h = findobj(gca,'Type','patch');
    [lh,lhicons] = legend(h([2 1]),'Original','Present');
    patchinlegend = findobj(lhicons,'type','patch');
    set(patchinlegend,'facea',0.3)
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    ylabel('A [\muW m^{-3}]','FontSize',10);
    xlabel('Age [Ma]','FontSize',10)
    lh.FontSize = 9;
    % Get current position (which I used to keep overall size the same)
    currentLegendPosition = lh.Position;
    % Define new position
    newLegendPosition = [0.709 0.865 currentLegendPosition([3 4])];
    % Set new position
    lh.Position = newLegendPosition;
% 
%     subplot(2,2,2)
%     sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,1],'Original')
%     xlim([age_div(1) age_div(end)]);
%     
%     set(gca,'Box','on');
%     hpax([-1 2]);
%     set(gca,'Units','normalized',...
%         'FontUnits','points',...
%         'FontWeight','normal',...
%         'FontSize',8);
%     ylabel('A [\muW m^{-3}]','FontSize',10);
%     xlabel('Age [Ma]','FontSize',10)
    

%     subplot(2,2,1)
%     sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,[1,0,0],'Present')
%     xlim([age_div(1) age_div(end)]);
%     
%     set(gca,'Box','on');
%     hpax([-1 2]);
%     set(gca,'Units','normalized',...
%         'FontUnits','points',...
%         'FontWeight','normal',...
%         'FontSize',8);
%     ylabel('A [\muW m^{-3}]','FontSize',10);
%     xlabel('Age [Ma]','FontSize',10)
%     
% 
    hppresent = data.hp_present;
    hporigin = data.hp_origin;
    
    textwidth = 16.99757;
%     set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth])
%     
%     pos = get(fig,'Position');
%     set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
%     'PaperSize',[pos(3), pos(4)])
    
    
    return
end
    
    

h = waitbar(0,'Please wait...');
for i = 1:length(ind)
    zero_vector = zeros(size(data.sm_ppm(ind(i),:)));
    [data.hp_origin(ind(i),:),~] = radtime(data.density_model(ind(i),:),data.k2o(ind(i),:),zero_vector...
    ,zero_vector,data.th_ppm(ind(i),:),data.u_ppm(ind(i),:),'K2O','Age',data.avg_age(ind(i),:),'Formula','r88');
    waitbar(i / length(ind))
end
close(h)





zero_vector = zeros(size(data.sm_ppm(ind,:)));
[data.hp_present(ind,:),~] = radtime(data.density_model(ind,:),data.k2o(ind,:),zero_vector...
    ,zero_vector,data.th_ppm(ind,:),data.u_ppm(ind,:),'K2O','Formula','r88');


h = waitbar(0,'Please wait...');
for i = 1:length(ind)
    zero_vector = zeros(size(data.sm_ppm(ind(i),:)));
    [data.hp_origin(ind(i),:),~] = radtime(data.density_model(ind(i),:),data.k2o(ind(i),:),zero_vector...
    ,zero_vector,data.th_ppm(ind(i),:),data.u_ppm(ind(i),:),'K2O','Age',data.avg_age(ind(i),:),'Formula','r88');
    waitbar(i / length(ind))
end
close(h)


fig = figure()
subplot(2,2,3:4)
for i = 1:length(age_div)-1
    ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1));
    n(i) = length(ind);
    agebin.ind{i} = ind;
    avg_age{i} = data.avg_age(ind);
    hp_origin{i} = data.hp_origin(ind);
    hp_present{i} = data.hp_present(ind);
end

[agebin.Qage,agebin.Qhp] = whisker(avg_age,hp_origin,'Color',[0.5 0.5 0.5],'Scale','log');
[agebin1.Qage,agebin1.Qhp] = whisker(avg_age,hp_present,'Color',[0.5 0.5 0.5],'Scale','log');
plot(0,0)
sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,1],'Original')
hold on
sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,[1,0,0],'Present')
hold off
xlim([age_div(1) age_div(end)]);
set(gca,'Box','on');
hpax([-1 2]);
h = findobj(gca,'Type','patch');
    [lh,lhicons] = legend(h([2 1]),'Original','Present');
    patchinlegend = findobj(lhicons,'type','patch');
    set(patchinlegend,'facea',0.3)
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    ylabel('A [\muW m^{-3}]','FontSize',10);
    xlabel('Age [Ma]','FontSize',10)
    lh.FontSize = 10;

subplot(2,2,2)
sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,1],'Original')
xlim([age_div(1) age_div(end)]);
set(gca,'Box','on');
hpax([-1 2]);
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    ylabel('A [\muW m^{-3}]','FontSize',10);
    xlabel('Age [Ma]','FontSize',10)

subplot(2,2,1)
sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,[1,0,0],'Present')
xlim([age_div(1) age_div(end)]);

set(gca,'Box','on');
hpax([-1 2]);
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    ylabel('A [\muW m^{-3}]','FontSize',10);
    xlabel('Age [Ma]','FontSize',10)


hppresent = data.hp_present;
hporigin = data.hp_origin;

textwidth = 16.99757;
    set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth])
    
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])
    %pdf
    

% %Present day heat production:
% data.hp_present = nan([length(data.sio2),1]);
% data.hp_present(ind,:) = data.density_model(ind,:).*(...
%     9.52.*(data.u_ppm(ind,:))...
%     +2.56.*(data.th_ppm(ind,:))...
%     +3.48.*(data.k2o(ind,:)*2*molecularwt('K')/molecularwt('K2O'))...
%     ).*10^(-5);
% 
% %Original heat production:
% data.hp_origin = nan([length(data.sio2),1]);
% Atot = zeros(0,0);
% 
% h = waitbar(0,'Please wait...');
% for i = 1:length(ind)
%     [data.hp_origin(ind(i),:),~] = radtime(data.density_model(ind(i),:),data.k2o(ind(i),:),...
%     data.u_ppm(ind(i),:),data.th_ppm(ind(i),:),data.rb_ppm(ind(i),:),data.avg_age(ind(i),:));
%     waitbar(i / length(ind))
% end
% close(h)
% hp = data.hp_origin;
% 
% figure()
% subplot(2,2,1:2)
% for i = 1:length(age_div)-1
%     ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1));
%     n(i) = length(ind);
%     agebin.ind{i} = ind;
%     avg_age{i} = data.avg_age(ind);
%     hp_origin{i} = data.hp_origin(ind);
%     hp_present{i} = data.hp_present(ind);
% end
% 
% [agebin.Qage,agebin.Qhp] = whisker(avg_age,hp_origin,'Color',[0.5 0.5 0.5],'Scale','log');
% [agebin1.Qage,agebin1.Qhp] = whisker(avg_age,hp_present,'Color',[0.5 0.5 0.5],'Scale','log');
% sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,'b','Original')
% hold on
% sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,'r','Present')
% hold off
% xlim([age_div(1) age_div(end)]);
% ylabel('Heat Production [\muW m^{-3}]');
% xlabel('Age [Ma]')
% set(gca,'Box','on');
% 
% subplot(2,2,3)
% sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,'b','Original')
% xlim([age_div(1) age_div(end)]);
% ylabel('Heat Production [\muW m^{-3}]');
% xlabel('Age [Ma]')
% set(gca,'Box','on');
% 
% subplot(2,2,4)
% sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,'r','Present')
% xlim([age_div(1) age_div(end)]);
% ylabel('Heat Production [\muW m^{-3}]');
% xlabel('Age [Ma]')
% set(gca,'Box','on');

return