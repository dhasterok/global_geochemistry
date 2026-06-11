function gciplot(data);

% since this gets used a lot, it is defined here.
ind0 = data.heat_production > 0 & data.p_velocity > 0;

figure;
subplot(231); hold on;
%p = plot(data.Fe_number,log10(data.heat_production),'.');
%set(p,'Color',[0.7 0.7 0.7]);
chemhpplot(data,'Fe_number',[0 1]);
%ftxt = {'ferroan','magnesian'};
%for i = 1:length(ftxt);
%    ind = (data.Fe_number > 0);
%    for j = 1:length(gfield{2})
%        ind(j,1) = ind(j) & strcmp(ftxt{i},gfield{1}{j});
%    end
%    
%    if isempty(find(ind))
%        continue;
%    end
%    p = plot(data.Fe_number(ind),log10(data.heat_production(ind)),'.');
%    set(p,'Color',[0.7 0.7 0.7]);
%    %set(p,'Color',[0.1 0.1 0.1]+0.2*(i-1));
%end

Fe_num = [0:0.05:1];
for i = 1:length(Fe_num)-1
    ind = ind0 & ...
        (Fe_num(i) < data.Fe_number & data.Fe_number <= Fe_num(i+1));

    Atemp{i} = data.heat_production(ind);
    Fetemp{i} = data.Fe_number(ind);
end

[X,Y] = whisker(Fetemp,Atemp,'Color',[0.5 0.5 0.5],'Scale','log');
xlabel('Fe Number');
ylabel('Heat Production');
hpax([-2 2],'y');
%golden;
    pbaspect([3 2.5 1]);
clear Atemp Fetemp

subplot(232); hold on;
mtxt = {'alkalic'; 'alkali-calcic'; 'calc-alkalic'; 'calcic'};

%p = plot(data.MALI,log10(data.heat_production),'.');
%set(p,'Color',[0.7 0.7 0.7]);
chemhpplot(data,'MALI',[-15 15]);
%ind = logical(zeros(size(data.MALI)));
%for i = 1:length(mtxt);
%    for j = 1:length(gfield{2})
%        ind(j,1) = strcmp(mtxt{i},gfield{2}{j});
%    end
%    
%    if isempty(find(ind))
%        continue;
%    end
%    p = plot(data.MALI(ind),log10(data.heat_production(ind)),'.');
%    set(p,'Color',[0.1 0.1 0.1]+0.2*(i-1));
%end
%legend(mtxt);

MALI = [-15:1.5:15];
for i = 1:length(MALI)-1
    ind = ind0 & (MALI(i) <= data.MALI & data.MALI < MALI(i+1));

    Atemp{i} = data.heat_production(ind);
    MALItemp{i} = data.MALI(ind);
end

[X,Y] = whisker(MALItemp,Atemp,'Color',[0.5 0.5 0.5],'Scale','log');
ylabel('Heat Production');
xlabel('MALI');
xlim([-15 15]);
hpax([-2 2],'y');
%golden;
    pbaspect([3 2.5 1]);
clear Atemp MALItemp

subplot(233);
%p = plot(data.ASI,log10(data.heat_production),'.');
%set(p,'Color',[0.7 0.7 0.7]);
chemhpplot(data,'ASI',[0.5 2]);
ASI = [0.5:0.1:2];
for i = 1:length(ASI)-1
    ind = ind0 & (ASI(i) <= data.ASI & data.ASI < ASI(i+1));

    Atemp{i} = data.heat_production(ind);
    ASItemp{i} = data.ASI(ind);
end

[X,Y] = whisker(ASItemp,Atemp,[0.5 0.5 0.5],'log');
ylabel('Heat Production');
xlabel('ASI');
xlim([0.5 2]);
hpax([-2 2],'y');
%golden;
    pbaspect([3 2.5 1]);
clear Atemp ASItemp


subplot(234);
%p = plot(data.maficity,log10(data.heat_production),'.');
%set(p,'Color',[0.7 0.7 0.7]);
chemhpplot(data,'maficity',[0 1]);
maficity = [0:0.05:0.8];
for i = 1:length(maficity)-1
    ind = ind0 & ...
        (maficity(i) <= data.maficity & data.maficity < maficity(i+1));

    Atemp{i} = data.heat_production(ind);
    Mtemp{i} = data.maficity(ind);
end

[X,Y] = whisker(Mtemp,Atemp,'Color',[0.5 0.5 0.5],'Scale','log');
ylabel('Heat Production');
xlabel('Maficity');
xlim([0 1]);
hpax([-2 2],'y');
%golden;
    pbaspect([3 2.5 1]);
clear Atemp Mtemp


subplot(235);
%p = plot(data.CIA_molar,log10(data.heat_production),'.');
%set(p,'Color',[0.7 0.7 0.7]);
chemhpplot(data,'CIA',[30 70]);
CIA = [30:2:70];
for i = 1:length(CIA)-1
    ind = ind0 & (CIA(i) <= data.CIA & data.CIA < CIA(i+1));

    Atemp{i} = data.heat_production(ind);
    CIAtemp{i} = data.CIA(ind);
end

[X,Y] = whisker(CIAtemp,Atemp,'Color',[0.5 0.5 0.5],'Scale','log');
ylabel('Heat Production');
xlabel('CIA (molar)');
xlim([30 70]);
hpax([-2 2],'y');
%golden;
pbaspect([3 2.5 1]);
clear Atemp CIAtemp

subplot(236);
%p = plot(data.WIP,log10(data.heat_production),'.');
%set(p,'Color',[0.7 0.7 0.7]);
chemhpplot(data,'WIP',[40 140]);
WIP = [40:5:140];
for i = 1:length(WIP)-1
    ind = ind0 & (WIP(i) <= data.WIP & data.WIP < WIP(i+1));

    Atemp{i} = data.heat_production(ind);
    WIPtemp{i} = data.WIP(ind);
end

[X,Y] = whisker(WIPtemp,Atemp,'Color',[0.5 0.5 0.5],'Scale','log');
ylabel('Heat Production');
xlabel('WIP (molar)');
xlim([40 140]);
hpax([-2 2],'y');
%golden;
    pbaspect([3 2.5 1]);
clear Atemp WIPtemp

return