function oxplot(data,varargin)
% oxplot(data)

% igneous axes (default)
ax = [30 100; 0 5; 0 30; 0 20; 0 20; 0 20; 0 10; 0 10; 0 2.5];
if nargin == 2
    if varargin{1} == 1
        % sedimentary (axes)
        ax = [30 100; 0 3; 0 30; 0 20; 0 20; 0 80; 0 10; 0 10; 0 2.5];
    end
end

%heat production plots
oxlist = {'SiO2','TiO2','Al2O3','FeO_TOT','MgO','CaO','Na2O','K2O','P2O5'};



fig1 = figure;
fig2 = figure;
fig3 = figure;
fig4 = figure;
for i = 1:length(oxlist)
    % heat production plots
    figure(fig1);
    subplot(3,3,i);
    chemhpplot(data,oxlist{i},ax(i,:));
    
    % velocity plots
    figure(fig2);
    subplot(3,3,i);
    chemplot(data,oxlist{i},ax(i,:),'wt.%','p_velocity',[5.8 8.4],'km s^{-1}');
    
    % density plots
    figure(fig3);
    subplot(3,3,i);
    chemplot(data,oxlist{i},ax(i,:),'wt.%','density',[2400 3600],'kg m^{-3}');
    
    % harker diagrams colored by heat production
    if i ~= 1
        figure(fig4);
        subplot(3,3,i);
        harker(data,oxlist{i},ax(1,:),ax(i,:));
        line([100,0],[100-ax(i,2),ax(i,2)])
    end
    ax(i,2)
end

return


function harker(data,el,xl,yl);

sio2 = data.sio2;
d = data{:,lower(el)};
ind = (sio2 > 0 & d > 0 & data.heat_production > 0);
%ec = linspace(xl(1),xl(2),40);
%ea = [-3:0.1:2];

%n = log10(hist2d(d(ind),log10(data.heat_production(ind)),ec,ea));
scatter(data.sio2(ind),d(ind),1,log10(data.heat_production(ind)),'.');
set(gca,'CLim',[-2 2]);
%[c,h] = contourf(ec(1:end-1) + diff(ec)/2,ea(1:end-1) + diff(ea)/2,n);
%colormap(flipud(gray));
xlabel(['SiO2 (wt. %)']);
ylabel([el,' (wt. %)']);
xlim(xl);
ylim(yl);


return
