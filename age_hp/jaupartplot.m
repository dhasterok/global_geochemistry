

%Median values
Density = 2800;
K2O = 2.1;
U = 1.6;
Th = 5.9;
Rb = 54;
Sm = 0;

%Crust line 1
Density_crust = 2600;
K2O_crust = K2O*2;
U_crust = U*2;
Th_crust = Th*2;

%Mantle line 1
Density_mantle = 3000;
K2O_mantle = K2O/1.25;
U_mantle = U/1.25;
Th_mantle = Th/1.25;


%Mantle above
% Density_mantle2 = 3000;
% K2O_mantle2 = K2O/2;
% U_mantle2 = U/2*(1.25/2);
% Th_mantle2 = Th/2*(4.6/2);



%Mantle below 1
factorm3 = 0.68;
Density_mantle3 = 3000;
K2O_mantle3 = K2O/0.6*factorm3;
U_mantle3 = U/1.25*factorm3;
Th_mantle3 = Th/1.25*factorm3;

%Mantle below 2
factorm2 = 0.29;
Density_mantle2 = 3000;
K2O_mantle2 = K2O/0.19*factorm2;
U_mantle2 = U/1.25*factorm2;
Th_mantle2 = Th/1.25*factorm2;






%Crust below
factor3 = 1.029;
Density_crust3 = 2600;
K2O_crust3 = K2O*2*factor3;
U_crust3 = U*3.5*(2.7/4)*factor3;
Th_crust3 = Th*1*(2.70/4)*factor3;

%Crust above
factor2 = 8.3;
Density_crust2 = 2600;
K2O_crust2 = K2O/2*factor2/100;
U_crust2 = U/2*factor2/1.3;
Th_crust2 = Th/2*factor2;








age = [0:100:4000];
hpval_crust = [];
hpval_mantle = [];
hpval_crust2 = [];
hpval_mantle2 = [];
hpval_crust3 = [];
hpval_mantle3 = [];

for i = 1:length(age)
    %MANTLE
    [hpval_mantle(i),~] = radtime(Density_mantle,K2O_mantle,0,0,Th_mantle,U_mantle,'K2O','Age',age(i),'Formula','r88');
    [hpval_mantle2(i),~] = radtime(Density_mantle2,K2O_mantle2,0,0,Th_mantle2,U_mantle2,'K2O','Age',age(i),'Formula','r88');
    [hpval_mantle3(i),~] = radtime(Density_mantle3,K2O_mantle3,0,0,Th_mantle3,U_mantle3,'K2O','Age',age(i),'Formula','r88');

    %CRUST
    [hpval_crust(i),~] = radtime(Density_crust,K2O_crust,0,0,Th_crust,U_crust,'K2O','Age',age(i),'Formula','r88');
    [hpval_crust3(i),~] = radtime(Density_crust3,K2O_crust3,0,0,Th_crust3,U_crust3,'K2O','Age',age(i));
    [hpval_crust2(i),~] = radtime(Density_crust2,K2O_crust2,0,0,Th_crust2,U_crust2,'K2O','Age',age(i));
    
end



%PLOT THE FIGURES
fig = figure()
subplot(2,1,1)
hold on


%'Color' and 'LineStyle'

%CRUST
plot(age,hpval_crust,'Color',[0 0 1],'LineStyle','-')
plot(age,hpval_crust3,'Color',[0 0.25 0.75],'LineStyle','--')
plot(age,hpval_crust2,'Color',[0 0.5 0.5],'LineStyle','--')




%MANTLE
plot(age,hpval_mantle,'Color',[1 0 0],'LineStyle','-')
plot(age,hpval_mantle2,'Color',[210/255 105/255 30/255],'LineStyle','--')
plot(age,hpval_mantle3,'Color',[205/255 0 0],'LineStyle','--')







%PLOT THE VERTICAL LINES
ind = (age == 1000);
plot([age(ind),age(ind)],[hpval_crust2(ind), hpval_mantle2(ind)],'-ok')

ind = (age == 900);
plot([age(ind),age(ind)],[hpval_crust(ind), hpval_mantle(ind)],'LineStyle','-','Marker','o','Color',[104/255 34/255 139/255])


ind = (age == 3000);
plot([age(ind),age(ind)],[hpval_crust2(ind), hpval_mantle2(ind)],'-ok')

ind = (age == 2900);
plot([age(ind),age(ind)],[hpval_crust(ind), hpval_mantle(ind)],'LineStyle','-','Marker','o','Color',[104/255 34/255 139/255])



hold off







lh = legend('Product','Product w/ temp reduction','Product w/ depletion',...
    'Source','Source w/ temp reduction','Source w/ depletion')
set(lh,'Location','northwest')

% 
% textwidth = 16.99757;
%     set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth/1.6])
%     set(gca,'Units','normalized',...
%         'FontUnits','points',...
%         'FontWeight','normal',...
%         'FontSize',8);
% lh.FontSize = 9;

xlabel('Time [Ma]','FontSize',10)
ylabel('Formation A [\muW m^{-3}]','FontSize',10)
set(gca,'box','on')
golden

% pos = get(fig,'Position');
%     set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
%     'PaperSize',[pos(3), pos(4)])
%export_fig '/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/jaupart.pdf' -q101

subplot(2,1,2)
x = [0:200:4000];
y = 3*exp(-x/3000);
plot(x,y,'ok')
hold on
%plot(x,y,'--k')
hold off
ylim([0 8])
xlabel('Age of sample [Ma]','FontSize',10)
ylabel('Present day A [\muW m^{-3}]','FontSize',10)
set(gca,'box','on')
golden

