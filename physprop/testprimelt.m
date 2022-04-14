close all;
clear all;

data = readtable('../data/PRIMELT3/Herzberg&Asimow2015.xlsx');

P = 1;
fe_ratio = [0.9; 0.9; 0.9; 0.9; 0.9; 0.9; 0.9; 0.88; NaN];
feti = nan(height(data),1);
feti(end) = 0.89547;

%profile on;
[afm,bm] = primelt3v(data,P,'fe_ratio',fe_ratio,'feti',feti);
%profile off;

figure;
plot([1100 1800],[1100 1800],'-');
hold on;
h(1) = scatter(data.Tp_bm,bm.potential_temperature,'filled');
h(2) = scatter(data.Tp_afm,afm.potential_temperature,'filled');
xlabel('PRIMELT3 T (ºC)');
ylabel('Matlab T (ºC)');
legend(h, 'batch','AFM','Location','southeast');
axis square;
set(gca,'Box','on');