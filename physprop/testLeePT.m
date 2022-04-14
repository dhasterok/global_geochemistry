close all;
clear all;

data = readtable('../data/LeePT/Lee_etal2009.xlsx');
prim = basaltPT(data);

figure
subplot(221);
scatter(data.T,prim.temperature,'filled');
hold on;
plot([1300 2000],[1300 2000]);
axis([1300 2000 1300 2000]);
xlabel('T_{LEE} (^\circC)');
ylabel('T_{MATLAB} (^\circC)');
axis square;
set(gca,'Box','on');

subplot(222);
scatter(data.P,prim.pressure,'filled');
hold on;
plot([0 18],[0 18]);
axis([0 18 0 18]);
xlabel('P_{LEE} (GPa)');
ylabel('P_{MATLAB} (GPa)');
axis square
set(gca,'Box','on');

subplot(223);
scatter(data.T,data.P,'filled');
hold on;
scatter(prim.temperature,prim.pressure);
legend('Lee et al.','Matlab','Location','southwest');
ylim([0 20]);
axis ij;
axis square;
ylabel('P (GPa)');
xlabel('T (^\circC)');
set(gca,'Box','on');