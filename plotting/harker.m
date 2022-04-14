function harker(data,varargin);

if nargin == 2
    C = varargin{1};
else
    C = [0.6 0.6 0.6];
end

switch length(C)
    case 3
        hsscatter(data,C);
    case 1
        if C == 1
            hshist(data);
        else
            hsscatter(data,C);
        end
end

return


function hsscatter(data)

sax = [0 100];

subplot(331); hold on;
plot(data.sio2,data.tio2, ...
    'ko','MarkerSize',4,'MarkerFaceColor',C);
xlabel('% SiO_2');
ylabel('% TiO_2');
xlim(sax);
ylim([0 7]);
set(gca,'Box','on');

subplot(332); hold on;
plot(data.sio2,data.al2o3, ...
    'ko','MarkerSize',4,'MarkerFaceColor',C);
xlabel('% SiO_2');
ylabel('% Al_2O_3');
xlim(sax);
ylim([0 30]);
set(gca,'Box','on');

subplot(333); hold on;
plot(data.sio2,data.feo_tot, ...
    'ko','MarkerSize',4,'MarkerFaceColor',C);
xlabel('% SiO_2');
ylabel('% FeO');
xlim(sax);
ylim([0 50]);
set(gca,'Box','on');

subplot(334); hold on;
plot(data.sio2,data.mgo, ...
    'ko','MarkerSize',4,'MarkerFaceColor',C);
xlabel('% SiO_2');
ylabel('% MgO');
xlim(sax);
ylim([0 60]);
set(gca,'Box','on');

subplot(335); hold on;
plot(data.sio2,data.cao, ...
    'ko','MarkerSize',4,'MarkerFaceColor',C);
xlabel('% SiO_2');
ylabel('% CaO');
xlim(sax);
ylim([0 100]);
set(gca,'Box','on');

subplot(336); hold on;
plot(data.sio2,data.mno, ...
    'ko','MarkerSize',4,'MarkerFaceColor',C);
xlabel('% SiO_2');
ylabel('% MnO');
xlim(sax);
ylim([0 1]);
set(gca,'Box','on');

subplot(337); hold on;
plot(data.sio2,data.na2o, ...
    'ko','MarkerSize',4,'MarkerFaceColor',C);
xlabel('% SiO_2');
ylabel('% Na_2O');
xlim(sax);
ylim([0 15]);
set(gca,'Box','on');

subplot(338); hold on;
plot(data.sio2,data.k2o, ...
    'ko','MarkerSize',4,'MarkerFaceColor',C);
xlabel('% SiO_2');
ylabel('% K_2O');
xlim(sax);
ylim([0 15]);
set(gca,'Box','on');

subplot(339); hold on;
plot(data.sio2,data.p2o5, ...
    'ko','MarkerSize',4,'MarkerFaceColor',C);
xlabel('% SiO_2');
ylabel('% P_2O_5');
xlim(sax);
ylim([0 3]);
set(gca,'Box','on');

return


function hshist(data,C);

sax = [0 100];
sbin = [1:2:99];

subplot(331); hold on;
yax = [0 7];
dy = [yax(2) - yax(1)]/50;
n = hist3([data.sio2,data.tio2],{sbin [yax(1)+dy/2:dy:yax(2)+dy/2]});
imagesc(sbin,[yax(1)+dy/2:dy:yax(2)-dy/2],log10(n(:,1:end-1))');
hold on;
fill([0 100 100 0],[100 0 100 100],0.7);
xlabel('% SiO_2');
ylabel('% TiO_2');
xlim(sax);
ylim(yax);
set(gca,'Box','on');
shading flat;
caxis([-0.1 3]);

subplot(332); hold on;
yax = [0 30];
dy = [yax(2) - yax(1)]/50;
n = hist3([data.sio2,data.al2o3],{sbin [yax(1)+dy/2:dy:yax(2)+dy/2]});
imagesc(sbin,[yax(1)+dy/2:dy:yax(2)-dy/2],log10(n(:,1:end-1))');
hold on;
fill([0 100 100 0],[100 0 100 100],0.7);
xlabel('% SiO_2');
ylabel('% Al_2O_3');
xlim(sax);
ylim(yax);
set(gca,'Box','on');
shading flat;
caxis([-0.1 3]);

subplot(333); hold on;
yax = [0 50];
dy = [yax(2) - yax(1)]/50;
n = hist3([data.sio2,data.feo_tot],{sbin [yax(1)+dy/2:dy:yax(2)+dy/2]});
imagesc(sbin,[yax(1)+dy/2:dy:yax(2)-dy/2],log10(n(:,1:end-1))');
hold on;
fill([0 100 100 0],[100 0 100 100],0.7);
xlabel('% SiO_2');
ylabel('% FeO');
xlim(sax);
ylim(yax);
set(gca,'Box','on');
shading flat;
caxis([-0.1 3]);

subplot(334); hold on;
yax = [0 60];
dy = [yax(2) - yax(1)]/50;
n = hist3([data.sio2,data.mgo],{sbin [yax(1)+dy/2:dy:yax(2)+dy/2]});
imagesc(sbin,[yax(1)+dy/2:dy:yax(2)-dy/2],log10(n(:,1:end-1))');
hold on;
fill([0 100 100 0],[100 0 100 100],0.7);
xlabel('% SiO_2');
ylabel('% MgO');
xlim(sax);
ylim([0 60]);
ylim(yax);
set(gca,'Box','on');
shading flat;
caxis([-0.1 3]);

subplot(335); hold on;
yax = [0 100];
dy = [yax(2) - yax(1)]/50;
n = hist3([data.sio2,data.cao],{sbin [yax(1)+dy/2:dy:yax(2)+dy/2]});
imagesc(sbin,[yax(1)+dy/2:dy:yax(2)-dy/2],log10(n(:,1:end-1))');
hold on;
fill([0 100 100 0],[100 0 100 100],0.7);
xlabel('% SiO_2');
ylabel('% CaO');
xlim(sax);
ylim(yax);
set(gca,'Box','on');
shading flat;
caxis([-0.1 3]);

subplot(336); hold on;
yax = [0 1];
dy = [yax(2) - yax(1)]/50;
n = hist3([data.sio2,data.mno],{sbin [yax(1)+dy/2:dy:yax(2)+dy/2]});
imagesc(sbin,[yax(1)+dy/2:dy:yax(2)-dy/2],log10(n(:,1:end-1))');
hold on;
fill([0 100 100 0],[100 0 100 100],0.7);
xlabel('% SiO_2');
ylabel('% MnO');
xlim(sax);
ylim([0 1]);
set(gca,'Box','on');
shading flat;
caxis([-0.1 3]);

subplot(337); hold on;
yax = [0 15];
dy = [yax(2) - yax(1)]/50;
n = hist3([data.sio2,data.na2o],{sbin [yax(1)+dy/2:dy:yax(2)+dy/2]});
imagesc(sbin,[yax(1)+dy/2:dy:yax(2)-dy/2],log10(n(:,1:end-1))');
hold on;
fill([0 100 100 0],[100 0 100 100],0.7);
xlabel('% SiO_2');
ylabel('% Na_2O');
xlim(sax);
ylim([0 15]);
ylim(yax);
set(gca,'Box','on');
shading flat;
caxis([-0.1 3]);

subplot(338); hold on;
yax = [0 15];
dy = [yax(2) - yax(1)]/50;
n = hist3([data.sio2,data.k2o],{sbin [yax(1)+dy/2:dy:yax(2)+dy/2]});
imagesc(sbin,[yax(1)+dy/2:dy:yax(2)-dy/2],log10(n(:,1:end-1))');
hold on;
fill([0 100 100 0],[100 0 100 100],0.7);
xlabel('% SiO_2');
ylabel('% K_2O');
xlim(sax);
ylim([0 15]);
ylim(yax);
set(gca,'Box','on');
shading flat;
caxis([-0.1 3]);

subplot(339); hold on;
yax = [0 3];
dy = [yax(2) - yax(1)]/50;
n = hist3([data.sio2,data.p2o5],{sbin [yax(1)+dy/2:dy:yax(2)+dy/2]});
imagesc(sbin,[yax(1)+dy/2:dy:yax(2)-dy/2],log10(n(:,1:end-1))');
hold on;
fill([0 100 100 0],[100 0 100 100],0.7);
xlabel('% SiO_2');
ylabel('% P_2O_5');
xlim(sax);
ylim(yax);
set(gca,'Box','on');
shading flat;
caxis([-0.1 3]);
cbar('Log (frequency)');

colormap(flipud(gray));

return
