function mmchemdd(data)
% MMCHEMDD - creates chemical plots used for protolith discrimination
%
%   mmchemdd(data) where data is a table of geochemistry.

iind = strcmp(data.rock_group,'igneous');
sind = strcmp(data.rock_group,'sedimentary');

data.type = iind;

% Ref: Winkler, Petrogenesis of Metamorphic rocks, 1979.
% ACF diagram
mclass(data(iind,:),{'Al2O3-Na2O-K2O','CaO-10/3*P2O5','FeO+MgO-TiO2'},1,'sio2');
colormap(flipud(bone));
drawnow;
mclass(data(sind,:),{'Al2O3-Na2O-K2O','CaO-10/3*P2O5','FeO+MgO-TiO2'},1,'sio2');
colormap(flipud(bone));
drawnow;

% Ref: Winkler, Petrogenesis of Metamorphic rocks, 1979.
% % A-CN-K diagram
% mclass(data(iind,:),{'Al2O3','CaO+Na2O','K2O'},1,'sio2');
% colormap(flipud(bone));
% drawnow;
% mclass(data(sind,:),{'Al2O3','CaO+Na2O','K2O'},1,'sio2');
% colormap(flipud(bone));
% drawnow;
% 
% % FeOT-MgO-CaO diagram
% mclass(data(iind,:),{'FeO','MgO','CaO'},1,'sio2');
% colormap(flipud(bone));
% drawnow;
% mclass(data(sind,:),{'FeO','MgO','CaO'},1,'sio2');
% colormap(flipud(bone));
% drawnow;

% figure;
% sio2 = [0:1:100];
% tio2 = [0:0.1:10];
% 
% subplot(121);
% [ni,out,CX,CY] = hist2d(data.sio2(iind),data.tio2(iind),sio2,tio2);
% imagesc(CX,CY,log10(ni));
% colorbar;
% hold on;
% plot(sio2,1.7/-42*(sio2 - 40) + 1.7,'-');
% plot([0 100],[2.5 2.5],'-');
% xlabel('SiO2 (wt.%)');
% ylabel('TiO2 (wt.%)');
% golden
% axis([0 100 0 8]);
% axis xy;
% caxis([0 3]);
% 
% subplot(122);
% [ns,out,CX,CY] = hist2d(data.sio2(sind),data.tio2(sind),sio2,tio2);
% imagesc(CX,CY,log10(ns));
% colorbar
% colormap(flipud(bone));
% hold on;
% plot(sio2,1.7/-42*(sio2 - 40) + 1.7,'-');
% plot([0 100],[2.5 2.5],'-');
% xlabel('SiO2 (wt.%)');
% ylabel('TiO2 (wt.%)');
% golden
% axis([0 100 0 8]);
% axis xy;
% caxis([0 3]);
% 
% fprintf('TiO2 tests\n');
% fprintf('Full vs. empty cell comparison\n');
% ui = ni(:) > 0 & ns(:) == 0;
% us = ns(:) > 0 & ni(:) == 0;
% 
% 
% fprintf('Ig.  %.1f%%\n',sum(sum(ni(ui)))/sum(iind)*100);
% fprintf('Sed. %.1f%%\n',sum(sum(ns(us)))/sum(sind)*100);
% 
% fprintf('Tarney 1979 SiO2-TiO2 protolith test\n');
% fprintf('Ig. TP   %.1f%%\n',sum(1.7/-42*(data.sio2(iind) - 40) + 1.7 > data.tio2(iind))/sum(iind)*100)
% fprintf('Ig. FN   %.1f%%\n',sum(1.7/-42*(data.sio2(iind) - 40) + 1.7 < data.tio2(iind))/sum(iind)*100)
% fprintf('Sed. TP  %.1f%%\n',sum(1.7/-42*(data.sio2(sind) - 40) + 1.7 < data.tio2(sind))/sum(sind)*100)
% fprintf('Sed. FN  %.1f%%\n',sum(1.7/-42*(data.sio2(sind) - 40) + 1.7 > data.tio2(sind))/sum(sind)*100)
% 
% fprintf('Simple >2.5 wt.%% TiO2 test\n');
% fprintf('Probably Ig. TP   %.1f%%\n',sum(data.tio2(iind) > 2.5)/sum(iind)*100);
% fprintf('Probably Ig. FP  %.1f%%\n',sum(data.tio2(sind) > 2.5)/sum(sind)*100);

data.al = data.al2o3/molecularwt('Al2O3');
data.c = data.cao/molecularwt('CaO');
data.alk = data.na2o/molecularwt('Na2O') + data.k2o/molecularwt('K2O');
data.fm = data.feo_tot/molecularwt('FeO') + data.mgo/molecularwt('MgO');
n = data.al + data.c + data.alk + data.fm;
data.al = 100*data.al./n;
data.c = 100*data.c./n;
data.alk = 100*data.alk./n;
data.fm = 100*data.fm./n;
data.si = 100*data.sio2/molecularwt('SiO2')./n;

% Niggli igneous field digitized from Li et al., (2018)
% doi:10.1016/j.precamres.2018.04.018
ni_ig = load('nigili_igneous.xy');
ni_ig = [ni_ig; ni_ig(1,:)];

figure;
S = [0:6:600];
NI = [-100:2:100];
subplot(121);
[ni,out,CX,CY] = hist2d(data.si(iind),data.al(iind)+data.fm(iind)-data.c(iind)-data.alk(iind),S,NI);
imagesc(CX,CY,log10(ni));
colorbar
colormap(flipud(bone));
hold on;
plot(ni_ig(:,1),ni_ig(:,2),'-');

xlabel('si');
ylabel('(al + fm) - (c + alk)');
golden
axis([0 600 -10 100]);
axis xy;
caxis([0 3]);

subplot(122);
[ns,out,CX,CY] = hist2d(data.si(sind),data.al(sind)+data.fm(sind)-data.c(sind)-data.alk(sind),S,NI);
imagesc(CX,CY,log10(ns));
colorbar
colormap(flipud(bone));
hold on;
plot(ni_ig(:,1),ni_ig(:,2),'-');

xlabel('si');
ylabel('(al + fm) - (c + alk)');
golden
axis([0 600 -10 100]);
axis xy;
caxis([0 3]);

fprintf('Niggli test\n');
fprintf('Full vs. empty cell comparison\n');
ui = ni(:) > 0 & ns(:) == 0;
us = ns(:) > 0 & ni(:) == 0;

fprintf('Ig.  %.1f%%\n',sum(sum(ni(ui)))/sum(iind)*100);
fprintf('Sed. %.1f%%\n',sum(sum(ns(us)))/sum(sind)*100);

fprintf('Ig. TP   %.1f\n',sum(inpolygon(data.si(iind),data.al(iind)+data.fm(iind)-data.c(iind)-data.alk(iind),ni_ig(:,1),ni_ig(:,2)))/sum(iind)*100);
fprintf('Ig. FN   %.1f\n',sum(~inpolygon(data.si(iind),data.al(iind)+data.fm(iind)-data.c(iind)-data.alk(iind),ni_ig(:,1),ni_ig(:,2)))/sum(iind)*100);
fprintf('Sed. TP  %.1f\n',sum(~inpolygon(data.si(sind),data.al(sind)+data.fm(sind)-data.c(sind)-data.alk(sind),ni_ig(:,1),ni_ig(:,2)))/sum(sind)*100);
fprintf('Sed. FN  %.1f\n',sum(inpolygon(data.si(sind),data.al(sind)+data.fm(sind)-data.c(sind)-data.alk(sind),ni_ig(:,1),ni_ig(:,2)))/sum(sind)*100);


return