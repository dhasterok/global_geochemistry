function hpvelhist(data)

dVp = 0.05;
dRho = 20;
Vp = [5.8:0.2:8.4];
Rho = [2500:50:3400];
SiO2 = [35:5:85];

f2 = figure;
%[m,merr,C] = hpVplot(data.p_velocity,data.heat_production,Vp,[0 0.447 0.741]);
%title(['log_{10}[A (\muW m^{-3})] = ',num2str(m(1)),'\pm',num2str(merr(1)), ...
%    '[V_P - 6 (km s^{-1})] + ',num2str(m(2)),'\pm',num2str(merr(2)),', r^2 = ',num2str(C^2)]);
[m,merr,C] = hpVplot2(data.p_velocity, data.heat_production, dVp, ...
    'XLim',[5.8 8.4], 'XFit',[6.0 7.4]);

f3 = figure;
%[m,merr,C] = hprhoplot(data.density_model,data.heat_production,Rho,[0 0.447 0.741]);
%title(['log_{10}[A (\muW m^{-3})] = ',num2str(m(1)),'\pm',num2str(merr(1)), ...
%    '[Rho - 2700 (kg m^{-3})] + ',num2str(m(2)),'\pm',num2str(merr(2)),', r^2 = ',num2str(C^2)]);
[m,merr,C] = hpVplot2(data.density_model, data.heat_production, dRho, ...
    'XShift',0, 'XLim',[2500 3400], 'XLabel','Density [kg m^{-3}]', ...
    'XFit',[2650 2950]);

return
