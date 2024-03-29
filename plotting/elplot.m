function elplot(data)

fields = {'F','Cl','Br','I','H','C','N', ...
    'S','As','B','Be','Bi','In','Sb', ...
    'Li','Cs','Rb','Ba','Cd', ...
    'Mn','Co','Zn','Cu','Ni', ...
    'V','Cr', ...
    'Ga','Ge','Hg', ...
    'Mo','Pa', ...
    'Pb','Pm', ...
    'Sc','Se', ...
    'Sn','Sr','Te', ...
    'Tl','W','Y', ...
    'Th','U','Zr','Hf','Nb','Ta', ...
    'La','Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu' ...
    'Ag','Au','Pt','Pd','Ir','Os','Re','Rh','Ru'};

for i = 1:length(fields)
    n(i) = sum(data{:,{[lower(fields{i}),'_ppm']}} > 0);
end

figure;
bar(n);
set(gca,'XTick',[1:length(fields)],'XTickLabel',fields);
xlim([0 length(fields)+1]);

return