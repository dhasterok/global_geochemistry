function hp_spatial(data,locales)
% HP_SPATIAL - spatial filtered heat production vs. P-velocity
%
%   hp_spatial(data,locales)
%
%   DATA is the input structure containing heat production and velocity
%   estimates and LOCALES is a cell array containing the regions to be plotted.
%   If locales contains a string, then it looks for the region identified in
%   data.country.  Special locales exist for 'global', 'continent', and 'ocean'.
%   Also continents 'NAM'.  If locales is an array of [lat,lon] pairs then
%   hp_spatial will find all the points which lie within the polygon.

if nargin ~= 2
    help hp_spatial
    return;
end

Vp = [5.8:0.1:8.2];
Rho = [2500:50:3400];
colour = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

n = length(locales);

f1 = figure;
f2 = figure;
for i = 1:n
    switch locales{i}
    case 'global'
        ind = ~isnan(data.heat_production);
    case 'NAM'
        ind = ((strcmp('US',data.country) ...
        | strcmp('CA',data.country) ...
        | strcmp('MX',data.country)) ...
        & ~isnan(data.heat_production));
    case 'continent'
        ind = (~strcmp('ocean',data.country) ...
        & ~isempty(data.country) ...
        & ~isnan(data.heat_production));
    otherwise
        ind = (strcmp(locales{i},data.country) & ~isnan(data.heat_production));
    end

    A{i} = data.heat_production(ind);
    V{i} = data.p_velocity(ind);
    rho{i} = data.density_model(ind);
    lat{i} = data.latitude(ind);
    lon{i} = data.longitude(ind);

    figure(f1);
    %[m{i},merr{i},stats{i}] = hpVplot2(V{i}, A{i}, 0.1, ...
    %    'Color',colour(i,:), 'XLim',[5.8 7.6], 'XFit',[6 7.4]);
    %subplot(3,1,1:2);
    [m{i},merr{i},C{i}] = hpVplot(V{i},A{i},Vp,colour(i,:));
    xlim([5.8 7.6]);

    figure(f2);
    [mr{i},mrerr{i},statsr{i}] = hprhoplot(rho{i},A{i},Rho,colour(i,:));
    %hpVhist(V{i},A{i},Vp,colour(i,:));
end


figure(f1);
%subplot(3,1,3); hold on;
legend(locales);

fprintf('\n\n');
for i = 1:n
    %fprintf('%s %i %f+/-%f, %f+/-%f, %f\n', ...
    %    locales{i},length(A{i}),m{i}(1),merr{i}(1),m{i}(2),merr{i}(2),stats{i}(1));
    fprintf('%s %i %f+/-%f, %f+/-%f, %f\n', ...
        locales{i},length(A{i}),m{i}(1),merr{i}(1),m{i}(2),merr{i}(2),C{i}.^2);
    
    fid1 = fopen(['outfiles/',locales{i},'.llv'],'w');
    fid2 = fopen(['outfiles/',locales{i},'.lla'],'w');
    for j = 1:length(A{i})
        fprintf(fid1,'%f %f %f\n',lat{i}(j),lon{i}(j),V{i}(j));
        fprintf(fid2,'%f %f %f\n',lat{i}(j),lon{i}(j),0.01125*log10(A{i}(j)) + 0.04375);
    end
    fclose(fid1);
    fclose(fid2);
end
fprintf('\n');
for i = 1:n
    fprintf('%s %i %f+/-%f %f+/-%f\n', ...
        locales{i},length(A{i}),mr{i}(1),mrerr{i}(1),mr{i}(2),mrerr{i}(2));
end
fprintf('\n\n');

return
