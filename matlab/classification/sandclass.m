function sed_name = sandclass(data,varargin);

% Errorlog:
%   26 Feb 2019 - fixed an error in feo_tot to fe2o3_tot conversion

plot_flag = 0;
if nargin == 2
    plot_flag = varargin{1};
end

sedgons2 = load_sedgons2;

cf = molecularwt('Fe2O3')/(2*molecularwt('FeO'));
fe2o3_tot = data.feo_tot*cf;

ind = find(data.sio2 > 0 & data.al2o3 > 0 & fe2o3_tot > 0 & data.k2o > 0 & ...
    ~strcmp('limestone',data.rock_type) & ...
    ~strcmp('dolomite',data.rock_type) & ...
    ~strcmp('quartzite',data.rock_type) & ...
    ~strcmp('oxide',data.rock_type) & ...
    ~strcmp('laterite',data.rock_type) & ...
    ~strcmp('bauxite',data.rock_type));

if plot_flag
    figure; hold on;
    
    SA = log10(data.sio2(ind)./data.al2o3(ind));
    FK = log10(fe2o3_tot(ind)./data.k2o(ind));
    eSA = [-0.5:0.02:3];
    eFK = [-2:0.05:5];
    
    n = hist2d(SA,FK,eSA,eFK);
    imagesc(eSA,eFK,log10(n));
    colorbar;
    caxis([-0.1 3]);
    cbar('No. data');
    axis xy;
    xlim([min(eSA) max(eSA)]);
    ylim([min(eFK) max(eFK)]);
    
    xlabel('log_{10}(SiO_2/Al_2O_3)');
    ylabel('log_{10}(Fe_2O_3/K_2O)');
    
    for i = 1:length(sedgons2)
        plot(sedgons2{i,1}(:,1),sedgons2{i,1}(:,2),'w-');
    end
end

[nsed,~] = size(sedgons2);
sed_name = data.rock_type;
%plot(log10(data.sio2(ind)./data.al2o3(ind)),log10(fe2o3_tot(ind)./data.k2o(ind)),'.');
for i = 1:nsed
    in = inpolygon(log10(data.sio2(ind)./data.al2o3(ind)),log10(fe2o3_tot(ind)./data.k2o(ind)), ...
        sedgons2{i,1}(:,1),sedgons2{i,1}(:,2));
    
    %plot(log10(data.sio2(ind(in))./data.al2o3(ind(in))),log10(fe2o3_tot(ind(in))./data.k2o(ind(in))),'.');
    
    sed_name(ind(in)) = {sedgons2{i,2}};
end
data.rock_type = sed_name;

%ind = strcmp('pelite',data.rock_type) & data.sio2 > 0 & data.al2o3 > 0 & fe2o3_tot > 0 & data.k2o > 0;
%plot(log10(data.sio2(ind)./data.al2o3(ind)),log10(fe2o3_tot(ind)./data.k2o(ind)),'k.');

rt = unique(data.rock_type);
if plot_flag
    figure;
    ternary('Q','F','L');
    hold on;
    for i = 1:length(rt)
        ind = strcmp(rt{i},data.rock_type);
        ternplot(data.quartz(ind),data.feldspar(ind),data.lithics(ind),'.');
    end

    sedgons = load_sedgons;
    for i = 1:length(sedgons)
        ternplot(sedgons{i,1}(:,1),sedgons{i,1}(:,2),sedgons{i,1}(:,3),'-');
    end
end

return