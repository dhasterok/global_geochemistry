close all;
clear all;

addpath ../maptools
addpath ../worldgrid
addpath ../plotting
addpath ../toolbox

eqmap = readtable('eq_depth_palette.csv');
volcmap = [0.8275    0.4745    0.6392; % interplate
    0.8745    0.7412    0.3686; % rifting
    0.2353    0.8275    0.9176; % subduction zone
    0.5882    0.5922    0.6275]; % unknown


eq = readtable('~/toolbox/geophysics/seismic/earthquakes/ncedc_3.0_1980_2015.csv');

eq.cdata = nan(height(eq),3);
for i = 1:height(eqmap)
    ind = eqmap.vmin(i) <= eq.Depth & eq.Depth < eqmap.vmax(i);
    eq.cdata(ind,1) = eqmap{i,'r'}/255;
    eq.cdata(ind,2) = eqmap{i,'g'}/255;
    eq.cdata(ind,3) = eqmap{i,'b'}/255;
end

s = shaperead('../../GIS/global_tectonics/point_data/GVP_Holocene+Pleistocene.shp');
volc = struct2table(s);
[~,~,ia] = unique(volc.TectonicS);
for i = 1:4
    volc.cdata(ia == i,1) = volcmap(i,1);
    volc.cdata(ia == i,2) = volcmap(i,2);
    volc.cdata(ia == i,3) = volcmap(i,3);
end
indh = strcmp(volc.Epoch,'Holocene');

sp = readtable('subduction_profiles.csv');
r = 100;

% load seismic velocity
%filename = '/Users/derrickhasterok/toolbox/geophysics/seismic/tomography/3D2018-08Sv-depth.nc';
%filename = '../seismic/tomography/TX2019slab_percent.nc';
%vslon = double(ncread(filename,'longitude'));
%vslat = double(ncread(filename,'latitude'));
%depth = double(ncread(filename,'depth'));
%dvs = ncread(filename,'dvs');

filename = '../../data/seismic/tomography/UU-P07_global.nc';
vplon = double(ncread(filename,'lon'));
vplat = flipud(double(ncread(filename,'lat')));
depth = double(ncread(filename,'depth'));
dvp = permute(ncread(filename,'dvp'),[2 1 3]);

% load crustal thickness
filename = '../worldgrid/MohoDepthSzwillus.nc';
mlon = ncread(filename,'lon');
mlat = ncread(filename,'lat'); 
moho = ncread(filename,'crustalthickness')';

%fig1 = figure;
%plotcoast;
%scatter(eq.Longitude,eq.Latitude,eq.Magnitude.^2,eq.cdata,'o','filled');

for i = 1:height(sp)
    % load elevation surrounding profile
    latlim = sort([sp.lat1(i) sp.lat2(i)]) + 5*[-1 1];
    if abs(sp.lon1(i) - sp.lon2(i)) > 180
        lonlim = sort([sp.lon1(i) sp.lon2(i)]) + 5*[1 -1];
    else
        lonlim = sort([sp.lon1(i) sp.lon2(i)]) + 5*[-1 1];
    end
    lonlim(lonlim < -180) = -180;
    lonlim(lonlim > 180) = 180;

    if min(lonlim) < -150 & max(lonlim) > 150
        [elev1,elat,elon1] = worldgrid(latlim,[max(lonlim) 179.9666],'etopo2');
        [elev2,elat,elon2] = worldgrid(latlim,[-180 min(lonlim)],'etopo2');

        elev = [elev1 elev2];
        elon = [elon1-360; elon2];
    else
        [elev,elat,elon] = worldgrid(latlim,lonlim,'etopo2');
    end

    p1 = [sp.lon1(i) sp.lat1(i)];
    p2 = [sp.lon2(i) sp.lat2(i)];

    % compute elevation along profile
    [delev,Elev] = gctrend(elon,elat,elev,p1,p2);

    % compute moho along profile
    [dmoho,Moho] = gctrend(mlon,mlat,moho,p1,p2,3);
    
    % compute seismic velcity cross section
    %[dvel,Depth,dVs] = gcslice(vslon,vslat,depth,dvs,p1,p2);
    [dvel,Depth,dVp] = gcslice(vplon,vplat,depth,dvp,p1,p2);

    % compute earthquake and volcano locations along profile
    [deq,inde] = gcxsection(eq.Longitude,eq.Latitude,r,p1,p2);
    [dv,indv] = gcxsection(volc.Longitude,volc.Latitude,r,p1,p2);
    
    %figure(fig1);
    %hold on;
    %scatter(eq.Longitude(inde),eq.Latitude(inde),eq.Magnitude(inde),'ro');
    %plot([p1(1) p2(1)],[p1(2) p2(2)],'g-');

    figure;
    %pcolor(dvel,Depth,dVs);
    pcolor(dvel,Depth,dVp);
    shading interp;
    hold on;
    colormap2('rwb');

    plot(delev,-Elev*1e-3,'k-');
    hold on;
    plot(dmoho,Moho,'k-');

    scatter(deq(inde),eq.Depth(inde),eq.Magnitude(inde).^2,eq.cdata(inde,:),'filled');
    
    scatter(dv(indv & indh),zeros(sum(indv & indh),1)-10,[],volc.cdata(indv & indh,:),'^','filled','MarkerEdgeColor',[0 0 0]);
    scatter(dv(indv & ~indh),zeros(sum(indv & ~indh),1)-10,[],volc.cdata(indv & ~indh,:),'v','filled','MarkerEdgeColor',[0 0 0]);
    
    title(sp.name(i));
    axis ij;
    xlabel('Distance (km)');
    ylabel('Depth (km)');
    ylim([-20,700]);
    axis equal;

end