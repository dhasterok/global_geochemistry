function data = sedclass(data,varargin);

% Errorlog:
%   26 Feb 2019 - fixed an error in feo_tot to fe2o3_tot conversion

plot_flag = 0;
if nargin == 2
    plot_flag = varargin{1};
end

%ind = data.sio2 >= 0 & data.feo_tot >= 0 & data.al2o3 >= 0 & data.cao >= 0 & data.mgo >= 0;
%data = data(ind,:);

cf = molecularwt('Fe2O3')/(2*molecularwt('FeO'));
fe2o3_tot = data.feo_tot*cf;

% quartz
Q = data.sio2/molecularwt('SiO2');
Q(Q < 0) = NaN;

% feldspar/clay/Fe-Al oxide
F = data.al2o3/molecularwt('Al2O3') + ...
    fe2o3_tot/molecularwt('Fe2O3');
F(data.al2o3 < 0 | fe2o3_tot < 0) = NaN;

% lithics (carbonates)
L = data.mgo/molecularwt('MgO') + data.cao/molecularwt('CaO');
L(data.mgo < 0 | data.cao < 0) = NaN;

T = Q + F + L;
Q = Q./T*100;
F = F./T*100;
L = L./T*100;

data.quartz = Q;
data.feldspar = F;
data.lithics = L;

% load QFL polygons
sedgons = load_sedgons;

[nsed,~] = size(sedgons);
sed_name = data.rock_type;
for i = 1:nsed
    in = inpolygon(Q,F,sedgons{i,1}(:,1),sedgons{i,1}(:,2));

    switch sedgons{i,2}
        case 'carbonate'
            ind = in & 0.4 < data.mgo/molecularwt('MgO') ...
                ./(data.mgo/molecularwt('MgO') + data.cao/molecularwt('CaO'));
            sed_name(ind) = {'dolomite'};
            ind = in & 0.4 >= data.mgo/molecularwt('MgO') ...
                ./(data.mgo/molecularwt('MgO') + data.cao/molecularwt('CaO'));
            sed_name(ind) = {'limestone'};
        case 'laterite'
            ind = in & ...
                data.al2o3/molecularwt('Al2O3') <= fe2o3_tot/molecularwt('Fe2O3');
            sed_name(ind) = {'laterite'};
            ind = in & ...
                data.al2o3/molecularwt('Al2O3') > fe2o3_tot/molecularwt('Fe2O3');
            sed_name(ind) = {'bauxite'};
        otherwise
            sed_name(in) = {sedgons{i,2}};
    end
end
data.rock_type = sed_name;

if plot_flag
    data.rock_type = sandclass(data,1);
else
    data.rock_type = sandclass(data);
end
%figure; hold on;
%ternary('SiO_2','(Al,Fe)_2O_3','(Ca,Mg)O');
%ternplot(Q,F,L,'.');
%for i = 1:nsed
%    ternplot(sedgons{i,1}(:,1),sedgons{i,1}(:,2),sedgons{i,1}(:,3),'-');
%end

if plot_flag
    figure;
    hold on;

    out = ternsurf(Q,F,L,ones(size(Q)),0.05);
    
    ternary('SiO_2','(Al,Fe)_2O_3','(Ca,Mg)O');
    trisurf(out.tri,out.xv,out.yv,out.nbin);
    %shading interp;
    shading flat;
    set(gca,'View',[0 90]);
    hold on;
    caxis([0 3]);
    cbar('No. sed.');

    for i = 1:nsed
        ternplot(sedgons{i,1}(:,1),sedgons{i,1}(:,2),sedgons{i,1}(:,3),'-');
    end
end

return
