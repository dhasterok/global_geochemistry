function carb_name = carbclass(data,varargin);

plottype = 'none';
opt = 1;
while opt < nargin
    switch lower(varargin{opt})
    case 'plot'
        plottype = lower(varargin{opt+1});
        % plottypes:
        % 'hist' - 2-D histogram
        % 'scatter' - points
    end
    opt = opt + 2;
end

carb_name = cell(size(data.sio2));   % TAS determined name
carb_name(:) = {''};

if any(strcmp(data.Properties.VariableNames,'co2'))
    carbind = data.co2 >= 20 | strcmpi(data.rock_name,'carbonatite') | abs(data.sio2) < 33;
    carb_name(carbind) = {'carbonatite'};
else
    carbind = strcmpi(data.rock_name,'carbonatite') | abs(data.sio2) < 33;
    carb_name(carbind) = {'carbonatite'};
end

if ~any(carbind)
    return
end

% silicocarbonatite
ind = carbind & 20 < data.sio2 & data.sio2 < 33;
carb_name(ind) = {'silicocarbonatite'};

C = data.cao;
M = data.mgo;
F = 0.8*data.feo_tot + ...
    0.2*data.feo_tot*molecularwt('Fe2O3')/(2*molecularwt('FeO')) + ...
    data.mno;
ind = data.fe2_fe_tot > 0;
F(ind) = data.fe2_fe_tot(ind).*data.feo_tot(ind) + ...
    (1 - data.fe2_fe_tot(ind)).*data.feo_tot(ind)*molecularwt('Fe2O3')/(2*molecularwt('FeO')) + ...
    data.mno(ind);

T = C + M + F;
C = C./T;
M = M./T;
F = F./T;

carbgons = load_carbgons;

[X,Y] = tern2xy(C,M,F);
for i = 1:length(carbgons(:,1))
    c = carbgons{i,1}(:,1);
    m = carbgons{i,1}(:,2);
    f = carbgons{i,1}(:,3);
    if ~strcmp(plottype,'none')
        if i == 1
            figure;
        end
        ternplot(c,m,f,'-k');
    end
    
    [x,y] = tern2xy(c,m,f);
    in = inpolygon(X,Y,x,y);
    
    ind = in & carbind;
    carb_name(ind) = {carbgons{i,2}};
    carb_name(ind & data.sio2 > 20) = {['silico-',carbgons{i,2}]};
end

scarbind = carbind & data.sio2 > 20 & ~isnan(C + M + F);
carbind = carbind & ~(data.sio2 > 20) & ~isnan(C + M + F);
if ~strcmpi(plottype,'none')

    ternary('CaO','MgO','FeO + Fe2O3 + MnO');
    hold on;
    
    if strcmp(plottype,'hist')
        ternsurf(C(carbind | scarbind),M(carbind | scarbind),F(carbind | scarbind),data.sio2(carbind | scarbind),0.05);
        hold on;
    elseif strcmp(plottype,'scatter');
        %[C(carbind),M(carbind),F(carbind)]
        ternplot(C(carbind),M(carbind),F(carbind),'.');
        ternplot(C(scarbind),M(scarbind),F(scarbind),'.');
    else
        warning('Unknown plot type.');
        return;
    end
    
    for i = 1:length(carbgons(:,1))
        c = carbgons{i,1}(:,1);
        m = carbgons{i,1}(:,2);
        f = carbgons{i,1}(:,3);

        ternplot(c,m,f,'-k');
        [x,y] = tern2xy(c,m,f);

        t = text(mean(x),mean(y),carbgons{i,2});
        set(t,'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

return