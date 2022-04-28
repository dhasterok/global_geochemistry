function tas_name = tas_sw(data,varargin);
% TAS - Total alkali vs. silica determination.
%
%   data = tas(data) determines the compositionally based rock names for
%   igneous rocks based on the total alkai silica scheme by Middlemost,
%   ESR, 1994.  Plutonic and volcanic rocks are estimated separately.
%   
%   data = tas(data,'Plot',plottype) where plottype is a string either
%   'hist' for a 2-D histogram or as 'scatter' producing an x-y scatter
%   plot, both with TAS fields.  For no plot (default), plottype = 'none'.
%
%   High-Mg volcanic rocks are further classified using the scheme by
%   Le Bas & Streckeisen, J. Geol. Soc. London, 1991.  Additional fields
%   include meimichite, komatiite, picrite, alkali picrite and boninite.
%
%   Ultramafic proutonic rocks are separated by Mg number into mantle (>80)
%   and cumulate peridotites (<80).
%
%   Carbonatites are specifically identified as such using given names
%   or high-CO2 content >20%.
%
%   tas(data,plot_flag) produces a 2D histogram of TA vs. S with negative
%   values indicating high-Mg lavas.

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

tas_name = cell(size(data.sio2));   % TAS determined name
tas_name(:) = {''};

% load TAS and UH-Mg fields
[tasgons,uhmgons] = load_tasgons;

[ntas,~] = size(tasgons);
[nuhm,~] = size(uhmgons);

volc = zeros([height(data) 1]);
volc(rockgroup(data,'all volcanic')) = 1;
volc(rockgroup(data,'all plutonic')) = 2;

TA = data.na2o + data.k2o;
S = data.sio2;
MgO = data.mgo;
TiO2 = data.tio2;

%Mg_number = mgnum(data.mgo,data.feo_tot,0);

% fix polygons for type of peridotite based on TiO2
for i = 1:ntas
    in = inpolygon(S,TA,tasgons{i,1}(:,1),tasgons{i,1}(:,2));

    ind = in & volc == 1;
    tas_name(ind) = {tasgons{i,2}};

    ind = in & volc == 2;
    tas_name(ind) = {tasgons{i,3}};

    % classification of ultramafic vs cumulate for plutonic rocks
    % Guided by Reverdatto et al. (Russ. Geol. & Geophys., 2008)
    if strcmp(tasgons{i,3},'peridotite') | strcmp(tasgons{i,3},'peridotgabbro')
        tas_name(ind & 0 < TiO2 & TiO2 <= 0.3 & 35 <= MgO) = {'mantle peridotite'};
        tas_name(ind & 0 < TiO2 & TiO2 <= 0.3 & 18 < MgO & MgO < 35) = {'mantle pyroxenite'};
        if ~strcmp(tasgons{i,3},'peridotgabbro')
            tas_name(ind & TiO2 > 0.3 & 22 <= MgO) = {'crustal peridotite'};
            tas_name(ind & TiO2 > 0.3 & 0 < MgO & MgO < 22) = {'crustal pyroxenite'};
        end
    end
end


% classification of high magnesian rocks
for i = 1:nuhm
    in = inpolygon(S,TA,uhmgons{i,1}(:,1),uhmgons{i,1}(:,2));
    
    if ~strcmp(uhmgons{i,2},'boninite')
        ind = in & volc == 1 & MgO >= 18;
        
        tas_name(ind & MgO > 18) = {uhmgons{i,2}};
        if strcmp(uhmgons{i,2},'komatiite')
            tas_name(ind & TiO2 > 1) = {uhmgons{i,3}};
        end
        
        if strcmp(uhmgons{i,2},'picrite') | ...
            strcmp(uhmgons{i,2},'alkali picrite')
            ind = in & volc == 2 & MgO >= 18 & S > 45;
            tas_name(ind) = {uhmgons{i,4}};
        end
    end
    tas_name(ind & i == 3 & TiO2 > 1) = {uhmgons{i,3}};
    
    if strcmp(uhmgons{i,2},'boninite')
        ind = in & volc == 1 & MgO > 8 & TiO2 < 0.5;
        tas_name(ind) = {uhmgons{i,2}};
        
        ind = in & volc == 2 & MgO > 8 & TiO2 < 0.5;
        tas_name(ind) = {uhmgons{i,4}};
    end
end

% assign carbonatite names to high CO2 rocks
carbname = carbclass(data,'Plot',plottype);

ind = ~strcmp(carbname,'');
tas_name(ind) = carbname(ind);



% print table of number of rocks in each field
[type,~,ind] = unique(tas_name);
fprintf('No. Samples         Rock Name (calc.)\n');
fprintf('-----------   ----------------------------\n');
for i = 1:length(type);
    num = length(find(ind == i));
    fprintf('%10i    %s\n',num,type{i});
end


% plots tas fields
if ~strcmpi(plottype,'none')
    figure;
    hold on;

    % 2-d histogram of tas values
    if strcmp(plottype,'hist')
        eS = [0:0.5:100];
        eTA = [-10:0.2:30];
        TA(MgO > 18) = -TA(MgO > 18);
        n = hist2d(S,TA,eS,eTA);
        imagesc(eS,eTA,log10(n));
        colormap([1 1 1; parula]);
        colorbar;
        caxis([-0.1 3]);
        axis xy;
    elseif strcmp(plottype,'scatter');
        hold on
        ind = strcmpi(data.QAP_name,'monzogranite');
        scatter(S(ind),TA(ind),[],'g','filled');
        
        ind = strcmpi(data.QAP_name,'syenogranite');
        scatter(S(ind),TA(ind),[],'b','filled');
        
        ind = strcmpi(data.QAP_name,'granodiorite');
        scatter(S(ind),TA(ind),[],'r','filled');
        
        ind = strcmpi(data.QAP_name,'tonalite');
        scatter(S(ind),TA(ind),[],'k','filled');
    else
        warning('Unknown plot type.');
        return;
    end

    for i = 1:ntas
        x = tasgons{i,1}(:,1);
        y = tasgons{i,1}(:,2);
        if strcmp(plottype,'hist')
            plot(x,y,'m-');
        else
            plot(x,y,'k-');
        end
        %fill(x,y,[0.7 0.7 0.7]);

        if length(volc) == length(find(volc == 1))
            text(mean(x),mean(y),tasgons{i,2});
            title('volcanic');
        elseif length(volc) == length(find(volc == 2))
            text(mean(x),mean(y),tasgons{i,3});
            title('plutonic');
        else
            text(mean(x),mean(y),tasgons{i,2});
            title('other');
        end
    end
    % plots uhmg fields
    for i = 1:nuhm
        x = uhmgons{i,1}(:,1);
        y = -uhmgons{i,1}(:,2);
        if strcmp(plottype,'hist')
            plot(x,y,'m-');
        else
            plot(x,y,'k-');
        end
        %fill(x,y,[0.7 0.7 0.7]);
        if isempty(find(volc ~= 2))
            text(mean(x),mean(y),uhmgons{i,2});
        else
            text(mean(x),mean(y),uhmgons{i,2});
        end
    end
    % plots TAS points from geochemistry
    %plot(S,TA,'.')
    xlabel('sio2 [wt.%]');
    ylabel('k2o + na2o [wt.%]');
    xlim([0 100]);
    ylim([-10 30]);
    set(gca,'Box','on');
    golden;
end

return
