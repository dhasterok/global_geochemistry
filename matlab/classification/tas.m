function tas_name = tas(data,varargin);
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
%   There are a few option value pairs that can be provided to tas:
%
%       'PlotType'      'hist' will produce a 2D histogram of TA vs. S with
%                       negative values indicating high-Mg lavas.
%                       'scatter' will produce a scatter plot.  The default
%                       is 'none'.
%
%       'SizeField'     if 'PlotType' is 'scatter', a field of data can be
%                       selected to set the size of symbols.
%
%       'ColorField'    if 'PlotType' is 'scatter', a field of data can be
%                       selected to set the color of symbols.
%
% See also: CARBCLASS

p = inputParser;
%checkdata = @(data) isa(data,'table');
addRequired(p,'data');
addParameter(p,'PlotType','none',@ischar);
addParameter(p,'SizeField','none',@ischar);
addParameter(p,'ColorField','none',@ischar);

parse(p,data,varargin{:});
p.Results
%data = p.data;
switch p.Results.PlotType
    case {'none','hist','scatter'}
        plottype = p.Results.PlotType;
    otherwise
        error('Unknown plot type.');
end

switch p.Results.SizeField
    case data.Properties.VariableNames
        ptsize = data{:,p.Results.SizeField};
    case 'none'
        ptsize = [];
    otherwise
        error('Unknown field for plotting size.');
end

switch p.Results.ColorField
    case data.Properties.VariableNames
        ptcolor = data{:,p.Results.ColorField};
    case 'none'
        ptcolor = [];
    otherwise
        error('Unknown field for plotting color.');
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

ind = data.na2o > 0 & data.k2o > 0 & data.sio2 > 0 & data.mgo > 0 & data.feo_tot > 0 & data.tio2 > 0;

TA = data.na2o + data.k2o;
S = data.sio2;
MgO = data.mgo;
TiO2 = data.tio2;
FeOT = data.feo_tot;

%Mg_number = mgnum(data.mgo,data.feo_tot,0);

% fix polygons for type of peridotite based on TiO2
for i = 1:ntas
    in = inpolygon(S,TA,tasgons{i,1}(:,1),tasgons{i,1}(:,2));

    ind = in & volc == 1;
    tas_name(ind) = {tasgons{i,2}};

    ind = in & volc == 2;
    tas_name(ind) = {tasgons{i,3}};
end


% classification of high magnesian rocks
in = inpolygon(S,TA,[33 33 65 65 33],[0 100 100 0 0]);
ind = in & S >= 52 & MgO >= 8 & TiO2 < 0.5;
tas_name(ind & volc == 1) = {'boninite'};
tas_name(ind & volc == 2) = {'sanukitoid'};

ind = in & S < 52 & MgO >= 12 & TA < 3 & FeOT >= 13;
tas_name(ind) = {'ferropicrite'};

ind = in & S < 52 & MgO >= 12 & TA < 3 & FeOT < 13;
tas_name(ind) = {'picrite'};

ind = in & S < 52 & MgO >= 12 & TA >= 3 & FeOT >= 13;
tas_name(ind) = {'alkali ferropicrite'};

ind = in & S < 52 & MgO >= 12 & TA >= 3 & FeOT < 13;
tas_name(ind) = {'alkali picrite'};

ind = in & S < 52 & MgO >= 18 & TiO2 > 1 & TA < 2;
tas_name(ind) = {'meimechite'};

ind = in & S < 45 & MgO >= 18 & TiO2 < 1 & TA < 2;
tas_name(ind & volc == 1) = {'komatiite'};
tas_name(ind & volc == 2) = {'intrusive komatiite'};

ind = in & 45 <= S & S < 52 & MgO >= 18 & TiO2 < 1 & TA < 2;
tas_name(ind & volc == 1) = {'basaltic komatiite'};
tas_name(ind & volc == 2) = {'gabbroic komatiite'};

% classification of ultramafic vs cumulate for plutonic rocks
% Guided by Reverdatto et al. (Russ. Geol. & Geophys., 2008,
% doi:10.1016/j.rgg.2008.01.002)
ind = volc == 2 & 33 <= S & S < 45;
% mafic
tas_name(ind & 0 < TiO2 & TiO2 <= 0.3 & 35 <= MgO) = {'mantle peridotite'};
tas_name(ind & 0 < TiO2 & TiO2 <= 0.3 & 18 <= MgO & MgO < 35) = {'mantle pyroxenite'};
% crustal (this will rename all the above rocks
%tas_name(ind & TiO2 > 0.3 & 22 <= MgO) = {'crustal peridotite'};
%tas_name(ind & TiO2 > 0.3 & 0 <= MgO & MgO < 22) = {'crustal pyroxenite'};

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

% plot uhmg data as negative
uhmglist = {'picrite','alkali picrite','ferropicrite','alkali ferropicrite', ...
    'komatiite','basaltic komatiite', ...
    'meimechite','boninite','sanukitoid','intrusive komatiite', ...
    'mantle peridotite', ...
    'crustal peridotite','mantle pyroxenite','crustal pyroxenite'};
iuhmg = logical(zeros([height(data),1]));
for i = 1:length(uhmglist)
    iuhmg(strcmp(data.rock_type,uhmglist{i})) = 1;
end

% plots tas fields
if ~strcmpi(plottype,'none')
    
    TA(iuhmg) = -TA(iuhmg);
    % 2-d histogram of tas values
    if strcmp(plottype,'hist')
        fig1 = figure;
        hold on;

        eS = [0:0.5:100];
        eTA = [-10:0.2:30];
        
        n = hist2d(S,TA,eS,eTA);
        imagesc(eS,eTA,log10(n));
        colormap([1 1 1; parula]);
        colorbar;
        caxis([-0.1 3]);
        axis xy;
        
        
        % MgO vs TA
        fig2 = figure;
        hold on;
        
        eTA = [0:0.2:10];
        eMg = [8:0.2:50];
        
        n = hist2d(MgO(iuhmg),-TA(iuhmg),eMg,eTA);
        imagesc(eMg,eTA,log10(n));
        colormap([1 1 1; parula]);
        colorbar;
        caxis([-0.1 3]);
        axis xy;
    elseif strcmp(plottype,'scatter');
        fig1 = figure;
        hold on;
        
        if isempty(ptcolor)
            scatter(S,TA,ptsize,'filled');
        else
            scatter(S,TA,ptsize,ptcolor,'filled');
        end
        
        fig2 = figure;
        hold on;
        if isempty(ptcolor) & isempty(ptsize)
            scatter(MgO(iuhmg),-TA(iuhmg),'filled');
        elseif isempty(ptcolor)
            scatter(MgO(iuhmg),-TA(iuhmg),ptsize(iuhmg),'filled');
        elseif isempty(ptsize)
            scatter(MgO(iuhmg),-TA(iuhmg),[],ptcolor(iuhmg),'filled');
        else
            scatter(MgO(iuhmg),-TA(iuhmg),ptsize(iuhmg),ptcolor(iuhmg),'filled');
        end
    else
        warning('Unknown plot type.');
        return;
    end

    figure(fig1);
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
    xlabel('SiO_2 [wt.%]');
    ylabel('K_2O + Na_2O [wt.%]');
    xlim([0 100]);
    ylim([-10 30]);
    set(gca,'Box','on');
    golden;
    
    % MgO vs TA plot
    figure(fig2);
    plot([12 12],[0 30],'w-');
    plot([12 50],[3 3],'w-');
    plot([18 18 50],[0 1 1],'w-');
    plot([18 18 50 50],[1 2 2 1],'w-');
    xlabel('MgO [wt.%]');
    ylabel('K_2O + Na_2O [wt.%]');
    xlim([8 50]);
    ylim([0 30]);
    set(gca,'Box','on');
    golden;
end

return
