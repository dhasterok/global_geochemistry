function tas_name = tas2(data,varargin);
% TAS - Total alkali vs. silica determination.
%
%   data = tas(data) determines the compositionally based rock names for
%   igneous rocks based on the total alkai silica scheme by Middlemost,
%   ESR, 1994.  Plutonic and volcanic rocks are estimated separately.
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

plot_flag = 0;
if nargin == 2
    plot_flag = varargin{1};
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
Mg = data.mgo;

Mg_number = mgnum(data.mgo,data.feo_tot,0);

% fix polygons for type of peridotite based on Mg-number
for i = 1:ntas
    in = inpolygon(S,TA,tasgons{i,1}(:,1),tasgons{i,1}(:,2));

    ind = in & volc == 1;
    tas_name(ind) = {tasgons{i,2}};

    ind = in & volc == 2;
    ind2 = ~strcmp(tasgons{i,3},'peridotite');
    tas_name(ind & ind2) = {tasgons{i,3}};

    % classification of ultramafic vs cumulate for plutonic rocks
    tas_name(ind & ~ind2 & Mg_number >= 80) = {'mantle peridotite'};
    tas_name(ind & ~ind2 & Mg_number < 80) = {'cumulate peridotite'};
end


% classification of high magnesian rocks
for i = 1:nuhm
    in = inpolygon(S,TA,uhmgons{i,1}(:,1),uhmgons{i,1}(:,2)) ...
        & volc == 1 & data.mgo > 18;

    tas_name(in) = {uhmgons{i,2}};
    tas_name(in & i == 3 & data.tio2 > 1) = {uhmgons{i,3}};
end

% assign carbonatite names to high CO2 rocks
data.rock_name = lower(data.rock_name);
ind = find(data.co2 >= 20 | strcmp('carbonatite',data.rock_name));
tas_name(ind) = {'carbonatite'};


% print table of number of rocks in each field
%AND save as a table
num_samples = cell(0);
sample_type = cell(0);

[type,~,ind] = unique(tas_name);
fprintf('No. Samples         Rock Name (calc.)\n');
fprintf('-----------   ----------------------------\n');
for i = 1:length(type);
    num = length(find(ind == i));
    fprintf('%10i    %s\n',num,type{i});
    num_samples{i,1} = num;
    sample_type{i,1} = type{i};
end

tas_table = table(num_samples,sample_type);
tas_table.Properties.VariableNames = {'num_samples','rock_name'};

%Print table of number of rocks in each field - merged for volc/plutonic
fprintf('No. Samples         Rock Name (calc.)\n');
fprintf('-----------   ----------------------------\n');
for i = 1:size(tasgons,1)
    if strcmpi(tasgons{i,2},'')
        ind = strcmpi(tas_table.rock_name(:,1),tasgons{i,3});
        if ~any(ind)
            continue;
        end
        fprintf('%10i    %s\n',tas_table.num_samples{ind},tas_table.rock_name{ind});
    elseif strcmpi(tasgons{i,3},'')
        ind = strcmpi(tas_table.rock_name(:,1),tasgons{i,2});
        if ~any(ind)
            continue;
        end
        fprintf('%10i    %s\n',tas_table.num_samples{ind},tas_table.rock_name{ind});
    else
        ind1 = strcmpi(tas_table.rock_name(:,1),tasgons{i,2});
        ind2 = strcmpi(tas_table.rock_name(:,1),tasgons{i,3});
        if (~any(ind1) && ~any(ind2))
            continue;
        end
        
        if (~any(ind1) && any(ind2))
            fprintf('%10i    %s\n',tas_table.num_samples{ind2},...
            strcat(tas_table.rock_name{ind2},'/',tas_table.rock_name{ind2}));
            continue;
        end
        
        if (any(ind1) && ~any(ind2))
            fprintf('%10i    %s\n',tas_table.num_samples{ind1},...
            strcat(tas_table.rock_name{ind1},'/',tas_table.rock_name{ind1}));
            continue;
        end
        fprintf('%10i    %s\n',tas_table.num_samples{ind1} + tas_table.num_samples{ind2},...
            strcat(tas_table.rock_name{ind1},'/',tas_table.rock_name{ind2}));
    end
end


% plots tas fields
if plot_flag
    figure;
    hold on
    eS = [20:0.5:100];
    eTA = [-6:0.2:30];
    TA(Mg > 18) = -TA(Mg > 18);
    plot(S,TA,'.r')
    xlim([20 100])
    ylim([-6 30])
    
    size(TA)
    size(S)
    bls = regress(TA,[ones(length(S),1) S])
    brob = robustfit(S,TA)
    
    plot(S,bls(1)+bls(2)*S,'r','LineWidth',2);
    plot(S,brob(1)+brob(2)*S,'g','LineWidth',2)
    
    hold off
 
    figure;
    hold on;

    % 2-d histogram of tas values
    eS = [20:0.5:100];
    eTA = [-6:0.2:30];
    TA(Mg > 18) = -TA(Mg > 18);
    n = hist2d(S,TA,eS,eTA);
    imagesc(eS,eTA,log10(n));
    colorbar;
    caxis([-0.1 3]);
    axis xy;

    for i = 1:ntas
        x = tasgons{i,1}(:,1);
        y = tasgons{i,1}(:,2);
        plot(x,y,'w-');
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
        plot(x,y,'w-');
        %fill(x,y,[0.7 0.7 0.7]);
        if isempty(find(volc ~= 2))
            text(mean(x),mean(y),uhmgons{i,2});
        else
            text(mean(x),mean(y),uhmgons{i,3});
        end
    end
    % plots TAS points from geochemistry
    %plot(S,TA,'.')
    xlabel('sio2 [wt.%]');
    ylabel('k2o + na2o [wt.%]');
    xlim([25 100]);
    ylim([-6 20]);
    set(gca,'Box','on');
    golden;

    % make histogram of magnesium number

    figure;
    ind = data.feo_tot > 0 & data.mgo > 0;
    histogram(Mg_number(ind));
    xlabel('Magnesum number');
    ylabel('No. observations');
end

return
