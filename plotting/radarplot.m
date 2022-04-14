function radar(R, varargin)
% Create a spider web or radar plot with an axes specified for each column
%
% spider_plot(P, P_labels, axes_interval, axes_precision) creates a spider
% web plot using the points specified in the array P. The column of P
% contains the data points and the rows of P contain the multiple sets of
% data points. Each point must be accompanied by a label specified in the
% cell P_labels. The number of intervals that separate the axes is
% specified by axes_interval. The number of decimal precision points is
% specified by axes_precision.
% 
% P - [vector | matrix]
% P_labels - [cell of strings]
% axes_interval - [integer]
% axes_precision - [integer]
%
% spider_plot(P, P_labels, axes_interval, axes_precision, line_spec) works
% the same as the function above. Additional line properties can be added
% in the same format as the default "plot" function in MATLAB.
%
% line_spec - [character vector]
%
% %%%%%%%%%%%%%%%%%%% Example of a Generic Spider Plot %%%%%%%%%%%%%%%%%%%
% % Clear workspace
% close all;
% clearlbls;
% clc;
% 
% % Point properties
% num_of_points = 6;
% row_of_points = 4;
%
% % Random data
% P = rand(row_of_points, num_of_points);
%
% % Scale points by a factor
% P(:, 2) = P(:, 2) * 2;
% P(:, 3) = P(:, 3) * 3;
% P(:, 4) = P(:, 4) * 4;
% P(:, 5) = P(:, 5) * 5;
%
% % Make random values negative
% P(1:3, 3) = P(1:3, 3) * -1;
% P(:, 5) = P(:, 5) * -1;
% 
% % Create generic labels
% P_labels = cell(num_of_points, 1);
% 
% for ii = 1:num_of_points
%     P_labels{ii} = sprintf('Label %i', ii);
% end
% 
% % Figure properties
% figure('units', 'normalized', 'outerposition', [0 0.05 1 0.95]);
% 
% % Axes properties
% axes_interval = 2;
% axes_precision = 1;
% 
% % Spider plot
% spider_plot(P, P_labels, axes_interval, axes_precision,...
%     'Marker', 'o',...
%     'LineStyle', '-',...
%     'LineWidth', 2,...
%     'MarkerSize', 5);
% 
% % Title properties
% title('Sample Spider Plot',...
%     'Fontweight', 'bold',...
%     'FontSize', 12);
% 
% % Legend properties
% legend('show', 'Location', 'southoutside');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = R.vals;
lbls = R.fields;

% default options
opt = 1;
axlim = [];
axes_interval = 2;
axes_precision = 1;
axes_scale = 'linear';
lloc = 'EastOutside';

% user defined options
while opt + 1 < nargin
    switch lower(varargin{opt})
        case 'axes'
            axlim = varargin{opt+1};
            axlim = [axlim(:,1) fliplr(axlim(:,2:end))];
            opt = opt + 2;
        case 'numticks'
            axes_interval = varargin{opt+1};
            opt = opt + 2;
        case 'precision'
            axes_precision = varargin{opt+1};
            opt = opt + 2;
        case 'legendlocation'
            lloc = varargin{opt+1};
            opt = opt + 2;
        otherwise
            error(['Unknown option, ',varargin{opt}]);
    end
end

%%% Point Properties %%%
% Number of points
[row_count, var_count] = size(P);

%reorder inputs so that variables are arranged clockwise around the plot
P = [P(:,1) fliplr(P(:,2:end))];
lbls = [lbls(:,1) fliplr(lbls(:,2:end))];

%%% Error Check %%%
% Check if axes properties are an integer
if floor(axes_interval) ~= axes_interval || floor(axes_precision) ~= axes_precision
    error('Error: Please enter in an integer for the axes properties.');
end

% Check if axes properties are positive
if axes_interval < 1 || axes_precision < 1
    error('Error: Please enter value greater than one for the axes properties.');
end

% Check if the labels are the same number as the number of points
if length(lbls) ~= var_count
    error('Error: Please make sure the number of labels is the same as the number of points.');
end

% Pre-allocation
max_values = zeros(1, var_count);
min_values = zeros(1, var_count);
axis_increment = zeros(1, var_count);

% Normalized axis increment
normalized_axis_increment = 1/axes_interval;

% Max and min value of each group
if isempty(axlim)
    max_values = nanmax(P,[],1);
    min_values = nanmin(P,[],1);
else
    max_values = axlim(2,:);
    min_values = axlim(1,:);
end
range = max_values - min_values;

ind = range == 0;
range(ind) = 0.2*P(1,ind);
max_values(ind) = P(1,ind) + 0.1*P(1,ind);
min_values(ind) = P(1,ind) - 0.1*P(1,ind);


% Iterate through number of variables
for ii = 1:var_count
    % Group of points
    group_points = [P(:, ii);min_values(ii);max_values(ii)];
    
    % Axis increment
    axis_increment(ii) = range(ii)/axes_interval;
    
    % Normalize points to range from [0, 1]
    P(:, ii) = (P(:, ii) - min(group_points))/range(ii);
    
    % Shift points by one axis increment
    P(:, ii) = P(:, ii) + normalized_axis_increment;   
end

%%% Polar Axes %%%
% Polar increments
polar_increments = 2*pi/var_count;

% Normalized  max limit of axes
axes_limit = 1;

% Shift axes limit by one axis increment
axes_limit = axes_limit + normalized_axis_increment;

[t,~] = cart2pol(0,axes_limit);

% Polar points
radius = [0; axes_limit];
% theta = 0:polar_increments:2*pi;
theta = t:polar_increments:2*pi+t;

% Convert polar to cartesian coordinates
[x_axes, y_axes] = pol2cart(theta, radius);


% Plot polar axes
grey = [1, 1, 1] * 0.5;
h = line(x_axes, y_axes,...
    'LineStyle', '--',...
    'LineWidth', 1,...
    'Color', grey);

% Iterate through all the line handles
for ii = 1:length(h)
    % Remove polar axes from legend
    h(ii).Annotation.LegendInformation.IconDisplayStyle = 'off';
end

%%% Polar Isocurves %%%
% Shifted axes interval
shifted_axes_interval = axes_interval+1;

% Incremental radius
radius = (0:axes_limit/shifted_axes_interval:axes_limit)';

% Convert polar to cartesian coordinates
[x_isocurves, y_isocurves] = pol2cart(theta, radius);

% Plot polar isocurves
hold on;
h = plot(x_isocurves', y_isocurves',...
    'LineStyle', ':',...
    'LineWidth', 1,...
    'Color', grey);

% Iterate through all the plot handles
for ii = 1:length(h)
    % Remove polar isocurves from legend
    h(ii).Annotation.LegendInformation.IconDisplayStyle = 'off';
end

%%% Data Points %%%
% Iterate through all the rows

for ii = 1:row_count
    % Convert polar to cartesian coordinates
    [x_points, y_points] = pol2cart(theta(1:end-1), P(ii, :));
    
    % Make points circular
    x_circular = [x_points, x_points(1)];
    y_circular = [y_points, y_points(1)];
    
    colors = [0, 0.4470, 0.7410;...
        0.8500, 0.3250, 0.0980;...
        0.9290, 0.6940, 0.1250;...
        0.4940, 0.1840, 0.5560;...
        0.4660, 0.6740, 0.1880;...
        0.3010, 0.7450, 0.9330;...
        0.6350, 0.0780, 0.1840];

    % Repeat colors is necessary
    repeat_colors = fix(row_count/size(colors, 1))+1;
    colors = repmat(colors, repeat_colors, 1);
    
    
%     color = rand(1,3);
    % Plot data points
    fill(x_circular, y_circular,colors(ii,:),...
        'FaceColor', colors(ii,:),...
        'EdgeColor', colors(ii,:),...
        'FaceAlpha', 0.1);   
end

%%% Axis Properties %%%
% Figure background
fig = gcf;
fig.Color = 'white';

% Shifted min value
shifted_min_value = min_values - axis_increment;

% Iterate through all the number of points
for hh = 1:var_count
    % Axis label for each row
    row_axis_labels = (shifted_min_value(hh):axis_increment(hh):max_values(hh))';
    

    % Iterate through all the isocurve radius
    for ii = 2:length(radius)
    
        % Display axis text for each isocurve
        text(x_isocurves(ii, hh), y_isocurves(ii, hh), ...
            sprintf(sprintf('%%.%if', axes_precision), ...
            row_axis_labels(ii)),...
            'Units', 'Data',...
            'Color', 'k',...
            'FontSize', 10,...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'middle');
    end
end

% Label points
x_label = x_isocurves(end, :);
y_label = y_isocurves(end, :);

% Shift axis label
shift_pos = 0.07;

% Iterate through each label
for ii = 1:var_count
    % Angle of point in radians
    theta_point = theta(ii);
    
    % Find out which quadrant the point is in
    if theta_point == 0
        quadrant = 0;
    elseif theta_point == pi/2
        quadrant = 1.5;
    elseif theta_point == pi
        quadrant = 2.5;
    elseif theta_point == 3*pi/2
        quadrant = 3.5;
    elseif theta_point == 2*pi
        quadrant = 0;
    elseif theta_point > 0 && theta_point < pi/2
        quadrant = 1;
    elseif theta_point > pi/2 && theta_point < pi
        quadrant = 2;
    elseif theta_point > pi && theta_point < 3*pi/2
        quadrant = 3;
    elseif theta_point > 3*pi/2 && theta_point < 2*pi
        quadrant = 4;
    end
    
    % Adjust text alignment information depending on quadrant
    switch quadrant
        case 0
            horz_align = 'left';
            vert_align = 'middle';
            x_pos = shift_pos;
            y_pos = 0;
        case 1
            horz_align = 'left';
            vert_align = 'bottom';
            x_pos = shift_pos;
            y_pos = shift_pos;
        case 1.5
            horz_align = 'center';
            vert_align = 'bottom';
            x_pos = 0;
            y_pos = shift_pos;
        case 2
            horz_align = 'right';
            vert_align = 'bottom';
            x_pos = -shift_pos;
            y_pos = shift_pos;
        case 2.5
            horz_align = 'right';
            vert_align = 'middle';
            x_pos = -shift_pos;
            y_pos = 0;
        case 3
            horz_align = 'right';
            vert_align = 'top';
            x_pos = -shift_pos;
            y_pos = -shift_pos;
        case 3.5
            horz_align = 'center';
            vert_align = 'top';
            x_pos = 0;
            y_pos = -shift_pos;
        case 4
            horz_align = 'left';
            vert_align = 'top';
            x_pos = shift_pos;
            y_pos = -shift_pos;
    end
    
    % Display text label
    text(x_label(ii)+x_pos, y_label(ii)+y_pos, lbls{ii},...
        'Units', 'Data',...
        'HorizontalAlignment', horz_align,...
        'VerticalAlignment', vert_align);
end

% Axis limits
axis square;
axis([-axes_limit, axes_limit, -axes_limit, axes_limit]);
axis off;

legend(R.types,'Location',lloc);
legend('boxoff');

return