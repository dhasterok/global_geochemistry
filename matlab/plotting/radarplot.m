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
%   The plot can be modified using the following name value pairs:
%
%       'AxisLim'           Set the axes limits along each dimension, a
%                           matrix 2 by N is required.
%
%       'NumTicks'          Define the number of tick marks to be placed on
%                           the axes.
%
%       'Precision'         
%
%       'Scale'             Change between 'linear' (default) and 'log' 
%                           scale
%
%       'Location'          Change the location of the legend using the
%                           same location options as legend().  The default
%                           is 'EastOutside'
%
%       'Axes'              axes graphics object
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

% user defined options
p = inputParser;

addParameter(p,'AxisLim',[],@isnumeric);
addParameter(p,'NumTicks',2,@isnumeric);
addParameter(p,'Precision',1,@isnumeric);
addParameter(p,'Scale','linear',@ischar);
addParameter(p,'Location','northeastoutside',@ischar);
addParameter(p,'Orientation','vertical',@ischar);
addParameter(p,'Axes',[],@isgraphics);     % axes
addParameter(p,'Colors',[]);

parse(p,varargin{:});
axlim = p.Results.AxisLim;
axes_interval = p.Results.NumTicks;
axes_precision = p.Results.Precision;
axes_scale = p.Results.Scale;
legend_loc = p.Results.Location;
ax = p.Results.Axes;
if isempty(ax);
    ax = gca;
end
C = p.Results.Colors;
orientation = p.Results.Orientation;

%%% Point Properties %%%
% Number of points
[group_count, field_count, q_count] = size(R.vals);

if isempty(C)
    % set colors
    C = [0, 0.4470, 0.7410;...
        0.8500, 0.3250, 0.0980;...
        0.9290, 0.6940, 0.1250;...
        0.4940, 0.1840, 0.5560;...
        0.4660, 0.6740, 0.1880;...
        0.3010, 0.7450, 0.9330;...
        0.6350, 0.0780, 0.1840];
    
    % Repeat colors is necessary
    repeat_colors = fix(group_count/size(C, 1))+1;
    C = repmat(C, repeat_colors, 1);
else
    % for some reason this will throw an error with some valid color maps
    %if ~all(validatecolor(C,'multiple'))
    %    error('Color input is invalid.');
    %end
end

%reorder inputs so that variables are arranged clockwise around the plot
%vals = [vals(:,1) fliplr(vals(:,2:end))];
%lbls = [lbls(:,1) fliplr(lbls(:,2:end))];

%%% Error Check %%%
% Check if axes properties are an integer
if floor(axes_interval) ~= axes_interval || floor(axes_precision) ~= axes_precision
    error('Please enter in an integer for the axes properties.');
end

% Check if axes properties are positive
if axes_interval < 1 || axes_precision < 1
    error('Please enter value greater than one for the axes properties.');
end

% Check if the labels are the same number as the number of points
if length(R.fields) ~= field_count
    error('Please make sure the number of labels is the same as the number of points.');
end

% Pre-allocation
%axis_increment = zeros(1, field_count);

% Determine min and max value of each field (for establishing axes limits)
if isempty(axlim)
    fieldmin = min(R.vals(:,:,1),[],1,'omitnan');
    fieldmax = max(R.vals(:,:,end),[],1,'omitnan');
else
    fieldmin = axlim(1,:);
    fieldmax = axlim(2,:);
end

ind = fieldmax - fieldmin == 0;
fieldmin(ind) = R.vals(1,ind,1) - 0.1*R.vals(1,ind,1);
fieldmax(ind) = R.vals(1,ind,end) + 0.1*R.vals(1,ind,end);
range = fieldmax - fieldmin;


% Normalized axis increment
normalized_axis_increment = 1/axes_interval;
axis_increment = range./axes_interval;

% Normalize data within each field [min,max] -> [0,1]+one tick
% the plus one tick to the normalized values moves the minimum away from the
% center of the plot.
unitvec = ones([group_count,1]);
for k = 1:q_count
    R.vals(:,:,k) = (R.vals(:,:,k) - unitvec*fieldmin)./(unitvec*range) + normalized_axis_increment;
end

%%% Polar Axes %%%
% Polar increments
polar_increments = 2*pi/field_count;

% Normalized  max limit of axes
axes_limit = 1;

% Shift axes limit by one axis increment
axes_limit = axes_limit + normalized_axis_increment;

[t,~] = cart2pol(0,axes_limit);

% Polar points
radius = [0; axes_limit];
% theta = 0:polar_increments:2*pi;
theta = t:polar_increments:2*pi+t;
%theta = [fliplr(theta)];

% Convert polar to cartesian coordinates
[x_axes, y_axes] = pol2cart(theta, radius);


% Plot polar axes
grey = [1, 1, 1] * 0.5;
h = line(ax,x_axes, y_axes,...
    'LineStyle', '--',...
    'LineWidth', 1,...
    'Color', grey);

% Iterate through all the line handles
for i = 1:length(h)
    % Remove polar axes from legend
    h(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
end

%%% Polar Isocurves %%%
% Shifted axes interval
shifted_axes_interval = axes_interval+1;

% Incremental radius
radius = (0:axes_limit/shifted_axes_interval:axes_limit)';

% Convert polar to cartesian coordinates
[x_isocurves, y_isocurves] = pol2cart(theta, radius);

% Plot polar isocurves
hold(ax,'on');
h = plot(ax,x_isocurves', y_isocurves',...
    'LineStyle', ':',...
    'LineWidth', 1,...
    'Color', grey);

% Iterate through all the plot handles
for i = 1:length(h)
    % Remove polar isocurves from legend
    h(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
end

%%% Data Points %%%
% Iterate through all the rows

if q_count == 1
    for i = 1:group_count
        mid(i) = radarline(ax, theta,R.vals(i,:),C(i,:));  
    end
else
    switch q_count
        case 2
            for i = 1:group_count
                mid(i) = radarpoly(ax, theta,[R.vals(i,:,1); R.vals(i,:,2)],C(i,:));  
            end
        case 3
            for i = 1:group_count
                radarpoly(ax, theta,[R.vals(i,:,1); R.vals(i,:,3)],C(i,:));  
                mid(i) = radarline(ax, theta,R.vals(i,:,2),C(i,:));  
            end
        case 5
            for i = 1:group_count
                low(i) = radarline(ax, theta,R.vals(i,:,1),C(i,:));
                high(i) = radarline(ax, theta,R.vals(i,:,5),C(i,:));
                set(low(i),'LineWidth',0.25);
                set(high(i),'LineWidth',0.25);

                radarpoly(ax, theta,[R.vals(i,:,1); R.vals(i,:,2)],C(i,:));  
                
                mid(i) = radarline(ax, theta,R.vals(i,:,3),C(i,:));
            end
    end
end

%%% Axis Properties %%%
% Figure background
ax.Color = 'white';
ax.XColor = 'none';
ax.YColor = 'none';

% Shifted min value
shifted_min_value = fieldmin - axis_increment;

% Iterate through all the number of points
for j = 1:field_count
    % Axis label for each row
    row_axis_labels = (shifted_min_value(j):axis_increment(j):fieldmax(j))';
    
    if isempty(row_axis_labels)
        continue;
    end

    % Iterate through all the isocurve radius
    for i = 2:length(radius)
    
        % Display axis text for each isocurve
        text(ax,x_isocurves(i, j), y_isocurves(i, j), ...
            sprintf(sprintf('%%.%if', axes_precision), ...
            row_axis_labels(i)),...
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
for i = 1:field_count
    % Angle of point in radians
    theta_point = theta(i);
    
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
    text(ax,x_label(i)+x_pos, y_label(i)+y_pos, R.fields{i},...
        'Units', 'Data',...
        'HorizontalAlignment', horz_align,...
        'VerticalAlignment', vert_align);
end

% Axis limits
axis(ax,'square');
axis(ax,1.1*[-axes_limit, axes_limit, -axes_limit, axes_limit]);
axis(ax,'off');
ax.XTick = [];
ax.YTick = [];


if isnumeric(R.types)
    for i = 1:height(R.types)
        types{i} = num2str(R.types(i));
    end
else
    types = R.types;
end

legend(ax, mid,types,'Location',legend_loc,'Orientation',orientation);
legend(ax,'boxoff');

return


function h = radarline(ax,theta,vals,C)

% Convert polar to cartesian coordinates
[x_points, y_points] = pol2cart(theta(1:end-1), vals);

% Make points circular
x_circular = [x_points, x_points(1)];
y_circular = [y_points, y_points(1)];

h = plot(ax,x_circular, y_circular,'Color',C,'LineWidth',1);

return


function h = radarpoly(ax,theta,vals,C)


% Convert polar to cartesian coordinates
[x_points, y_points] = pol2cart(repmat(theta(1:end-1),height(vals),1), vals);

% Make points circular
x_circular = [x_points, x_points(:,1)];
y_circular = [y_points, y_points(:,1)];

if height(vals) == 1
    h = fill(ax, x_circular,y_circular,C, ...
        'FaceColor', C,...
        'EdgeColor', C,...
        'FaceAlpha', 0.1);
else
    outer = polyshape([x_circular(2,:); y_circular(2,:)]');
    inner = polyshape([x_circular(1,:); y_circular(1,:)]');
    
    poly = subtract(outer,inner);

    h = plot(ax, poly, ...
        'FaceColor', C, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.1);
end

return