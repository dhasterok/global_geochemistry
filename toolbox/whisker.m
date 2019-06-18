function varargout = whisker(varargin);
% WHISKER - computes quantiles for a whisker plot
%
%   q = whisker(w,x,y) computes the quantiles of a data vector y
%   where w is the width of the boxes centered at x's.
%
%   [qx,qy] = whisker(x,y) where x and y are cells with bins of data
%   vectors, and computes quantiles of both x and y.
%
%   By default a box and whisker plot is produced, but may be explicitly
%   requested whisker(...,'Plot'), or as no plot whisker(...,'NoPlot').
%
%   Options:
%       whisker(..., 'Option',value)
%
%       'Scale' - 'linear' (default) or 'log' (log10)
%       'Color' - color triplet (e.g., [0.5 0.5 0.5])
%       'Quantiles' - five element vector in ascending order
%           [0.05 0.25 0.5 0.75 0.95] (default) or
%           [0.05 0.16 0.5 0.84 0.95] for 1 and 2 standard deviations
%

% 28 Jul 2018 by D. Hasterok

% assign input data
opt = 0;
if nargin < 2
    error('Incorrect number of inputs.');
elseif iscell(varargin{1})
    x = varargin{1};
    y = varargin{2};
    if nargin > 2
        opt = 3;
    end
else
    w = varargin{1};
    x = varargin{2};
    y = varargin{3};
    if nargin > 3
        opt = 4;
    end
end

% defaults
colour = [0 0.447 0.741];
scale = 'linear';
Q = [0.05 0.25 0.5 0.75 0.95];
%Q = [0.05 0.16 0.5 0.84 0.95];
plotflag = 1;
if opt ~= 0
    while opt <= nargin
        %lower(varargin{opt})
        switch lower(varargin{opt})
            case 'color'
                colour = varargin{opt+1};
                opt = opt + 2;
            case 'scale'
                scale = varargin{opt+1};
                opt = opt + 2;
            case 'quantiles'
                Q = varargin{opt+1};
                opt = opt + 2;
            case 'plot'
                plotflag = 1;
                opt = opt + 1;
            case 'noplot'
                plotflag = 0;
                opt = opt + 1;
        end
    end
end


% check for y as cell
if iscell(y)
    np = length(y);
else
    tmp = y;
    [~,np] = size(tmp);
    for i = 1:np
        y{i} = tmp;
    end
end

% set scale
switch scale
case 'log'
    for i = 1:length(y)
        y{i} = log10(y{i});
    end
case 'linear'
    % do nothing
otherwise
    warning('(whisker) Unknown scale option.');
end

for i = 1:np % for number of x's or groups of x
    % compute quantiles
    Y(i,:) = quantile(y{i},Q);

    xbar = [];
    if iscell(x)
        % compute quantiles of x
        X(i,:) = quantile(x{i},Q);
        varargout{1} = X;
        varargout{2} = Y;
    else
        varargout{1} = Y;
    end
    
    % plotting
    if plotflag
        hold on;
        if iscell(x)
            % 95% bars
            xbar = plot([X(i,1),X(i,5)],[Y(i,3),Y(i,3)],'-');
            ybar = plot([X(i,3),X(i,3)],[Y(i,1),Y(i,5)],'-');

            % 25 to 75% box
            fill([X(i,2) X(i,4) X(i,4) X(i,2) X(i,2)], ...
            [Y(i,2) Y(i,2) Y(i,4) Y(i,4) Y(i,2)],colour);

            % 50% point
            point = plot(X(i,3),Y(i,3),'o');
        else
            [nxr,~] = size(x);
            if nxr == 2
                % errorbar
                xbar = plot([x(1,i)-x(2,i) x(1,i)+x(2,i)],[Y(i,3) Y(i,3)],'-');
            end

            % 95% bar
            ybar = plot([x(1,i),x(1,i)],[Y(i,1),Y(i,5)],'-');
            % 25 to 75% box
            if length(w) == 1
                fill([x(1,i)-w/2 x(1,i)+w/2 x(1,i)+w/2 x(1,i)-w/2 x(1,i)-w/2], ...
                    [Y(i,2) Y(i,2) Y(i,4) Y(i,4) Y(i,2)],colour);
            else
                fill([x(1,i)-w(i) x(1,i)+w(i) x(1,i)+w(i) x(1,i)-w(i) x(1,i)-w(i)], ...
                    [Y(i,2) Y(i,2) Y(i,4) Y(i,4) Y(i,2)],colour);
            end

            point = plot(x(1,i),Y(i,3),'o');
        end
    
        if ~isempty(xbar)
            set(xbar,'Color',colour);
        end
        set(ybar,'Color',colour);
        set(point,'MarkerFaceColor','w','Color',colour);
        
        hold off;
    end
end

return
