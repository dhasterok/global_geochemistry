function varargout = barbline(xp,yp,dl,varargin)
% BARBLINE - add barbs to lines
%
%   h = BARBLINE(x,y) will plot barbs (triangles) along the path of a line
%   at regular intervals.
%
%   A number of option value pairs can be added to change properties of the
%   barbs.
%
%   Options:
%       'Color'     supply a normalized color triplet to change the color
%
%       'Scale'     change the size of the arrowhead
%
%       'Reverse'   change side of line that barbs are plotted on
%
% D. Hasterok, U. Adelaide, 2021

% parse inputs
% ------------------------
p = inputParser;

addParameter(p,'Color',[]);
addParameter(p,'LineWidth',0.5,@isnumeric);
addParameter(p,'Scale',1,@isnumeric);
addParameter(p,'Reverse',false);

parse(p,varargin{:});

xp = xp(:);
yp = yp(:);

C = p.Results.Color;
lw = p.Results.LineWidth;
S = p.Results.Scale;

if p.Results.Reverse
    alpha = pi;
else
    alpha = 0;
end
% ------------------------

%plot line
if isempty(C)
    h(1) = plot(xp,yp,'LineWidth',lw);
else
    h(1) = plot(xp,yp,'Color',C,'LineWidth',lw);
end
C = h(1).Color;
hold on;

% if multisegment line, do each individually
ind = find(isnan(xp) | isnan(yp));
if isempty(ind)
    h = [h; plotbarbs(xp,yp,dl,alpha,S,C)];
else
    if ind(1) ~= 1
        xp = [NaN; xp];
        yp = [NaN; yp];
        ind = [1; ind+1];
    end
    
    if ind(end) ~= length(xp)
        xp = [xp; NaN];
        yp = [yp; NaN];
        ind = [ind; length(xp)];
    end

    for i = 2:length(ind)
        h = [h; plotbarbs(xp(ind(i-1)+1:ind(i)-1),yp(ind(i-1)+1:ind(i)-1),dl,alpha,S,C)];
    end
end

if nargout == 1
    varargout{1} = h;
end

return


function h = plotbarbs(xp,yp,dl,alpha,S,C)

xb = [-0.015 0.015 0 -0.015];
yb = [0 0 0.02598 0];

% compute line length
dp = [0; cumsum(sqrt(diff(xp).^2 + diff(yp).^2))];

dmid = dp(end)/2; % middle distance

% determine the number of barbs to add
n = floor(dp(end)/dl);

% limit the number of barbs per line to 20
if n > 20
    n = 20;
    dl = dp(end)/n;
end

% determine center distances of barbs along line
if n == 0
    dc = dmid;
    n = 1;
else
    if mod(n,2) == 0
        dc = linspace(dmid-n/2*dl, dmid+n/2*dl, n);
    else
        dc = linspace(dmid-(n-1)/2*dl, dmid+(n-1)/2*dl, n);
    end
end

% initialize vectors
xc = zeros([n,1]);
yc = xc;
dx = xc;
dy = xc;
theta = xc;

% compute center positions of barbs along line and angles to rotate barbs
for i = 1:n
    [xp,yp]
    dp
    dc(i)
    [xc(i),yc(i),ind] = barbposition(xp,yp,dp, dc(i));
    if ind == 0 || ind > length(xp)
        continue;
    end
    % determine scale and angle of barb
    xl = get(gca,'XLim');
    yl = get(gca,'YLim');

    pbr = pbaspect;
    R = pbr(2)/pbr(1);

    dx(i) = R*diff(xl);
    dy(i) = diff(yl);
    
    dX = (xp(ind+1) - xp(ind))/dx(i);
    dY = (yp(ind+1) - yp(ind))/dy(i);

    theta(i) = atan2(dY,dX) + alpha;
end

% find locations of arrow verticies by scaling and rotating
xh = S*( repmat(cos(theta),1,4).*repmat(xb,n,1) - ...
    repmat(sin(theta),1,4).*repmat(yb,n,1) ) .* repmat(dx,1,4) + repmat(xc,1,4);
yh = S*( repmat(sin(theta),1,4).*repmat(xb,n,1) + ...
    repmat(cos(theta),1,4).*repmat(yb,n,1) ) .* repmat(dy,1,4) + repmat(yc,1,4);

% plot arrows
h = fill(xh',yh',C,'EdgeColor','none');

return


function [xc,yc,i] = barbposition(xp,yp,dp, dc)

i = find(dp <= dc,1,'last');

if dp(i) == dc
    xc = xp(i);
    yc = yp(i);
    
    if i == length(dp)
        i = i - 1;
    end
else
    xc = xp(i) + (xp(i+1) - xp(i))*(dc - dp(i))/(dp(i+1) - dp(i));
    yc = yp(i) + (yp(i+1) - yp(i))*(dc - dp(i))/(dp(i+1) - dp(i));
end

return