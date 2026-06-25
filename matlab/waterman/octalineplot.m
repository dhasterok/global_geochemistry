function varargout = octalineplot(varargin)

% parse inputs
p = inputParser;
addRequired(p,'Latitude',@isnumeric);
addRequired(p,'Longitude',@isnumeric);
addParameter(p,'Projection','Waterman',@ischar);
addParameter(p,'Style','Butterfly',@ischar);
addParameter(p,'Shift',0,@isnumeric);
addParameter(p,'Transform','forward',@ischar);
addParameter(p,'Inset',false,@islogical);
addParameter(p,'MaxLat',-60,@isnumeric);
addParameter(p,'Save','',@ischar);

parse(p,varargin{:});

[x,y,xa,ya] = cutpolyconvert(p.Results.Latitude,p.Results.Longitude,p);

ax = gca;
hold(ax,'on');
prepmap;
h(1) = plot(ax,x,y,'k-');

if p.Results.Inset
    h = [h, plot(ax,xa,ya,'k-')];
end

if nargout == 1
    varargout{1} = h;
end

if ~strcmp(p.Results.Save,'')
    if p.Results.Inset
        x = [x; NaN; xa];
        y = [y; NaN; ya];
    end
    s = geoshape(y,x,'Geometry','line');

    filename = [p.Results.Save,'_',p.Results.Projection,'_',p.Results.Style,'_Shift_',num2str(p.Results.Shift)];
    shapewrite(s,filename);
end

return