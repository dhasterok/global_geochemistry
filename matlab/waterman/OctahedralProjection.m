function varargout = OctahedralProjection(varargin)
% OCTAHEDRALPROJECTION - computes octahedral coordinates for maps
%
%   Computes coordinates for Waterman and Cahill projections that result in
%   minimal distortion.  The sections can be cut out to produce octohedral
%   representations of the globe.
%
%   [x,y] = OctahedralProjection(lat,lon) converts latitude and longitude
%   into octahedral coordinates for plotting on a map.  The default map is
%   a Waterman Butterfly.
%
%   To convert from octahedral cooridnates to geographic, use the transform
%   option: [lat,lon] = OctahedralProjection(x,y,'Transform','inverse').
%
%   There are several option value pairs:
%
%       'Projection'        choose between Waterman and Cahill
%                           Cahill not yet implemented
%
%       'Style'             'Butterfly' or 'M'
%
%       'Shift'             shifts the map latitude to the east by the 
%                           specified angle in degrees
%
%       'Transform'         either 'forward': (lat,lon) to (x,y) or
%                           'inverse': (x,y) to (lat,lon)
%
%       'Inset'             plot Antarctica as contiguous, true or false
%                           default is false        
%
%       'MaxLat'            if an inset is included, then a maximum
%                           latitude extent for the inset may be included.
%                           default is -60, but -50 will capture the entire
%                           Scotia Plate
%   
% see also OCTAGRATICULE, OCTABORDER, OCTACOAST, OCTARASTER

p = inputParser;
addRequired(p,'Coord1',@isnumeric);
addRequired(p,'Coord2',@isnumeric);
addParameter(p,'Projection','Waterman',@ischar);
addParameter(p,'Style','Butterfly',@ischar);
addParameter(p,'Shift',0,@isnumeric);
addParameter(p,'Transform','forward',@ischar);
addParameter(p,'Inset',false,@islogical);
addParameter(p,'MaxLat',-60,@isnumeric);

parse(p,varargin{:});

% ensure shift is within allowable range
if -180 > p.Results.Shift || p.Results.Shift > 180
    error('Shift must be between -180º and 180º');
end

switch p.Results.Projection
    case 'Waterman'
        config.name = "Waterman";
    case 'Cahill'
        config.name = "Cahill";
    otherwise
        error('Unknown projection.');
end

switch p.Results.Style
    case 'Butterfly'
        config.altitude = 2*sqrt(3);
        config.fullHeight = 4/sqrt(3);
        %config.cutSize = (sqrt(3)-1)/2;
        %config.fullWidth = 4;
        %config.cutWidth = 2;
        %config.cutHeight = sqrt(3);
        %config.hasAspect = true;

        config.octants = [ ...
            0, -1/sqrt(3), -pi/2, -3*pi/4, -pi/2,  pi/2;
            0, -1/sqrt(3), -pi/6,   -pi/4, -pi/2,  pi/2;
            0, -1/sqrt(3),  pi/6,    pi/4, -pi/2,  pi/2;
            0, -1/sqrt(3),  pi/2,  3*pi/4, -pi/2,  pi/2];

        xb = [-0.366025430609034  -0.183012738535297   0.183012692073737   0.366025430609034;
            -0.865267261154332  -0.864989752084107   0.000277509070225   0.865267261154332;
            -2.594521387005865  -2.594521358177850   0.000000028828015   2.594521387005865;
            -3.281088886420939  -3.098076157704126   0.183012728716814   3.281088886420939;
            -3.464101615137754  -2.598076178394547   0.866025436743206   3.464101615137754;
            -3.464101615137754  -2.598076178394547   0.866025436743206   3.464101615137754;
            -3.647114343854569  -3.281088886420941   0.366025457433629   3.647114343854569;
            -4.333681843269643  -3.464101586309739   0.869580256959904   4.333681843269643;
            -6.062935969121177  -3.463824106067529   2.599111863053646   6.062935969121177;
            -6.562177799666475  -3.281088923064018   3.281088876602457   6.562177799666475;
            -6.562177799666475  -3.281088876602457   3.281088923064018   6.562177799666475;
            -6.062935969121177  -2.599111863053646   3.463824106067529   6.062935969121177;
            -4.333681843269643  -0.869580256959904   3.464101586309739   4.333681843269643;
            -3.647114343854569  -0.366025457433629   3.281088886420941   3.647114343854569;
            -3.464101615137754  -0.866025436743206   2.598076178394547   3.464101615137754;
            -3.464101615137754  -0.866025436743206   2.598076178394547   3.464101615137754;
            -3.281088886420939  -0.183012728716814   3.098076157704126   3.281088886420939;
            -2.594521387005865  -0.000000028828015   2.594521358177850   2.594521387005865;
            -0.865267261154332  -0.000277509070225   0.864989752084107   0.865267261154332;
            -0.366025430609034  -0.183012692073737   0.183012738535297   0.366025430609034];
        
        yb = [2.000000026824595   1.683012692073736   1.683012665249142   1.999999973175405;
            2.499241846275631   1.500277493915180   1.001035647639549   1.500758153724369;
            3.497947587918353   0.502052362150060  -0.995895225768293   0.502052412081647;
            3.683012655430659   0.000000000000001  -1.683012655430657   0.316987344569341;
            2.999999962321849  -0.500000019028754  -1.499999980971246   1.000000037678151;
            2.999999962321849  -0.500000019028754  -1.499999980971246   1.000000037678151;
            3.683012655430658  -0.316987344569343  -2.000000000000000   0.316987344569341;
            3.497947587918353  -1.004104774231706  -2.502052362150059   0.502052412081646;
            2.499241846275630  -3.001035647639549  -3.500277493915180   1.500758153724369;
            2.000000026824594  -3.683012665249142  -3.683012692073736   1.999999973175405;
            1.999999973175405  -3.683012692073736  -3.683012665249142   2.000000026824594;
            1.500758153724369  -3.500277493915180  -3.001035647639549   2.499241846275630;
            0.502052412081646  -2.502052362150059  -1.004104774231706   3.497947587918353;
            0.316987344569341  -2.000000000000000  -0.316987344569343   3.683012655430658;
            1.000000037678151  -1.499999980971246  -0.500000019028754   2.999999962321849;
            1.000000037678151  -1.499999980971246  -0.500000019028754   2.999999962321849;
            0.316987344569341  -1.683012655430657   0.000000000000001   3.683012655430659;
            0.502052412081647  -0.995895225768293   0.502052362150060   3.497947587918353;
            1.500758153724369   1.001035647639549   1.500277493915180   2.499241846275631;
            1.999999973175405   1.683012665249142   1.683012692073736   2.000000026824595];
    case 'M'
        config.altitude = 2*sqrt(3);
        config.fullHeight = sqrt(3);
        %config.cutSize = 0;
        %config.fullWidth = 4;
        %config.cutWidth = 0;
        %config.cutHeight = sqrt(3);
        %config.hasAspect = true;

        config.octants = [ ...
            -1, 0, -pi/6, -3*pi/4, -pi/2,  pi/2;
            -1, 0,  pi/6,   -pi/4, -pi/2,  pi/2;
             1, 0, -pi/6,    pi/4, -pi/2,  pi/2;
             1, 0,  pi/6,  3*pi/4, -pi/2,  pi/2];

        xb = [-3.647114353673051  -3.281088923064018   3.281088876602457   3.647114307211491;
            -4.329091367221862  -3.463824106067530   2.599111863053647   3.464379124207979;
            -6.058622973315604  -3.464101586309740   0.869580256959904   3.464101643965769;
            -6.562177772841880  -3.281088886420942   0.366025457433629   3.647114343854567;
            -6.062177793860834  -2.598076178394548   0.866025436743206   4.330127051552428;
            -6.062177793860834  -2.598076178394548   0.866025436743206   4.330127051552428;
            -6.745190501558696  -3.098076157704126   0.183012728716814   3.830127072571383;
            -6.928203201447493  -2.594521358177851   0.000000028828016   4.333681872097658;
            -6.927925721205284  -0.864989752084108   0.000277509070225   6.063213478191400;
            -6.745190538201771  -0.183012738535297   0.183012692073737   6.745190491740211;
            -6.745190491740211  -0.183012692073737   0.183012738535297   6.745190538201771;
            -6.063213478191400  -0.000277509070225   0.864989752084108   6.927925721205284;
            -4.333681872097658  -0.000000028828016   2.594521358177851   6.928203201447493;
            -3.830127072571383  -0.183012728716814   3.098076157704126   6.745190501558696;
            -4.330127051552428  -0.866025436743206   2.598076178394548   6.062177793860834;
            -4.330127051552428  -0.866025436743206   2.598076178394548   6.062177793860834;
            -3.647114343854567  -0.366025457433629   3.281088886420942   6.562177772841880;
            -3.464101643965769  -0.869580256959904   3.464101586309740   6.058622973315604;
            -3.464379124207979  -2.599111863053647   3.463824106067530   4.329091367221862;
            -3.647114307211491  -3.281088876602457   3.281088923064018   3.647114353673051];


        yb = [2.683012692073735   2.683012665249141   2.683012692073735   2.683012665249141;
            2.500277493915179   2.001035647639549   2.500277493915179   2.001035647639549;
            1.502052362150060   0.004104774231707   1.502052362150059   0.004104774231707;
            1.000000000000000  -0.683012655430658   1.000000000000000  -0.683012655430658;
            0.499999981160924  -0.499999980971247   0.499999980971246  -0.499999981160925;
            0.499999981160924  -0.499999980971247   0.499999980971246  -0.499999981160925;
            0.683012655430657  -1.000000000000001   0.683012655430657  -1.000000000000001;
            -0.004104774231707  -1.502052362150061  -0.004104774231708  -1.502052362150061;
            -2.001035647639549  -2.500277493915180  -2.001035647639549  -2.500277493915180;
            -2.683012665249142  -2.683012692073736  -2.683012665249142  -2.683012692073736;
            -2.683012692073736  -2.683012665249142  -2.683012692073736  -2.683012665249142;
            -2.500277493915180  -2.001035647639549  -2.500277493915180  -2.001035647639549;
            -1.502052362150061  -0.004104774231708  -1.502052362150061  -0.004104774231707;
            -1.000000000000001   0.683012655430657  -1.000000000000001   0.683012655430657;
            -0.499999981160925   0.499999980971246  -0.499999980971247   0.499999981160924;
            -0.499999981160925   0.499999980971246  -0.499999980971247   0.499999981160924;
            -0.683012655430658   1.000000000000000  -0.683012655430658   1.000000000000000;
            0.004104774231707   1.502052362150059   0.004104774231707   1.502052362150060;
            2.001035647639549   2.500277493915179   2.001035647639549   2.500277493915179;
            2.683012665249141   2.683012692073735   2.683012665249141   2.683012692073735];
    case 'ArcticFlower'
        config.altitude = 2*sqrt(3);
        config.fullHeight = 4/sqrt(3);
        %config.cutSize = (sqrt(3)-1)/2;
        %config.fullWidth = 4;
        %config.cutWidth = 2;
        %config.cutHeight = sqrt(3);
        %config.hasAspect = true;

        % octant coordinates
        % xp shift, yp shift, rotation angle, shift
        config.octants = [ ...
            0, 0, -3*pi/4, -3*pi/4, -pi/2,  pi/2;
            0, 0, -pi/4,   -pi/4, -pi/2,  pi/2;
            0, 0,  pi/4,    pi/4, -pi/2,  pi/2;
            0, 0,  3*pi/4,  3*pi/4, -pi/2,  pi/2];
    case 'AntarcticFlower'
    otherwise
        error('Unknown style.');
end

switch p.Results.Transform
    case 'forward'
        lat = p.Results.Coord1;
        lon = p.Results.Coord2;

        [nr,nc] = size(lon);
        lon = lon(:);
        lat = lat(:);
        
        % copy lat and lon for Antarctic region
        if p.Results.Inset
            ind = lat <= p.Results.MaxLat + 1e-4;
            if any(ind)
                % if the file contains segments separated by NaN's then add
                % NaN's to list of copied points
                ind2 = find(ind) + 1;
                if ind2(end) > length(lat)
                    ind2 = ind2(1:end-1);
                end
                
                ind(ind2(isnan(lat(ind2)))) = true;
    
                alat = lat(ind);
                alon = lon(ind);
    
                alat = alat*pi/180;
                alon = alon*pi/180;
    
                % convert to inset coordinates
                xa = nan(size(alon));
                ya = xa;
    
                ind = ~isnan(alat) & ~isnan(alon);
                [xa(ind),ya(ind)] = southpoleForward(alat(ind),alon(ind), config);
    
                % shift Antarctica to below map
                switch p.Results.Style
                    case 'Butterfly'
                        ya = ya - 4;
                    case 'M'
                        ya = ya - 3;
                        xa = xa - 2*sqrt(3);
                    case 'ArcticFlower'
                    case 'AntarcticFlower'
                end
            else
                xa = [];
                ya = [];
            end
        end

        % shift coordinates if central meridian is not zero
        if p.Results.Shift ~= 0
            lon = lon + p.Results.Shift;
            ind = lon > 180;
        
            lon(ind) = lon(ind) - 360;
        
            ind = lon < -180;
            lon(ind) = 360 + lon(ind);
        end
        
        % convert to radians
        lon = lon*pi/180;
        lat = lat*pi/180;
        
        x = nan(size(lon));
        y = x;

        % convert to octagonal coordinates
        ind = ~isnan(lat) & ~isnan(lon);
        [x(ind),y(ind)] = forward(lat(ind),lon(ind), config);

        if p.Results.Inset
            
        end
        
        x = reshape(x,nr,nc);
        y = reshape(y,nr,nc);
        
        if nargout == 2
            varargout = {x,y};
        elseif nargout == 4
            if exist('xa','var')
                varargout = {x,y,xa,ya};
            else
                varargout = {x,y,[],[]};
            end
        end
    case 'inverse'
        x = p.Results.Coord1;
        y = p.Results.Coord2;
        
        [nr,nc] = size(x);
        x = x(:);
        y = y(:);

        tol = 1e-4;
        quadrant = nan(size(x));
        for q = 1:4
            for n = 1:2
                for m = 1:4
                    if n == 1 && m == 1
                        [in,on] = inpolygon(x,y,xb(:,q),yb(:,q));
                        
                        quadrant(in | on) = q;
                    end
                    k = isnan(quadrant);
                    [in,on] = inpolygon(x+tol*(-1)^n,y+tol*(-1)^m,xb(:,q),yb(:,q));
                    
                    quadrant(k & (in | on)) = q;
                end
            end
        end
        %quadrant
        %find(isnan(quadrant))
        
        lon = nan(size(x));
        lat = lon;
        ind = ~isnan(quadrant);
        [lat(ind),lon(ind)] = inverse(x(ind), y(ind), quadrant(ind), config);
        
        % convert to degrees
        lon = lon*180/pi;
        lat = lat*180/pi;
        
        % shift coordinates if central meridian is not zero
        if p.Results.Shift ~= 0
            lon = lon - p.Results.Shift;

            ind = lon > 180;
            lon(ind) = lon(ind) - 360;
        
            ind = lon < -180;
            lon(ind) = 360 + lon(ind);
        end

        lon = reshape(lon,nr,nc);
        lat = reshape(lat,nr,nc);

        if p.Results.Inset
            maxlat = p.Results.MaxLat;

            blon2 = [linspace(-179.9999,-90.0001,90) -90.0001 -90.0001 -89.9999 -89.9999 ...
                linspace(-89.9999,-0.0001,90) -0.0001 -0.0001 0.0001 0.0001 ...
                linspace(0.0001,89.9999,90) 89.9999 89.9999 90.0001 90.0001 ...
                linspace(90.0001,179.9999,90) 179.9999 179.9999 -179.9999 -179.9999];
            blat2 = [repmat(maxlat,1,90) maxlat -71.38865281 -71.38865 maxlat ...
                repmat(maxlat,1,90) maxlat -71.38865281 -71.38865281 maxlat ...
                repmat(maxlat,1,90) maxlat -71.38865281 -71.38865281 maxlat ...
                repmat(maxlat,1,90) maxlat -71.38865281 -71.38865281 maxlat];
        
            [~,~,xv,yv] = OctahedralProjection(blat2,blon2, ...
                'Projection',p.Results.Projection,'Style',p.Results.Style,'Inset',true,'MaxLat',maxlat);
            
            %plot(x,y,'o');
            %hold on;
            %plot(xv,yv);
            [in,on] = inpolygon(x,y,xv,yv);

            ind = in | on;

            xa = x;
            ya = y;

            switch p.Results.Style
                case 'Butterfly'
                    ya = ya + 4;
                case 'M'
                    ya = ya + 3;
                    xa = xa + 2*sqrt(3);
            end

            % quadrants numbered from clockwise from lower left
            quadrant = 2*ones(size(xa));
            ind = xa < 0 & ya < 0;
            quadrant(ind) = 1;
            ind = xa > 0 & ya < 0;
            quadrant(ind) = 4;
            ind = xa > 0 & ya > 0;
            quadrant(ind) = 3;
            
            theta = (pi/2*quadrant - pi/4);
            
            [xP,yP,x0,y0] = forward(-90*pi/180,-180*pi/180,config);
            
            xMj = x0 + sin(theta).*xa + cos(theta).*ya;
            yMj = y0 - cos(theta).*xa + sin(theta).*ya;

            % finish forward calculation so that inverse function can be used
            % vertex coordinates
            xP = config.octants(quadrant,1)*config.altitude;
            yP = config.octants(quadrant,2)*config.altitude;
            
            theta = config.octants(quadrant,3);

            % rotate points onto the correct octant
            xtmp = xP + sin(theta).*xMj + cos(theta).*yMj;
            ytmp = yP - cos(theta).*xMj + sin(theta).*yMj + config.fullHeight*config.altitude/2;

            % compute inverse with Antarctica inset
            [alat,alon] = inverse(xtmp,ytmp, quadrant, config);

            alat = 180/pi*alat';
            alon = 180/pi*alon';
        end
        
        if nargout == 2
            varargout = {lat,lon};
        else
            if exist('alat','var')
                varargout = {lat,lon,alat,alon};
            else
                varargout = {lat,lon,[],[]};
            end 
            
        end

    otherwise
        error('Unknown transform.');
end

return


% FORWARD - converts form (lat,lon) to (x,y) coordinates
function [x,y,varargout] = forward(lat, lon, config)

x = nan(size(lat));
y = x;

quadrant = floor(2*(lon+pi)/pi) + 1;
ind = quadrant > 4;
quadrant(ind) = quadrant(ind) - 4;
quadrant(isnan(quadrant) | quadrant < 1) = 1;

% relative longitude
%ind = sub2ind(size(config.octants),quadrant,repmat(4,size(quadrant)));
lonr = mod(lon - config.octants(quadrant,4) + pi, 2*pi) - pi;

% vertex coordinates
xP = config.octants(quadrant,1)*config.altitude;
yP = config.octants(quadrant,2)*config.altitude;

% rotation angle
theta = config.octants(quadrant,3);

% project the coordinates to a generic octant
%[lat*180/pi lon*180/pi abs(lat)*180/pi, abs(lonr)*180/pi]
[x,y] = faceForward(abs(lat), abs(lonr));

%relative octant coordinates (Mj stands for "Mary Jo Graca")
xMj = x;
yMj = y;
             
%reflect the southern hemisphere over the equator
ind = (lat < 0);
xMj(ind) = 2*config.altitude - xMj(ind);

ind = (lonr < 0);
yMj(ind) = -yMj(ind);

if nargout == 4
    varargout = {xMj yMj};
end

% rotate points onto the correct octant
x = xP + sin(theta).*xMj + cos(theta).*yMj;
y = yP - cos(theta).*xMj + sin(theta).*yMj + config.fullHeight*config.altitude/2;

return


% FACEFORWARD - converts a generic face (lat,lon) to (x,y)
function [x,y] = faceForward(lat, lon)

x = nan(size(lat));
y = x;

[xELD,yELD] = jointHeights(lon);

% now that we've established the meridian, lets place our point on it
l = sqrt(diff(xELD,1,2).^2 + diff(yELD,1,2).^2);

desirLength = repmat(linInterp(lat, pi/2, 0, 0, totalLength(xELD, yELD)),1,3) - ...
    [zeros(size(lat)) cumsum(l(:,1:end-1),2)];

% if it fits on this segment
ind = l >= desirLength;

i = 4 - sum(ind,2);
i(i == 4) = 3;

ind = sub2ind(size(desirLength),[1:size(desirLength,1)]',i);
ind2 = sub2ind(size(xELD),[1:size(xELD,1)]',i);
ind3 = sub2ind(size(xELD),[1:size(xELD,1)]',i+1);

% interpolate and return
x = linInterp(desirLength(ind), 0, l(ind), xELD(ind2), xELD(ind3));
y = linInterp(desirLength(ind), 0, l(ind), yELD(ind2), yELD(ind3));

return


% SOUTHPOLEFORWARD - converts (lat,lon) to (x,y) for Antarctic inset
function [x,y] = southpoleForward(lat,lon,config)

% with Antarctica Inset
[~,~,xMj,yMj] = forward(lat,lon, config);

% rotate Antarctic points to correct orientation
if exist('xMj','var')
    xMj = xMj;
    yMj = -yMj;

    quadrant = floor(2*(lon+pi)/pi) + 1;
    theta = (pi/2*quadrant - pi/4);

    % reset south pole to origin
    [~,~,x0,y0] = forward(-90*pi/180,-180*pi/180,config);
    %plot(x0-x0,y0-y0,'b+')
    x = sin(theta).*(xMj-x0) + cos(theta).*(yMj-y0);
    y = -(-cos(theta).*(xMj-x0) + sin(theta).*(yMj-y0));
    %hold on;
    %plot(xa,ya);
end

return


% INVERSE - converts (x,y) to (lat,lon) coordinates
function [lat,lon] = inverse(x, y, quadrant, config)

y = y - config.fullHeight*config.altitude/2;

lat = nan(size(x));
lon = lat;

% try each octant

% vertex coordinates
xV = config.octants(quadrant,1)*config.altitude;
yV = config.octants(quadrant,2)*config.altitude;

% rotation angle
theta = config.octants(quadrant,3);

% convert to generic coordinates
xMj = sin(theta).*(x-xV) - cos(theta).*(y-yV);
yMj = cos(theta).*(x-xV) + sin(theta).*(y-yV);

%if (sqrt(3)*abs(yMj) > min(xMj, 2*config.altitude-xMj) + 1e-12) % if the angle is wrong,
    %continue % check the next one
%end

% invert generic coordinates
%[min([xMj, 2*config.altitude-xMj],[],2), abs(yMj)]
[lat,lon] = faceInverse(min([xMj, 2*config.altitude-xMj],[],2), abs(yMj));
% now find the correct octant and shift accordingly
% if you got not found, keep looking
%if isnan(lat) && isnan(lon)
    %continue
%end
 
lat = lat.*sign(config.altitude-xMj); % undo the reflections
lon = lon.*sign(yMj);

% if the resulting coordinates are wrong
%if (lat < config.octants(i,5) - 1e-6 || lat > config.octants(i,6) + 1e-6) 
    %continue % move on
%end

lon = mod(lon + config.octants(quadrant,4) + pi, 2*pi) - pi;

return


% FACEINVERSE - converts a generic face (x,y) to (lat,lon)
function [lat,lon] = faceInverse(x,y)

nx = length(x);

lat = nan(size(x));
lon = lat;

xELD = [(sqrt(3)-1)/2, sqrt(3)/2, 3*sqrt(3)/2, NaN];
dYdL = [0, 2/pi, 6/pi, (4 + sqrt(8))/pi];
sin15 = (sqrt(3) - 1)/sqrt(8);
cos15 = (sqrt(3) + 1)/sqrt(8);
lonDiag = pi/(4+sqrt(8));

% generic quadrant coordinates
% v = [(sqrt(3) - 1)/2,  0;
%    sqrt(3)/2,         1/2;
%    3*sqrt(3)/2,       3/2;
%    (1 + 7*sqrt(3))/4, (5 + sqrt(3))/4;
%    2*sqrt(3),          1];
% v = [v; flipud(v(:,1)) -flipud(v(:,2))];

% this describes the footprint of the octant
inside = (y <= x-(sqrt(3)-1)/2 & y <= x/sqrt(3) & y <= (2-sqrt(3))*(x+3) & ...
    y <= (7+4*sqrt(3))-(2+sqrt(3))*x  & x <= 2*sqrt(3));

if all(~inside)
    return
end

% figure out in which segment it is
ind = repmat(xELD,nx,1) < repmat(x,1,4);
i = sum(ind,2) + 1;
i(i > 4) = 4;
i(i == 1) = 2;

% well-behaved segments?
ind = i <= 3;
lon(ind) = y(ind)./linInterp(x(ind), xELD(i(ind)-1)', xELD(i(ind))', dYdL(i(ind)-1)', dYdL(i(ind))');

% the well-behaved part of the last segment?
ind = i > 3;
lon(ind) = y(ind)./linInterp(x(ind), xELD(i(ind)-1)', 2*sqrt(3), dYdL(i(ind)-1)', dYdL(i(ind))');
    
% the diagonal part of the last segment?
ind = ind & lon > lonDiag;
% surprisingly, the equation becomes quadratic here
a = dYdL(3)*dYdL(4)*sin15; 
b = (dYdL(4)*cos15-dYdL(3))*(1.5*sqrt(3) - x(ind)) - dYdL(4)*sin15*y(ind) - dYdL(3)*(sqrt(3)/2+dYdL(4)*lonDiag*sin15);
c = (sqrt(3)/2+dYdL(4)*lonDiag*sin15)*y(ind) + (1-dYdL(4)*lonDiag*cos15)*(1.5*sqrt(3)-x(ind));
lon(ind) = (-b - sqrt(b.*b - 4*a*c))/(2*a);


[xELD,yELD] = jointHeights(lon);
ind1 = sub2ind(size(xELD),[1:nx]',i);
ind2 = sub2ind(size(xELD),[1:nx]',i-1);
pointLength = sqrt((x - xELD(ind2)).^2 + (y - yELD(ind2)).^2);

% add in the lengths of all previous segments
sxELD = cumsum([zeros(size(x)) sqrt(diff(xELD,1,2).^2 + diff(yELD,1,2).^2)],2);

pointLength = pointLength + sxELD(ind2);
totlen = totalLength(xELD, yELD);

lat = linInterp(pointLength, 0, totlen, pi/2, 0);

%lat(~inside) = NaN;
%lon(~inside) = NaN;

return


% SOUTHPOLEINVERSE - converts (x,y) to (lat,lon) for Antarctic inset
function [lat,lon] = southpoleInverse(x,y,config)

return


% JOINTHEIGHTS - computes coordinates
function [xELD,yELD] = jointHeights(lon)

nl = length(lon);
LON = repmat(lon,1,4);

xELD = repmat([(sqrt(3)-1)/2, sqrt(3)/2, 3*sqrt(3)/2, NaN],nl,1);
dYdL = [0, 2/pi, 6/pi, (4 + sqrt(8))/pi];
lonDiag = pi/(4 + sqrt(8));

% sine and cosine of 15º
sin15 = (sqrt(3) - 1)/sqrt(8);
cos15 = (sqrt(3) + 1)/sqrt(8);

% The height of the intersections of this meridian with the Equal Line
% Delineations (ELD)
yELD = repmat(dYdL,nl,1).*LON;

% if it hits the equatorial ELD on the vertical (hexagon) part
ind = abs(lon) <= lonDiag;
xELD(ind,4) = 2*sqrt(3);

% if it hits it on the diagonal (square) part
ind = ~ind;
xELD(ind,4) = -(abs(lon(ind)) - lonDiag)*dYdL(4)*sin15 + 2*sqrt(3);
yELD(ind,4) = (abs(lon(ind)) - lonDiag)*dYdL(4)*cos15 + 1;

return


function len = totalLength(xs, ys)

% compute meridian length
len = sum(sqrt(diff(xs,1,2).^2 + diff(ys,1,2).^2),2);

return
