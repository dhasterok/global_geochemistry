function varargout = adjbyquantile2(xv,yv,varargin)
% ADJBYQUANTILE - normalizes to a reference position along a quantile
% contour
%
%   y_at_ref = adjbyquantile(xv,yv,xref,x,y) adjusts (normalizes for a
%   set) the position of y from x to xref along the quantile where y sits
%   with respect to the global set of observations (xv,yv)
%
%   Inputs:
%       xv - position of all x's
%       yv - position of all y's
%       xref - reference position to normalize to
%       x - observation point x
%       y - observation point y
%
%   Most of the output are used for plotting in tracenorm. If no xref is
%   provided, only outputs with a '*' are returned.
%   Outputs:
%       y_at_ref        normalize location of y w.r.t. xref
%       Qy              estimated quantile level for each y
%       ymodel          model of scale parameters for normalized y, if no
%                       subset is given, ymodel uses the full dataset (yv)
%           .mu,            estimated mean in log-normal space
%           .sigma,         estimated standard deviation in log-normal
%                           space
%           .sigma_mu,      uncertainty in mu
%           .rms,           misfit between empirical and modeled CDF
%       yref            values y projected to xref
%       Yref            values of yv projected to xref 
%       Q               quantile levels used for fitting
%       xq *            x bin midpoints for plotting yq
%       yq *            yv bin quantiles at levels in Q
%       Yq *            smoothed yv quantiles at levels in Q 
%       model *         empirical model of scale parameters produced from
%                       yv within each x bin
%           .x              x bin midpoints
%           .mu             mu for each x bin
%           .sigma          sigma for each x bin
%           .sigma_mu       uncertainty of mu for each x bin
%           .rms            rms for each empirical and modeled CDF for x bins
%       smooth_model *  parameters for producing smooth quantiles
%           .mm,            quadratic coefficients for mu scale parameter
%           .mm_alpha95,    confidence interval on .mm coefficients
%           .mm_rms,    	misfit for model .mm
%           .ms,            linear coefficients for sigma scale parameter
%           .ms_alpha95,    confidence interval on .ms coefficients
%           .ms_rms,    	misfit for model .ms
%
% See spidernorm, tracenorm, gausscensor

% parse inputs
xref = [];
if nargin > 2
    xref = varargin{1};
    if nargin == 5
        x = varargin{2};
        y = varargin{3};
    else
        x = xv;
        y = yv;
    end
end


% Step 1 - create x bins
minx = min(xv);
maxx = max(xv);

xedges = linspace(minx,maxx,21);
xmid = midpt(xedges);


% Step 2 & 3 - determine empirical CDF for each bin and estimate scale
% parameters
Q = [0.0001 0.5 2.5 10:10:90 97.5 99.5 99.999]'/100;
F = sqrt(2)*erfinv(2*Q - 1);

% compute quantiles of yv w.r.t. binned ranges of xv
xq = repmat(xmid,length(Q),1);
yq = zeros([length(Q) length(xedges)-1]);
for i = 1:length(xedges)-1
    ind = xedges(i) <= xv & xv <= xedges(i+1);
    
    % compute quantiles and log-normal scale factors (model)
    [tmp,yqtmp] = gausscensor(yv(ind),'scale','log','quantiles',Q);
    
    model(i).x = xmid(i);
    model(i).mu = tmp.mu;
    model(i).sigma = tmp.sigma;
    model(i).sigma_mu = tmp.sigma_mu;
    model(i).rms = tmp.rms;
    
    yq(:,i) = yqtmp(:,2);
        
    % number of data
    model(i).N = sum(ind);
end

% turn the structure of scale parameters into a table for easier data handling
model = struct2table(model);


% Step 4 - produce smooth models for mu and sigma with field x to reduce noise in the
% distribution by using low-order polynomals
mind = ~isnan(model.mu);

% weight the fits by the number of data in each set of scale parameter
% determinations (this reduces the effects of bias from bins with little
% data)
W = diag(sqrt(model.N(mind)));

% quadratic fit to mu parameter
Aw2 = W*[ones(size(model.N(mind))) xmid(mind)' xmid(mind)'.^2];
Cw2 = inv(Aw2'*Aw2); 
smooth_model.mm = Cw2*Aw2'*W*model.mu(mind);
smooth_model.mm_alpha95 = 1.96*diag(Cw2).^(1/2);
smooth_model.mm_rms = sqrt((W*model.mu(mind) - Aw2*smooth_model.mm)'*(W*model.mu(mind) - Aw2*smooth_model.mm)/sum(mind));

% linear fit to sigma parameter
Aw = W*[ones(size(model.N(mind))) xmid(mind)'];
Cw = inv(Aw'*Aw); 
smooth_model.ms = Cw*Aw'*W*model.sigma(mind);
smooth_model.ms_alpha95 = 1.96*diag(Cw).^(1/2);
smooth_model.ms_rms = sqrt((W*model.sigma(mind) - Aw*smooth_model.ms)'*(W*model.sigma(mind) - Aw*smooth_model.ms)/sum(mind));

for i = 1:length(xedges)-1
    Yq(:,i) = smooth_model.mm(1) + smooth_model.mm(2)*xmid(i) + smooth_model.mm(3)*xmid(i)^2 ...
        + (smooth_model.ms(1) + smooth_model.ms(2)*xmid(i))*F;    
end

if isempty(xref)
    varargout{1} = xq;
    varargout{2} = yq;
    varargout{3} = Yq;
    varargout{4} = model;
    varargout{5} = smooth_model;
    return
end


% Step 5 - interpolate yq quantiles to reference xref value
for i = 1:length(Q)
    yref(i,1) = interp1(xq(i,:),yq(i,:),xref); % this is for plotting
    Yref(i,1) = interp1(xq(i,:),Yq(i,:),xref);
end

% find at y quantile for given x
warning off;
ind = ~isnan(Yq);

% produce a function to interpolate quantiles in Q with known positions
% (xq, Yq) to any x position.
QQ = repmat(Q,length(xedges)-1,1);
F = scatteredInterpolant(xq(ind),Yq(ind),QQ(ind));
warning on;

% find censored data and assume a below detection limit of 10 times the
% smallest detected value;
ind = y <= 0;
ya = y;
ya(ind) = -y(ind);
ya(y == 0) = 10*min(y(y > 0));

% interpolate data to the smoothed quantiles
Qy = F(x,log(ya));

% interpolate y quantile to xref along a quantile
y_at_ref = exp(interp1(Q,Yref,Qy));

% fix censored values
y_at_ref(ind) = -y_at_ref(ind);


% Step 6 - determine scale factors of adjusted data
ymodel = gausscensor(y_at_ref,'scale','log');

varargout{1} = y_at_ref;
varargout{2} = Qy;
varargout{3} = ymodel;
varargout{4} = yref;
varargout{5} = Yref;
varargout{6} = Q;
varargout{7} = xq;
varargout{8} = yq;
varargout{9} = Yq;
varargout{10} = model;
varargout{11} = smooth_model;

return