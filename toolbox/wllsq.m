function [m,r,rsq,Cw] = wllsq(x,y,sigx,sigy)
% WLLSQ - weighted least squares fitting
%
%   [m,r,rsq] = wllsq(x,y,sigx,sigy) will use weighted least squares to fit
%   a linear model to a set of (x,y) data with uncertainties in x and y.
% 
%   Returns model parameters m = (slope alpha95; intercept alpha95), r
%   signed model residuals, rsq is the r^2 value.

% D. Hasterok 2018

% first iteration, compute linear least squares with y weights... ignore
% the x weights.
if nargin == 2
    x = x(:);
    y = y(:);
    
    A = [x(:) ones(size(x(:)))];
    Cw = inv(A'*A);
    whos
    m = Cw*A'*y;
    r = y - A*m;
   
    m_a95 = 1.96*diag(Cw).^(1/2);
    
    % de-mean y
    %l = ones(size(y));
    %yd = y - l*inv(l'*l)*l'*y;
        
    % compute r^2 value
    %rsq = 1 - r'*r/(yd'*yd);
    rsq = sum((A*m - mean(y)).^2)/sum((y - mean(y)).^2);
        
    m = [m_new m_a95];
end

if (size(x,1) == 1 || size(x,2) == 1) && (size(y,1) == 1 || size(y,2) == 1)
    x = x(:);
    y = y(:);
    
    sigx = sigx(:);
    sigy = sigy(:);
    
    A = [x(:) ones(size(x(:)))];
elseif (size(x,1) ~= 1 && size(x,2) ~= 1) && (size(y,1) == 1 || size(y,2) == 1)
    A = x;
    y = y(:);
    
    sigx = sigx(:);
    sigy = sigy(:);

else
    ndims(x)
    ndims(y)
    error('Not sure what to do, check input dimensions.');
end

Wy = diag(1./sigy);

Aw = Wy*A;
yw = Wy*y;

m_new = inv(Aw'*Aw)*Aw'*yw;
r0 = yw - Aw*m_new;

tol = 1e-4;
c = 1;
maxit = 1;
while 1
    m_old = m_new;
    
    Wyx = diag(1./(sigy + m_old(1)^2*sigx));
    
    Aw = Wyx*A;
    yw = Wyx*y(:);
    
    % covariance matrix
    Cw = inv(Aw'*Aw);
    
    % compute weighted linear model
    m_new = Cw*Aw'*yw;
    
    % compute residuals
    r = yw - Aw*m_new;
    
    %fprintf('%f %f %f %f\n',m_new(1),m_new(2),norm(r),norm(yw));
    
    if norm(r)/norm(r0) < tol | c > maxit;
        
        m_a95 = 1.96*diag(Cw).^(1/2);
        % de-mean yw
        %l = ones(size(yw));
        %ywd = yw - l*inv(l'*l)*l'*yw;
        
        % compute r^2 value
        %rsq = 1 - r'*r/(ywd'*ywd);
        rsq = sum((Aw*m_new - mean(yw)).^2)/sum((yw - mean(yw)).^2);
        break;
    end
    c = c + 1;
end

m = [m_new m_a95];

return