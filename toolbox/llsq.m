function [m,r,rms,rsq] = llsq(yd,A);

% covariance matrix
C = inv(A'*A);

% model
m = C*A'*yd;
m_a95 = 1.96*diag(C).^(1/2);
m = [m m_a95];

y = A*m(:,1);

% residual
r = yd - y;
l = ones(size(yd));

% RMS misfit
rms = sqrt(sum(r.^2/length(yd)));

% zero mean
rsq = sum((y - mean(yd)).^2)/sum((yd - mean(yd)).^2);
%yd0 = y - mean(y);

% compute r^2 value
% this didn't work
%rsq = 1 - r'*r/(yd0'*yd0);
%mss = sum((A*m(:,1) - mean(A*m(:,1))).^2);
%rss = sum(r.^2);

%rsq = mss/(mss + rss);

return