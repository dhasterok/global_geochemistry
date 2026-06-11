function [m,mint,stats] = hpVplot2(V,A,dV,varargin);
% HPVPLOT2 - produces a smoothed windowed quantile plot.
%
%   [m,mint,stats] = hpVplot2(X,Y,dX) produces a plot of X vs. log10(Y) in
%   equal steps in x.  Returns a linear model M with 95% confidence intervals
%   on fit (MINT) and STATS (r^2 f-statistic p-value.
%
%   [m,mint,stats] = hpVplot2(X,Y,dX,OPTION1,VALUE1,OPTION2,VALUE2,...)
%   allows setting a number of options:
%   
%   Defaults:
%       XLim   : [5.8 8.2]
%       XLabel : 'P-Velocity [km s^{-1}]'
%       XScale : 'linear'
%       XShift (for fitting) : -6
%
%       YLim   : [-2 1]
%       YLabel : 'Heat Production [\muW m^{-3}]'
%       YScale : 'log'
%
%       XFit   : none (range for fitting)...no fitting
%
%       Color  : [0 0.447 0.741]

% set default values
xrange = [5.8 8.2];
colour = [0 0.447 0.741];

xtxt = 'P-Velocity [km s^{-1}]';
xscale = 'linear';
xshift = -6;

ytxt = 'Heat Production [\muW m^{-3}]';
yscale = 'log';
yrange = [-2 1];
% read input options
opt = 1;
while opt + 3 < nargin
    switch lower(varargin{opt})
        case 'color'
            colour = varargin{opt+1};
        case 'xlim'
            xrange = varargin{opt+1};
        case 'xscale'
            xscale = varargin{opt+1};
        case 'xlabel'
            xtxt = varargin{opt+1};
        case 'xshift'
            xshift = varargin{opt+1};
        case 'yscale'
            yscale = varargin{opt+1};
        case 'ylabel'
            ytxt = varargin{opt+1};
        case 'ylim'
            yrange = varargin{opt+1};
        case 'xfit'
            xfit = varargin{opt+1};
        otherwise
            error(['Unknown option ',varargin{i}]);
    end
    opt = opt + 2;
end

if strcmpi('log',xscale)
    V = log10(V);
end

if strcmpi('log',yscale)
    A = log10(A);
end

ind = xrange(1) <= V & V <= xrange(2);
[Vm,Aq,N] = qrunavg(V(ind),A(ind),dV,colour)

if length(N) < 2
    m = NaN;
    merr = NaN;
    C = NaN;
    return
end

hold on;
ind = xfit(1) <= Vm & Vm <= xfit(2);

tmp = corrcoef(Vm(ind),Aq(ind,3));
X = [Vm+xshift ones(size(Vm))];
if length(tmp(:)) == 4
    C = tmp(1,2);
    [m,mint,r,rint,stats] = regress(Aq(ind,3),X(ind,:),0.95);

    stats(3); % p-value
    
else
    m  = NaN;
    merr = NaN;
end

subplot(3,1,1:2); hold on;
p = plot(Vm(ind),m(1)*(Vm(ind)+xshift) + m(2),'-');

C = colour-0.3;
C(C<0) = 0;
set(p,'Color',C);

if strcmpi(xscale,'log');
    hpax(xrange,'x');
else
    xlim(xrange);
end
if strcmpi(yscale,'log');
    hpax(yrange,'y');
else
    ylim(yrange);
end
ylabel(ytxt);
text(xrange(1)+0.1,log10(0.045),['r^2 = ',num2str(stats(1))]);
text(xrange(1)+0.1,log10(0.03),['log_{10} A [\muW m^{-3}]= ',num2str(m(1)),' (Vp[km s^{-1}] + ',num2str(xshift),') + ',num2str(m(2))]);
text(xrange(1)+0.1,log10(0.015),['m_{\alpha_{95}} = (',num2str(rint(1,1)),', ',num2str(rint(1,2)),'), b_{\alpha_{95}} (',num2str(rint(2,1)),', ',num2str(rint(2,2)),')']);

subplot(3,1,3); hold on;
ylabel('N');
if strcmpi(xscale,'log');
    hpax(xrange,'x');
else
    xlim(xrange);
end
xlabel(xtxt);

return