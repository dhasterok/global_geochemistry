function [m,merr,C,Vedges] = hpVplot(V,A,varargin)
% [m,merr,C] = hpVplot(V,A,Vp,colour);

field = 'p_velocity';
colour = [54 169 225]/255;
opt = 1;
Vedges = [];
while opt + 2 < nargin
    switch lower(varargin{opt})
        case 'field'
            field = varargin{opt+1};
        case 'color'
            colour = varargin{opt+1};
        case 'binedges'
            Vedges = varargin{opt+1};
        otherwise
            error('Unknown option.');
    end
    opt = opt + 2;
end

switch field
    case 'p_velocity'
        xlbl = 'Estimated P-Velocity [km s^{-1}]';
        xl = [5.8 8.2];
        xfit = [6.0 7.4];
        Vshift = 6;
        if isempty(Vedges)
            Vedges = [5.8:0.1:8.2];
        end
    case 's_velocity'
        xlbl = 'Estimated S-Velocity [km s^{-1}]';
        xl = [3.65 4.6];
        xfit = [3.65 4.05];
        Vshift = 3.65;
        if isempty(Vedges)
            Vedges = [3.6:0.05:4.6];
        end
    case 'vmoho'
        xlbl = 'Crustal Thickness [km]';
        xl = [0 75];
        xfit = xl;
        Vshift = 0;
        if isempty(Vedges)
            Vedges = [0:5:75];
        end
end

for i = 1:length(Vedges)-1
    ind = Vedges(i) <= V & V < Vedges(i+1);
    
    Atemp{i} = A(ind);
    Vtemp{i} = V(ind);
    M(i,1) = sum(ind);
end

N = 0;
for i = 1:length(Atemp)
    N = N + length(Atemp{i});
end

if N < 2
    m = NaN;
    merr = NaN;
    C = NaN;
    return
end
[X,Y] = whisker(Vtemp,Atemp,'Color',colour,'Scale','log');

% interquartile range
iqr_x = X(:,4) - X(:,2);
iqr_y = Y(:,4) - Y(:,2);

% estimated 1 standard error about mean
stderr_x = iqr_x./(1.35*sqrt(M));
stderr_y = iqr_y./(1.35*sqrt(M));

hold on;
%ind = (6 <= X(:,3) & X(:,3) <= 7.4) & M >= 9;
ind = (xfit(1) <= X(:,3) & X(:,3) <= xfit(2)) & ~isnan(Y(:,3)) & M >= 9;

%tmp = corrcoef(X(ind,3),Y(ind,3));
if sum(ind) >= 3
    %C = tmp(1,2);

    [mtmp,r,rsq] = wllsq(X(ind,3)-Vshift,Y(ind,3),stderr_x(ind),stderr_y(ind));
    
    % compute correlation coefficient
    C = sign(mtmp(1,1))*sqrt(rsq);
    
    
    
    merr = mtmp(:,2);
    m = mtmp(:,1);
    %[m,s] = polyfit(X(ind,3)-6,Y(ind,3),1);

    %invR = inv(s.R);
    %merr = sqrt(diag(invR*invR')*s.normr^2/s.df);
else
    m  = NaN;
    merr = NaN;
end

p = plot(Vedges,polyval(m,Vedges-Vshift),'-');
set(p,'Color',colour);

hpax([-2 1],'y');
xlabel(xlbl);
xlim(xl);
set(gca,'Box','on');
golden;

return
  
