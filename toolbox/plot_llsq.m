function p = plot_llsq(x,y,m,varargin);
% plots a linear regression model from LLSQ and confidence intervals

alpha = 0.95;
mlim = [min(x) max(x)];

if size(x,2) == 2
    w = x(:,2);
    x = x(:,1);
else
    w = ones(size(x));
end

c = 1;
if nargin > 3
    while c <= nargin-3
        switch length(varargin{c})
            case 1
                alpha = varargin{c};
            case 2
                mlim = varargin{c};
            otherwise
                error('Input not recognized.');
        end
        c = c + 1;
    end
end

N = length(x);

xtmp = linspace(mlim(1),mlim(2),50)';
ytmp = [ones(size(xtmp)) xtmp]*m(:,1);

xmu = sum(w.*x)/sum(w);

se = sqrt(sum((y - [ones(size(x)) x]*m(:,1)).^2)/(N - 2)) * ...
    sqrt(1/N + (xtmp - xmu).^2/sum((x - xmu).^2));

yci = repmat(ytmp,1,2) + ...
    tinv((1 - alpha)/2,N-2)*repmat([-1 1],length(xtmp),1).*repmat(se,1,2);

p{2} = fill([xtmp; flipud(xtmp)],[yci(:,1); flipud(yci(:,2))],[0.7 0.7 0.7], ...
    'EdgeColor','none');
p{1} = plot(mlim', [1 1; mlim]'*m(:,1), 'k-');

%plot(xtmp,m(1,1)+m(1,2) + (m(2,1)+m(2,2))*xtmp);
%plot(xtmp,m(1,1)+m(1,2) + (m(2,1)-m(2,2))*xtmp);
%plot(xtmp,m(1,1)-m(1,2) + (m(2,1)+m(2,2))*xtmp);
%plot(xtmp,m(1,1)-m(1,2) + (m(2,1)-m(2,2))*xtmp);

return