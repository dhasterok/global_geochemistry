function [xv,yv] = gridtern(dt);

g = [0:dt:1];
n = length(g);
nv = 1;

extend = 0;
if nargin == 6
    extend = varargin{1};
end

for i = 1:n
    for j = 1:n
        for k = 1:n
            if abs(g(i) + g(j) + g(k) - 1) > 0.001*dt
                continue;
            else
                Tv(nv,:) = [g(i),g(j),g(k)];
                nv = nv + 1;
            end
        end
    end
end

if extend > 0
    ga = [-extend:dt:0-dt];
    gbc = [0:dt:1+extend];
    na = length(ga);
    nbc = length(gbc);
    for i = 1:na
        for j = 1:nbc
            for k = 1:nbc
                if abs(ga(i) + gbc(j) + gbc(k) - 1) > 0.001*dt
                    continue;
                else
                    Tv(nv,:) = [ga(i),gbc(j),gbc(k)];
                    nv = nv + 1;
                end
            end
        end
    end
end
Tv = flipud(Tv);
nv = nv - 1;

% convert ternary verticies to Cartesian
[xv,yv] = tern2xy(Tv(:,1),Tv(:,2),Tv(:,3));

return