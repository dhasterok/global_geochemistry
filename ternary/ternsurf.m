function varargout = ternsurf(a,b,c,val,dt,varargin)
% TERNSUF - Creates surface plot for ternary diagram.
%
%    TERNSURF(A,B,C,val,DT) computes a surface using Delaunay triangulation
%    with a discritization of DT on the ternary axis.  A, B, and C
%    represent the compositional triplets.
%
%    [T,V] = TERNSURF(A,B,C,val,DT) returns the centers of the triangular
%    coordinates, T, and the estimated values of the mesh, V.

% create ternary verticies to be used in triangulation
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

% determine vertex connections by Delaunay triangluation
tri = delaunay(xv,yv);

% The triangles should be all be the same area (within tolerance). However,
% some additional triangles may be created as a result of roundoff errors.
% It looks like these occur on the edges so we'll use area to remove them
% since they should have near zero area.
%
% If this code breaks it is most likely in this section (try changing
% tolerance?).
%----------------------------
area = heron(dt,dt,dt);

triarea = heron([xv(tri(:,1)), yv(tri(:,1))], ...
    [xv(tri(:,2)), yv(tri(:,2))], ...
    [xv(tri(:,3)), yv(tri(:,3))]);

ind = (abs(area - triarea) < (1e-3*dt)^3);
tri = tri(ind,:);

nt = length(tri(:,1));
if nt ~= round(1/dt^2)
    warning(['Incomplete tesselation (',num2str(1/dt^2),'/',num2str(nt),').']);
end

%figure(100);
%triplot(tri,xv,yv);
%----------------------------

% convert ternary points to Cartesian ignoring points with nan values
ind = (~isnan(val) & ~isnan(a) & ~isnan(b) & ~isnan(c));% ...
%    & (a < 100 & b < 100 & c < 100 & a > 0 & b > 0 & c > 0);
val = val(ind);
[xp,yp] = tern2xy(a(ind),b(ind),c(ind));

% compute statistics for each triangle
%figure;
for i = 1:nt
    ind = inpolygon(xp,yp,xv(tri(i,:)),yv(tri(i,:)));
    hold on;
    
    q(i,:) = quantile(val(ind),[0.05 0.25 0.5 0.75 0.95]);
    % uncomment below to check values
    plot([xv(tri(i,:));xv(tri(i,1))],[yv(tri(i,:));yv(tri(i,1))],'-');
    n(i,1) = sum(ind);
    %pause;
    %text(mean(xv(tri(i,:))),mean(yv(tri(i,:))),num2str(sum(ind)));
    if isnan(q(i,3))
        continue;
    end
    text(mean(xv(tri(i,:))),mean(yv(tri(i,:))),num2str(q(i,3)));
    
end
%n(n < 2) = NaN;

% place statistic at centroid
xc = sum(xv(tri),2)/3;
yc = sum(yv(tri),2)/3;
%[q(:,3) (q(:,4) - q(:,2))*1.35]
m = scatteredInterpolant(xc,yc,q(:,3),'linear');
s = scatteredInterpolant(xc,yc,(q(:,4) - q(:,2))*1.35,'linear');
N = scatteredInterpolant(xc,yc,log10(n),'linear');

out.tri = tri;
out.xv = xv;
out.yv = yv;
out.median = m(xv,yv);
out.sd = s(xv,yv);
out.nbin = N(xv,yv);

if nargout > 0
    varargout{1} = out;
    return
end

figure;
subplot(221);
trisurf(tri,xv,yv,m(xv,yv));
%shading interp;
shading flat;
set(gca,'View',[0 90]);
hold on;
if extend > 0;
    ternextend(extend);
end
ternary;
%plot(xc,yc,'c.');
%for i = 1:length(xc)
%    text(xc(i),yc(i),num2str(m(xc(i),yc(i))));
%end
hold off;

subplot(222);
trisurf(tri,xv,yv,s(xv,yv));
%shading interp;
shading flat;
set(gca,'View',[0 90]);
hold on;
if extend > 0;
    ternextend(extend);
end
ternary;
hold off;

subplot(223);
trisurf(tri,xv,yv,N(xv,yv));
%shading interp;
shading flat;
set(gca,'View',[0 90]);
hold on;
if extend > 0;
    ternextend(extend);
end
ternary;
hold off;


return