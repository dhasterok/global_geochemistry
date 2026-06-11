function [A1,A2,A3] = mclass(data,verticies,wt2mol,varargin)

a1 = verticies{1};
a2 = verticies{2};
a3 = verticies{3};

% check number of inputs and create an initial index of logicals if
% necessary
field = 'heat_production';
if nargin == 4
    field = varargin{1};
    if ~any(strcmpi(field,data.Properties.VariableNames)) & ~strcmpi(field,'none')
        error('Field not found in input data.');
    end
elseif nargin ~= 3
    help mclass;
    error('Incorrect number of arguments.');
end

if strcmp(field,'heat_production')
    ind0 = data.heat_production > 0;
else
    ind0 = logical(ones([height(data) 1]));
end

% parse data
A1 = parsea(data,a1,wt2mol);
A2 = parsea(data,a2,wt2mol);
A3 = parsea(data,a3,wt2mol);
%B = A2 - (data.NA2O/molecularwt('Na2O') + data.K2O/molecularwt('K2O'));
%[min(B) max(B) mean(B(~isnan(B))) std(B(~isnan(B)))]


figure;
% Create AFM diagram
T = sum(A1 + A2 + A3,2);

ind = ~isnan(T) & ind0;
A1 = A1(ind)./T(ind);
A2 = A2(ind)./T(ind);
A3 = A3(ind)./T(ind);

if strcmp(field,'heat_production')
    D = log10(data.heat_production(ind));
elseif strcmpi(field,'none')
    D = [];
else
    D = data{ind,field};
end

if sum(A1 < 0) > 0
    ternextend(0.5);
end
hold on;
ternary(a1,a2,a3);
if isempty(D)
    t = ternplot(A1,A2,A3,'.');
    return
else
    t = ternplot(A1,A2,A3,D);
end
set(t,'MarkerSize',2*ones(size(get(t,'MarkerSize'))));
caxis([-1 1]);
cbar('log_{10} A');


if sum(A1 < 0) > 0
    ternsurf(A1,A2,A3,D,0.05,0.5);
else
    ternsurf(A1,A2,A3,D,0.05);
end
subplot(221);
ternary(a1,a2,a3);
caxis([-1 1]);
cbar('log_{10} A');
%colormap(hot);

subplot(222);
ternary(a1,a2,a3);
caxis([0 2]);
cbar('\sigma_{\mu}');
%cbar('\sigma_{\mu}(N-1)^{-1/2}');

subplot(223);
ternary(a1,a2,a3);
caxis([0 3]);
cbar('N');
%cbar('\sigma_{\mu}(N-1)^{-1/2}');

return