function [A1,A2,A3] = mclass_diff(data,verticies,wt2mol,indp,indq)

a1 = verticies{1};
a2 = verticies{2};
a3 = verticies{3};

% check number of inputs and create an initial index of logicals if
% necessary
ind0 = data.heat_production > 0;

% parse data
A1 = parsea(data,a1,wt2mol);
A2 = parsea(data,a2,wt2mol);
A3 = parsea(data,a3,wt2mol);
%B = A2 - (data.NA2O/molecularwt('Na2O') + data.K2O/molecularwt('K2O'));
%[min(B) max(B) mean(B(~isnan(B))) std(B(~isnan(B)))]


figure;
subplot(221);
% Create ternary diagram
T = sum(A1 + A2 + A3,2);

A1 = A1./T;
A2 = A2./T;
A3 = A3./T;

ind = ~isnan(T) & ind0;

if sum(A1 < 0) > 0
    ternextend(0.5);
end
hold on;
ternary(a1,a2,a3);
t = ternplot(A1(ind & indp),A2(ind & indp),A3(ind & indp),'.');
set(t,'MarkerSize',2);
t = ternplot(A1(ind & indq),A2(ind & indq),A3(ind & indq),'.');
set(t,'MarkerSize',2);


if sum(A1 < 0) > 0
    outp = ternsurf(A1(ind & indp),A2(ind & indp),A3(ind & indp),log10(data.heat_production(ind & indp)),0.05,0.5);
    outq = ternsurf(A1(ind & indq),A2(ind & indq),A3(ind & indq),log10(data.heat_production(ind & indq)),0.05,0.5);
else
    outp = ternsurf(A1(ind & indp),A2(ind & indp),A3(ind & indp),log10(data.heat_production(ind & indp)),0.05);
    outq = ternsurf(A1(ind & indq),A2(ind & indq),A3(ind & indq),log10(data.heat_production(ind & indq)),0.05);
end

% Number of data
subplot(222);
trisurf(outp.tri,outp.xv,outp.yv,outp.nbin);
%shading interp;
shading flat;
set(gca,'View',[0 90]);
hold on;
ternary(a1,a2,a3);
caxis([0 3]);
cbar('No. sed.');

subplot(223);
trisurf(outq.tri,outq.xv,outq.yv,outq.nbin);
%shading interp;
shading flat;
set(gca,'View',[0 90]);
hold on;
ternary(a1,a2,a3);
caxis([0 3]);
cbar('No. metased.');
%cbar('\sigma_{\mu}(N-1)^{-1/2}');

% differences in median values
figure;
subplot(223);
histogram(outp.median - outq.median,'BinEdges',[-2:0.2:2]);
xlim([-2 2]);
xlabel('P - Q');


subplot(224);
trisurf(outp.tri,outp.xv,outp.yv, outp.median - outq.median);
%shading interp;
shading flat;
set(gca,'View',[0 90]);
hold on;
ternary;
ternary(a1,a2,a3);
caxis([-1 1]);
colormap(rwb);
cbar('log_{10} sed. - log_{10} metased.');



return


