close all;
clear all;

kcf = 1e4*2*molecularwt('K')/molecularwt('K2O');

k2o0 = 4.52; % wt.%
th0 = 17.7;  % ppm
u0 = 4.0;    % ppm

[~,~,Ciso,abundance] = radtime(1,k2o0,0,0,th0,u0,'K2O','Age',4000);

k2o = (k2o0*kcf*(1 - 0.000117) + Ciso(:,:,1))/kcf;
th = Ciso(:,:,4);
u = Ciso(:,:,5) + Ciso(:,:,6);

[~,~,Ciso] = radtime(1,k2o,0,0,th,u,'K2O','Age',-4000,'Abundance',4000);

tmp(1) = (k2o0*kcf*(1 - 0.000117) + Ciso(:,:,1))/kcf;

tmp(2) = Ciso(:,:,4);
tmp(3) = Ciso(:,:,5) + Ciso(:,:,6);

% this doesn't take into account u234 created by u238 decay properly, the
% error is tmp(3) - u0
(tmp - [k2o0 th0 u0])./[k2o0 th0 u0]

% take present day-zero age (eh...last 200 Ma) concentrations
k2o0 = 4.52; % wt.%
th0 = 17.7;  % ppm
u0 = 4.0;    % ppm

% compute concentrations today if a rock created X Ma ago had PD-ZA
% concentration at its point of crystallization
t = [-4000:50:0];
for i = 1:length(t)
    [A(i),~,Ciso] = radtime({'granite'},k2o0,0,0,th0,u0,'K2O','Age',t(i),'Abundance',-t(i));
    
    k(i) = (k2o0*kcf*(1 - 0.000117) + Ciso(:,:,1));
    th(i) = Ciso(:,:,4);
    u(i) = Ciso(:,:,5) + Ciso(:,:,6);
end
t = -t;

figure;
subplot(321);
plot(t,k);

subplot(323);
plot(t,u,t,th);


subplot(322);
plot(t,log10(k./u));
hpax([-0.5 1.5]);
hpax([2 6]);
subplot(324);
plot(t,log10(th./u));
hpax([0 1]);

subplot(325);
plot(t,A);