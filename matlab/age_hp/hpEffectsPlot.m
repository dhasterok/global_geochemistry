clear all
%close all
clc

figure()
x = [0:10:4000];

% Reference line:
y_ref = ones(size(x));

% Erosion
y_ero = 1 - 0.75.*(1 - exp(-0.001038605.*x));
y_ero = 1 - 0.65.*(1 - exp(-0.0025.*x));

% Secular Cooling
y_sec = 1 - 0.65.*(1 - exp(-0.001038605.*x));

% Reworking and tectonics
y_tec = (1 - 0.675.*(1 - exp(-0.0003.*x))) .* (1+(0.05 + 0.1.*(x./4000)).*sind(0.3.*(2-0.8.*x./4000).*x));
y_tec = (1 - 0.675.*(1 - exp(-0.0003.*x))) .* (1+(0.05 + 0.1.*(x./4000)).*sind(0.37.*x));


% Thermal stability
y_them = -476026.2 + (0.9771206 - -476026.2)./(1 + (x./138593.1).^3.780061);
y_them = -198829.1 + (0.9852472 - -198829.1)./(1 + (x./258457.4).^3.053539);
y_them = -198829.088 + (0.9852472 - -198829.1)./(1 + (x./258457.4).^3.053539);

y_them = -198829.088 + (0.9852472 - -198829.1)./(1 + (x./258457.4).^3.2);


% Mantle depletion
y_dep = 0.007.*exp(0.001.*x + 0) + 1 - 0.008;

% Decay adjust curve:
y_dec = 0.03.*exp(0.0007.*x + 0) + 1 - 0.03;
%y_dec = 0.04.*exp(0.0005.*x + 0) + 1 - 0.03;

subplot(1,2,1)
hold on
plot(x,y_ref,'--k')
%plot(x,y_dec,'--k')
plot(x,y_ero,'-g')
plot(x,y_sec,'-r')
plot(x,y_tec,'-k')
plot(x,y_them,'-m')
plot(x,y_dep,'-b')
%plot(x,y_dec,'--b')
hold off

xlim([0 4000])
ylim([0 2])
axis square


subplot(1,2,2)
% Decay adjusted
hold on
plot(x,y_ref.*y_dec,'--k')
plot(x,y_ero.*y_dec,'-g')
plot(x,y_sec.*y_dec,'-r')
plot(x,y_tec.*y_dec,'-k')
plot(x,y_them.*y_dec,'-m')
plot(x,y_dep.*y_dec,'-b')
hold off

xlim([0 4000])
ylim([0 2])
axis square
