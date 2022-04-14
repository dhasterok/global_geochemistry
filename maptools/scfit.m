function [mu,rms] = scfit(E,P);

delta = sphangle(P(:,1),P(:,2),E(1),E(2));

mu = mean(delta);

N = length(P(:,1));
rms = sqrt(sum((mu*ones(size(delta)) - delta).^2)/N);

return
