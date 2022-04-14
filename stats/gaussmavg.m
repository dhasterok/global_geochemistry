function [t,y,n] = gaussmavg(tobs,yobs,dt,nstep)

if mod(nstep,2) ~= 0
    nstep = nstep + 1;
end

t = [floor(min(tobs)/dt)*dt:dt:ceil(max(tobs)/dt)*dt]';
y = nan(length(t)-nstep+1,5);
n = zeros([length(t)-nstep+1 1]);
for i = nstep/2:length(t)-nstep/2
    ind = t(i-nstep/2+1) <= tobs & tobs <= t(i+nstep/2);
    n(i-nstep/2+1) = sum(ind);
    
    %tt(i-nstep/2+1,1) = (t(i-nstep/2+1) + t(i+nstep/2))/2;
    if sum(ind) < 10
        continue;
    end
    
    [mu,sigma,mu_95,sigma_95] = normfit(yobs(ind),0.05);
    [~,~,mu_68,sigma_68] = normfit(yobs(ind),0.32);

    y(i-nstep/2+1,:) = [mu_95(1), mu_68(1), mu, mu_68(2), mu_95(2)];
end

t = midpt(t);
t = t(nstep/2:end-nstep/2+1);

%ind = ~isnan(y(:,3));
%y = y(ind,:);


return