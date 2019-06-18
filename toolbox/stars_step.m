function [] = stars_step(data,len)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = length(data);

% compute t-statistic for 2*len-2 degrees of freedom
t = tstat(2*len - 2);

i = 1;
j = i + len;
pts = []
while j < N
    % compute the average running variance for windows of size len
    varl = avg_var(data,i,j,len);

    % compute difference considered significant
    sig_diff = t*sqrt(2*varl^2/len);

    R1 = mean(data(i:j-1));
    R2 = mean(data(j:j+len-1))
    
    T = R1 + [-1 1]*sig_diff;
    
    if R2 > T(2)
        
    elseif R2 < T(1)
        
        RSI = cumsum(data() - T
        j = i + len;
        pts = [pts; j];
    else
        i = i + 1;
        j = j + 1;
    end
end

return

function varl = avg_var(data,i,j,len);

for k = i:j-len+1;
    var(k) = std(data(k:k+len)).^2;
end

varl = mean(var);

return

