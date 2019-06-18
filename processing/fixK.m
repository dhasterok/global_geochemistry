function data = fixK(data);

mK2O = 2*39.0986 + 15.9994;
mK = 39.0986;

% check K vs K2O
if ~isfield(data,'K2O') & ~isfield(data,'K_PPM')
    return;
elseif ~isfield(data,'K_PPM');
    data.K_PPM = data.K2O * 2*mK/mK2O*1e4;
    return;
end

for i = 1:length(data.K_PPM)
    K = data.K_PPM(i);
    if isnan(K)
        continue;
    end
    if ~isnan(data.K2O(i))
        Ko = data.K2O(i) * 2*mK/mK2O*1e4;
        if K/Ko < 0.5 | K/Ko > 2
            fprintf('@ index (%i) [Kox,K] = [%f,%f]\n',i,Ko,K);
        else
            data.K_PPM(i) = K;
        end
    end    
end

ind = find(data.K2O > 0 & data.K_PPM > 0 & ~isnan(data.K2O) & ~isnan(data.K_PPM));

figure;
subplot(121);
plot(log10(data.K2O(ind)),log10(data.K_PPM(ind)),'.');
hold on;
xlabel('K_2O [wt.%]');
ylabel('K [ppm]');

K2O = [1e-3 100];
plot(log10(K2O),log10(K2O * 2*mK/mK2O*1e4),'-');

subplot(122);
histogram(log10(data.K2O(~isnan(data.K2O) & data.K2O > 0) * 2*mK/mK2O*1e4));
hold on;
histogram(log10(data.K_PPM(~isnan(data.K_PPM) & data.K_PPM > 0)));
xlabel('log_{10} K [ppm]');
legend('K_2O as K','K');

return
