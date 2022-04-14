function stats = basicstats(v);

stats.ndata = length(v);
stats.min = min(v);
stats.max = max(v);
stats.mean = mean(v);
stats.stdev = std(v);
q = quantile(v,[0.025 0.25 0.5 0.75 0.975]);
stats.q025 = q(1);
stats.q250 = q(2);
stats.q500 = q(3);
stats.q750 = q(4);
stats.q975 = q(5);

return