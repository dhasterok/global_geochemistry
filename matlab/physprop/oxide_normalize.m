function t = oxide_normalize(t)

global OXIDES

oxides = lower(OXIDES);

t{:,oxides} = 100 * t{:,oxides}./repmat(sum(t{:,oxides},2),1,length(oxides));

return