function vq = prepdata_avg(data,el,scale)
% PREPDATA_AVG - prepares average value for plotting
%
%   vq = prepdata_avg(data,el,scale) computes the quantile levels for
%   plotting the 'average' data using a gaussian CDF model taking into
%   account left-censored data.  el is a field of data and scale is either
%   'linear' or 'log'.

if iscell(el) & length(el) == 2
    % to handle censoring it would require dealing with BDL values in the
    % denominator resulting in right censoring of ratios, or what to do
    % with BDL/BDL values?  Best to skip handling them for now.
    ind = data{:,el{1}} > 0 & data{:,el{2}} > 0;
    v = data{ind,el{1}}./data{ind,el{2}};
else
    ind = data{:,el} > 0;
    v = data{ind,el};
end

switch scale
    case 'linear'
        [model,vq] = gausscensor(v);
    case 'log'
        [model,vq] = gausscensor(v,'Scale','log10');
    otherwise
        error(['Unknown scale (',scale,').']);
end

return