function [heat,heat_rate] = radheat;

isotope = {'k40','rb87','sm148','th232','u235','u238'};
%isotope = {'th232'};


for i = 1:length(isotope)
    T = readtable([isotope{i},'_chain.xlsx']);
    for j = 1:length(T.half_life)
        switch T.unit{j}
        case 'y'
            T.half_life(j) = T.half_life(j)*365.2422*24*3600; % s/yr
            T.half_life_uncertainty(j) = T.half_life_uncertainty(j)*365.2422*24*3600;
        case 'd'
            T.half_life(j) = T.half_life(j)*24*3600; % s/day
            T.half_life_uncertainty(j) = T.half_life_uncertainty(j)*24*3600;
        case 'h'
            T.half_life(j) = T.half_life(j)*3600; % s/hr
            T.half_life_uncertainty(j) = T.half_life_uncertainty(j)*3600;
        case 'm'
            T.half_life(j) = T.half_life(j)*60; % s/min
            T.half_life_uncertainty(j) = T.half_life_uncertainty(j)*60;
        case 'ms'
            T.half_life(j) = T.half_life(j)*10^-3; % s/ms
            T.half_life_uncertainty(j) = T.half_life_uncertainty(j)*10^-3;
        case 'us'
            T.half_life(j) = T.half_life(j)*10^-6; % s/us
            T.half_life_uncertainty(j) = T.half_life_uncertainty(j)*10^-6;
        case 'ns'
            T.half_life(j) = T.half_life(j)*10^-9; % s/ns
            T.half_life_uncertainty(j) = T.half_life_uncertainty(j)*10^-9;
        end
    end

    T.probability = T.probability/100;
    T.probability_uncertainty = T.probability_uncertainty/100;

    ind = T.probability == 1;
    T.probability_uncertainty(ind) = 0;

    heat(i,:) = construct_chain(0,T);
    heat_rate(i,:) = construct_chain(0,T,[0 0]);
end

heat(:,2) = sqrt(heat(:,2));
heat_rate(:,2) = sqrt(heat_rate(:,2));

return

function heat = construct_chain(i,T,varargin)

% start with parent
% determine number of daughter products
% for each daughter make parent
%   repeat

if i == 0
    indp = find(strcmp(T.parent,T.parent(1)));

    for j = 1:length(indp)
        if nargin == 2
            tmpheat(j,:) = heat_lookup(indp(j), T, construct_chain(indp(j),T));
        else
            disp('start')
            half_life(1) = T.half_life(indp(j));
            half_life(2) = T.half_life_uncertainty(indp(j));
            tmpheat(j,:) = heat_lookup(indp(j), T, construct_chain(indp(j),T,half_life), half_life);
        end
    end

    heat = sum(tmpheat,1);
    return
end

ind = find(strcmp(T.parent,T.daughter(i)));
if length(ind) == 0
    heat = [0 0];
    return
end

for j = 1:length(ind)
    if nargin == 2
        tmpheat(j,:) = heat_lookup(ind(j), T, construct_chain(ind(j),T));
    else
        half_life(1) = varargin{1}(1) + T.half_life(ind(j))
        half_life(2) = sqrt(varargin{1}(2)^2 + T.half_life_uncertainty(ind(j))^2);
        tmpheat(j,:) = heat_lookup(ind(j), T, construct_chain(ind(j),T,half_life), half_life);
    end
end
heat = sum(tmpheat,1);

return


function [Q,dQ] = heat_lookup(i, T, Q, varargin)

if nargin == 4
    t = varargin{1};
    lambda(1) = log(2)/t(1);
    lambda(2) = log(2)*t(2)/t(1)^2;
end

if ~isnan(T.average_beta(i))
    tmp(1) = T.average_beta(i);
    tmp(2) = T.average_beta_uncertainty(i);
else
    tmp(1) = T.total_energy(i);
    tmp(2) = T.total_energy_uncertainty(i);
end
if nargin < 4
    Q(1) = T.probability(i) * (tmp(1) + Q(1));
    Q(2) = T.probability_uncertainty(i)^2 * (tmp(1) + Q(1))^2 + ...
        T.probability(i)^2 * (tmp(2)^2 + Q(2));
else
    Q(1) = T.probability(i) * (tmp(1) * lambda(1) + Q(1)) ;
    Q(2) = T.probability_uncertainty(i)^2 * (tmp(1) * lambda(1) + Q(1))^2 + ...
        T.probability(i)^2 * (tmp(2)^2 * lambda(1)^2 + tmp(1)^2 * lambda(2)^2 + Q(2));
end

return
