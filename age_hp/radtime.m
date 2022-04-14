function [A,A_iso,C_iso,abundance,C_el] = radtime(density,K,Rb,Sm,Th,U,varargin)
% RADTIME - Computes radioactivity given chemistry.
%
%   A = radtime(density,K,Rb,Sm,Th,U) computes the volumetric heat production
%   of a sample with K in wt.% and Rb, Sm, Th, and U in ppm.
%
%   A = radtime(rock_type,K,Rb,Sm,Th,U) if density is unknown a rock type may
%   be entered from the standard list below.  The rock type must be entered
%   as a cell array, (e.g. rock_type = {'granite','subalkalic basalt', ...}).
%   Note porosity is not taken into account in the estimated densities.  To
%   account for sample porosity, multiply porosity by the computed heat
%   production.
%
%   To obtain the contribution due to each isotope,
%   [A,A_iso,C_iso,abundance,C_el] = radtime(density,K,Rb,Sm,Th,U);
%
%
%   Several additional options are available
%   (note options are case insensitive):
%
%       'K2O': converts an input K2O in wt.% to K (i.e.
%           radtime(density,K2O,Rb,Sm,Th,U,'K2O');
%
%       'Age': compute heat production using the following chemistry
%              projected to age of formation (in Ma)
%              radtime(density,K,Rb,Sm,Th,U,'Age',t);
%
%       'TimeZero': using concentrations at TimeZero as the initial
%              concetration at time, t0 and then projected to Age, t.  If
%              not age is given, present day t = 0 is assumed.
%
%       'Formula': computes heat production using coefficients determined
%              by one of the following studies
%              'hg17' - Heat Generation Model, revision 2017 (HG.r2017) by
%                   Hasterok & Gard (EPSL, submitted)
%              'r88'  - Rybach (Handbook on Terrestrial Heat Flow Density, 1988)
%              'ts14' - Turcotte & Schubert (Geodynamics, 3rd ed., 2014)
%              'd15' - Dye (Geodynamics, 3rd ed., 2014)
%
%   Standard rock type list:
%       plutonic                     Volcanic
%           'alkalic gabbro' 'carbonatite' 'cumulate peridotite'
%           'diorite' 'foid gabbro' 'foid monzodiorite' 'foid monzosyenite'
%           'foid syenite' 'gabbroic diorite' 'granite' 'granodiorite'
%           'intermediate foidolite' 'mafic foidolite' 'mantle peridotite'
%           'monzodiorite' 'monzogabbro' 'monzonite'
%           'peridotgabbro' 'quartz monzonite' 'quartzolite'
%           'subalkalic gabbro' 'syenite' 'ultra-high alkali plutonic'
%           'ultramafic foidolite' 'alkali picrite' 'alkalic basalt'
%           'andesite' 'basaltic andesite' 'basaltic trachyandesite'
%           'boninite' 'carbonatite' 'dacite' 'intermediate' 'komatiite'
%           'mafic foidite' 'meimechite' 'phonolite' 'phonotephrite'
%           'picrite' 'picrobasalt' 'rhyolite' 'silexite' 'subalkalic basalt'
%           'tephriphonolite' 'tephrite' 'trachyandesite' 'trachybasalt'
%           'trachydacite' 'trachyte' 'ultra-high alkali volcanic'
%           'ultramafic foidite' 
%
%       sedimentary:
%           'quartzite'         'arkose'
%           'subarkose'         'wacke'
%           'limestone'         'dolomite'
%           'iron-rich sand'    'iron-rich shale'
%           'laterite'          'oxide'         
%           'quartz arenite'    'litharenite'
%           'sublitharenite'    'shale'
%

if nargin == 0
    help radtime
    return
end

% if rock types are given instead of densities, use lookup table.
if iscell(density) | ischar(density)
    rock_type = density;
    density = lookup_density(rock_type);
end

% set options for computing K, age, and the coefficients used to compute
% heat production
kopt = 0;
age = 0;
popt = 'hg17';
c = 1;
aflag = 0;
if nargin > 6
    while c <= nargin-6
        switch lower(varargin{c})
            case 'k' % K as K(wt.%)
                kopt = 0;
                c = c + 1;
            case 'k2o' % K as K2O(wt.%)
                kopt = 1;
                c = c + 1;
            case 'age' % age in Ma
                age = varargin{c+1};
                c = c + 2;
            case 'timezero' % time in Ma
                age0 = varargin{c+1};
                c = c + 2;
            case 'abundance' % time in Ma
                aflag = 1;
                t0 = varargin{c+1};
                c = c + 2;
            case 'formula' % sets heat production coefficients
                popt = varargin{c+1};
                if ~any(strcmp({'hg17','r88','ts14','d12'},popt))
                    warning('Unknown formula requested, using hg17.');
                    popt = 'hg17';
                end
                c = c + 2;

            otherwise
                error(['An unknown option was entered in argument ,', ...
                    num2str(6 + c),'.']);
        end
    end
end

if kopt
    % Convert from wt.% K2O to ppm K
    K = 2*molecularwt('K')/molecularwt('K2O') * K * 1e4;
else
    % convert wt.% K to ppm K
    K = K*1e4;
end

% heat production coefficients
% isotope order: K40, Rb87, Sm148, Th232, U235, U238
switch popt
case 'hg17' % Hasterok and Gard (2017)
    half_life = [1.248e3 4.97e4 7e6 1.40e4 7.04e2 4.468e3]; % Ma
    abundance = [0.000117 0.2783 0.1124 1 0.007204 0.992742];
    heat_prod = [2.7906e-5 4.0071e-8 7.3e-9 2.50289e-5 56.6777e-5 9.2840e-5];
case 'r88' % Rybach (1988)
    half_life = [1.23e3 4.97e4 7e6 1.39e4 7.13e2 4.51e3]; % Ma
    abundance = [0.000117 1 1 1 0.00711 0.9928];
    heat_prod = [2.9744e-5 0 0 2.56e-5 57.5e-5 9.17e-5];
case 'ts14' % Turcotte and Schubert (Geodynamics, 3rd ed., 2014)
    half_life = [1.25e3 4.97e4 7e6 1.40e4 7.04e2 4.47e3]; % Ma
    abundance = [0.000119 1 1 1 0.0071 0.9928];
    heat_prod = [2.92e-5 0 0 2.64e-5 56.9e-5 9.46e-5];
case 'd12' % Dye (2012)
    half_life = [1.28e3 4.97e4 7e6 1.40e4 7.04e2 4.46e3]; % Ma
    abundance = [0.000117 1 1 1 0.007204 0.992796];
    heat_prod = [2.847e-5 0 0 2.628e-5 56.847e-5 9.513e-5];
end

lambda = log(2)./half_life;

if length(K) > 1 & length(age) > 1
    age = age(:)';
end

element = [K(:), Rb(:), Sm(:), Th(:), U(:), U(:)];

density = density(:);

% correct abundances if necessary
if aflag
    for i = 1:3
         x = abundance(i)*exp(t0*lambda(i));
         abundance(i) = x/(x + (1 - abundance(i)));
    end
    
    x = abundance(5:6).*exp(t0*lambda(5:6));
    %x(3) = 0.004446*exp(t0*log(2)/2.455e5); % u234
    abundance(5:6) = x(1:2)/sum(x);
end

for i = 1:size((element),2)
    % compute temporal term
    T = exp(age*lambda(i));              
    
    % determine heat production contribution of each isotope
    for j = 1:length(age)
        A_iso(:,j,i) = density.*heat_prod(i).*element(:,i).*abundance(i)*T(:,j);
        C_iso(:,j,i) = element(:,i).*abundance(i)*T(:,j);
    end
end

% total heat production
A = sum(A_iso,3);

% element concentration
for i = 1:size((element),2)
    if i == 5
        T1 = exp(age*lambda(5));
        T2 = exp(age*lambda(6));
        C_el(:,i) = element(:,i)*((1 - abundance(5) + abundance(6)) + ...
            (abundance(5)*T1 + abundance(6)*T2));
    elseif i == 6
        break;
    else
        T = exp(age*lambda(i));
        C_el(:,i) = element(:,i)*((1 - abundance(i)) + ...
            abundance(i)*T);
    end
end

return


% density lookup table.  if the file does not exist, it will create it on
% the first run.  This table can of course be customized by the user by
% changing rock names and/or density values.
function density = lookup_density(input_name)

if ischar(input_name)
    input_name = {input_name};
end
input_name = input_name(:);

try
    model = readtable('density_of_common_rocks.csv');
catch
    rock_type = {'alkalic gabbro' 'carbonatite' 'cumulate peridotite' ...
        'diorite' 'foid gabbro' 'foid monzodiorite' 'foid monzosyenite' ...
        'foid syenite' 'gabbroic diorite' 'granite' 'granodiorite' ...
        'intermediate foidolite' 'mafic foidolite' 'mantle peridotite' ...
        'monzodiorite' 'monzogabbro' 'monzonite' ...
        'peridotgabbro' 'quartz monzonite' 'quartzolite' ...
        'subalkalic gabbro' 'syenite' 'ultra-high alkali plutonic' ...
        'ultramafic foidolite' 'alkali picrite' 'alkalic basalt' ...
        'andesite' 'basaltic andesite' 'basaltic trachyandesite' ...
        'boninite' 'carbonatite' 'dacite' 'intermediate' 'komatiite' ...
        'mafic foidite' 'meimechite' 'phonolite' 'phonotephrite' ...
        'picrite' 'picrobasalt' 'rhyolite' 'silexite' 'subalkalic basalt' ...
        'tephriphonolite' 'tephrite' 'trachyandesite' 'trachybasalt' ...
        'trachydacite' 'trachyte' 'ultra-high alkali volcanic' ...
        'ultramafic foidite' 'arkose' 'dolomite' 'iron-rich sand' ...
        'iron-rich shale' 'laterite' 'limestone' 'litharenite' ...
        'oxide' 'quartz arenite' 'quartzite' 'shale' 'subarkose' ...
        'sublitharenite' 'wacke'}';
    density = [2954 2825 3160 2787 3001 2831 2757 2631 2897 2660 ...
        2718 2760 2944 3364 2838 2906 2762 3245 2690 2700 ...
        2969 2664 2598 3081 3095 2957 2805 2899 2844 3020 ...
        2846 2736 2789 3164 2903 3288 2644 2830 3111 3062 ...
        2675 2722 2961 2731 2972 2770 2905 2704 2672 2685 ...
        3049 2641 2782 2733 2756 3313 2726 2691 3194 2854 ...
        2677 2682 2731 2745 2669]';

    model = table(rock_type,density);
    writetable(model,'density_of_common_rocks.csv');
end

ind = zeros(size(input_name));
for i = 1:length(input_name)
    tmp = find(strcmpi(model.rock_type,input_name{i}) == 1);
    if ~isempty(tmp)
        ind(i) = tmp;
    end
end
if any(ind == 0)
    warning(['Did not understand some rock names, please check ', ...
        'standard list by running m-file without options.']);
    fprintf('Assuming unknown rocks are granite.');
    
    ind(ind == 0) = 10;
end

density = model.density(ind);

return
