function prim = basaltPT(data,varargin)
% BASALTPT - estimates PT from basalt composition.
%
%   prim = basaltPT(data) estimates the P and T of a mafic rock 43 < SiO2 <
%   60 wt.%.  The function returns the P and T conditions prim.pressure and
%   prim.temperature along with the estimated primative composition.  The
%   table prim also contains the molar percentage of cations and the the
%   oxygen 8 normalized compositions (see original paper for details).
%
%   There are several options that can be specified as option value pairs,
%   i.e., prim = basaltPT(data,option,value)
%   Options:
%       'MgO_min'           Minimum allowable MgO for calculation.  The
%                           default value is 10 wt%.
%
%       'Fe_ratio'          Molar Fe2+/(Fe2+ + Fe3+) ratio.  The default
%                           value is 0.95.  Note in the original paper the
%                           value was specified as Fe3+/(Fe2+ + Fe3+)
%
%       'Fo'                Target forsterite (Fo) number.  The default
%                           value is 0.9.
%
%       'KD'                The equilibrium constant is the ratio of the
%                           activity of opx to the product of ol and SiO2
%                           in the melt. The default value is 0.32, but
%                           values typically range (0.3,0.35) between 10
%                           and 30 wt% MgO.  If you specify a value here,
%                           you will also need to specify variable_KD = 0,
%                           otherwise it will be ignored.
%
%       'variable_KD'       Use a fixed (0) or variable (1) KD value.  The
%                           default is variable_kd = 1
%
% Original method by Lee et al., EPSL, 2009, doi:10.1016/j.epsl.2008.12.020
% Adapted by D. Hasterok (Univ. Adelaide) with fixes

addpath toolbox

% Corrected errors in Lee's original Excel code:
%   1. Cr2O3 was inadvertently set to 0, this affects the normalization,
%   but since Cr2O3 is small it shouldn't affect the result much.
%   2. In the computation of cation percentages from the original data the
%   fe2o3 percentage use a value from the normalized rather than
%   unnormalized data...this mistake isn't propagated through any
%   subsequent calculations, but will give the wrong cation numbers for the
%   original melt.
%   3. The mass of Fe2O3 and FeO is incorrect/inconsistent in a couple of
%   places in the original spreadsheet.  The way I've formulated it, the
%   molcular weights are only specified once in the specified once to
%   prevent such errors.
%   4. In the iterative primitive melt calculation, KD was fixed rather
%   than updated iteratively (I've contacted Lee to see if this was
%   unintended, but have not heard back.  The text makes it appear that it
%   was incorrect).  The result of the change produces less scatter in the
%   output of the original supplied data and slightly higher temperatures
%   on average.

% default values
% KD = (Fe/Mg)ol/(Fe/Mg)melt
variable_kd = true;
KD = 0.32;

% Mg number of mantle source
Fo_max = 0.9;

% Minimum MgO wt%  of evolved magma
mgo_min = 10;

% molar Fe2+/(Fe2+ + Fe3+)
% Lee defines this as Fe3+/(Fe2+ + Fe3+), but I've converted here because
% this seems to be the more common way to specify the ratio and makes it
% compatible with the function I wrote that is used by other functions
% similar to this one.
fe_ratio = 1 - 0.05;

% step size for observed to primative melt adjustment
step_size = 0.005;
if nargin > 1
    opt = 1;
    while opt + 1 <= nargin
        switch lower(varargin{opt})
            case 'variable_kd'
                variable_kd = varargin{opt+1};
                opt = opt + 2;
            case 'kd'
                KD = varargin{opt+1};
                opt = opt + 2;
            case 'fo'
                Fo_max = varargin{opt+1};
                opt = opt + 2;
            case 'mgo_min'
                mgo_min = varargin{opt+1};
                opt = opt + 2;
            case 'fe_ratio'
                fe_ratio = varargin{opt+1};
                opt = opt + 2;
            otherwise
                error(['Unknown option (',varargin{1},').']);
        end
    end
end

global CATIONS OXIDES SPECIES NCAT MW

CATIONS = {'Si','Ti','Al','Cr','FeIII','FeII','Mn','Mg','Ca','Na','K','H'};
OXIDES = {'SiO2','TiO2','Al2O3','Cr2O3','Fe2O3','FeO','MnO','MgO','CaO','Na2O','K2O','H2O'};
SPECIES = {'Si4O8','Ti4O8','Al16_3O8','Cr16_3O8','Fe16_3O8','Fe4Si2O8', ...
    'Mn4Si2O8','Mg4Si2O8','Ca4Si2O8','Na2Al2Si2O8','K2Al2Si2O8','H16O8'};
NCAT = [1 1 2 2 2 1 1 1 1 2 2 2];

% molecular weights used by Lee
%MW = [60.08 79.86 101.96 151.99 159.69 71.84 70.94 40.3 56.08 61.98 94.2 18.02];
% molecular weights using atomic weights table
MW = molecularwt(OXIDES);

% add H2O column if needed
if ~any(strcmp(data.Properties.VariableNames,'h2o'))
    data.h2o = zeros(height(data),1);
end

for i = 1:length(OXIDES)
    switch OXIDES{i}
        case 'Fe2O3'
            continue
        case 'FeO'
            ind = isnan(data{:,'feo_tot'});
            data{ind,'feo_tot'} = 0;
        otherwise
            ind = isnan(data{:,lower(OXIDES{i})});
            data{ind,lower(OXIDES{i})} = 0;
    end
end

% estimate Fe3+
data = fe_conversion(data,fe_ratio);
%sum(data{1:10,lower(OXIDES)},2)

% normalize oxides to 100%
data = oxide_normalize(data);
%sum(data{1:10,lower(OXIDES)},2)
%data(1:10,lower(OXIDES))

% compute cations % (anhydrous)
data = cation_normalize(data);
%sum(data{1:10,CATIONS},2)
%data(1:10,CATIONS) % different from Lee

[data.Fo,data.KD] = compute_fo(data.Mg,data.FeII,variable_kd,KD);
%data.Fo(1:10) % same as Lee

%data(1:10,:)

% compute primative composition iteratively
tmp = lower(OXIDES);
%prim{:,2:end} = nan([height(prim),length(OXIDES)+2+length(CATIONS)]);

ind = 43 > data.sio2 | data.sio2 > 60 | data.mgo < mgo_min;
prim = obs_to_prim( ...
    data(:,{'sample_id',tmp{1:end},'Fo','KD',CATIONS{1:end}}), ...
    step_size, Fo_max, variable_kd, KD);
prim{ind,{tmp{1:end},'Fo','KD','ol','iter',CATIONS{1:end}}} = NaN;
disp([num2str(sum(ind)),' samples were rejected']);

%prim = obs_to_prim(prim,step_size,Fo_max,variable_kd,KD);
    
%prim(1:10,:)
%sum(prim{1:10,lower(OXIDES)},2)

prim = cation_normalize(prim);
prim = oxygen8_normalization(prim);
%prim(1:10,:)
prim.temperature = 916.45 + 13.68 * prim.Mg4Si2O8 + 4580./prim.Si4O8 - ...
    0.509 * prim.H16O8 .* prim.Mg4Si2O8;

%P = (log(species.Si) - 4.019 + 0.0165 * species.FeII + 0.0005 * species.Ca.^2) ./ ...
%    (-770*T.^-1 + 0.0058*sqrt(T) - 0.003*species.H);

%PLee parameters
b0 = 4.019;     %basic constant
b1 = -770;      %temperature 1/T
b2 = 0.0058;    %temperature T^0.5
b3 = 0.0165;	%Fe
b4 = 0.0005;	%Ca
b5 = 0.003;     %water

prim.pressure = (log(prim.Si4O8) - b0 + b3 * prim.Fe4Si2O8 + b4 * prim.Ca4Si2O8.^ 2) ./ ...
    (b1 ./ (prim.temperature + 273.15) + b2 * (prim.temperature + 273.15).^0.5 - b5 * prim.H16O8);

return


% Computes Forsterite fraction
function [Fo,KD] = compute_fo(Mg,Fe,variable_kd,KD)

% determine kd value
if variable_kd == 1
    % Tamura et al. (J. Petrol., 2000)
    KD = 0.25324 + 0.0033663 * (Mg + 0.33 * Fe);
end

% Compute forsterite fraction
Fo = 1 ./ (1 + KD .* Fe ./ Mg );

return


% Estimates primitive melt composition
function t = obs_to_prim(t,step_size,Fo_max,variable_kd,KD)

global OXIDES

tol = 0.0005;
maxit = 200;

dm = repmat(step_size,height(t),1);
ind = t.Fo < Fo_max;
dm(~ind) = -step_size;

t = addvars(t, ones(height(t),1), 'After','KD','NewVariableNames','ol');
t = addvars(t, zeros(height(t),1), 'After','ol','NewVariableNames','iter');
for c = 1:maxit
    ind = abs(t.Fo - Fo_max) >= tol;
    if ~any(ind)
        break
    end
    
    %disp([c t.KD t.Fo t.ol t{1,lower(OXIDES)}]);
    MWol = 2 * (t.Fo(ind) * 40.3044 + (1 - t.Fo(ind)) * 71.8444) + 60.0843;
    FeOol = 100 * (2 * (1 - t.Fo(ind)) * 71.8444) ./ MWol;
    MgOol = 100 * 2 * t.Fo(ind) * 40.3044 ./ MWol;
    SiO2ol = 100 * 60.0843 ./ MWol;

    t.sio2(ind) = (SiO2ol .* dm(ind) + t.sio2(ind)) ./ (1 + dm(ind));
    t.feo(ind) = (FeOol .* dm(ind) + t.feo(ind)) ./ (1 + dm(ind));
    t.mgo(ind) = (MgOol .* dm(ind) + t.mgo(ind)) ./ (1 + dm(ind));
    
    t{ind,{'tio2','al2o3','cr2o3','fe2o3','mno','cao','na2o','k2o','h2o'}} = ...
        t{ind,{'tio2','al2o3','cr2o3','fe2o3','mno','cao','na2o','k2o','h2o'}} ./ (1 + dm(ind));
    
    % if KD changes each iteration
    t(ind,:) = cation_normalize(t(ind,:));
    [t.Fo(ind),t.KD(ind)] = compute_fo(t.Mg(ind),t.FeII(ind),variable_kd,KD);
    
    % if KD is constant during the iterative process
    %[t.Fo(ind),t.KD(ind)] = compute_fo(t.mgo(ind)/molecularwt('MgO'),t.feo(ind)/molecularwt('FeO'),0,t.KD(ind));
    
    t.iter(ind) = c;
    t.ol(ind) = (1 + dm(ind)).^c;
end
%disp([c t.KD t.Fo t.ol t{1,lower(OXIDES)}]);

return


function t = oxygen8_normalization(t)

global SPECIES

t.Si4O8 = 0.25 * (t.Si - 0.5 * (t.FeII + t.Mg + t.Ca + t.Mn) - 0.5*(t.Na + t.K));
t.Ti4O8 = 0.25 * t.Ti;
t.Al16_3O8 = 0.1875 * (t.Al - t.Na);
t.Cr16_3O8 = 0.1875 * t.Cr;
t.Fe16_3O8 = 0.1875 * t.FeIII;
t.Fe4Si2O8 = 0.25 * t.FeII;
t.Mn4Si2O8 = 0.25 * t.Mn;
t.Mg4Si2O8 = 0.25 * t.Mg;
t.Ca4Si2O8 = 0.25 * t.Ca;
t.Na2Al2Si2O8 = 0.5 * t.Na;
t.K2Al2Si2O8 = 0.5 * t.K;
t.H16O8 = 0.625 * t.H;

t{:,SPECIES} = 100 * t{:,SPECIES} ./ repmat(sum(t{:,SPECIES},2),1,length(SPECIES));

return

% % CALCULATION OF PRESSURES AND TEMPERATURES
% 
% % PLee parameters
% b0 = 4.019;     %basic constant
% b1 = -770;      %temperature 1/T
% b2 = 0.0058;    %temperature T^0.5
% b3 = 0.0165;    %Fe
% b4 = 0.0005;    %Ca
% b5 = 0.003;     %water
% 
% % Putirka 2005 with Na, K, H2O
% % T is in Celsius, P is GPa
% Tp = 3063.2 / (Log(Fo * 100 * 2 / 3 / Mg) + 2.106193 ...
%     - 0.019 * tmp.sio2 - 0.08 * (tmp.na2o + tmp.k2o) + 0.028 * tmp.h2o);
% Pa = 0.1 * (Exp(0.00252 * Tp - 0.12 * tmp.sio2 + 5.027));
% PLee = (log(species.Si) - b0 + b3 * species.FeII + b4 * species.Ca.^ 2) ./ ...
%     (b1 ./ (Tp + 273.15) + b2 * (Tp + 273.15).^ 0.5 - b5 * species.H);
% 
% 
% % Putirka 2005 without compositional effects
% % T is in Celsius, P is GPa
% Tpx = 4490.5 / (log(Fo * 100 * 2 / 3 / Mg) + 2.02);
% Pax = 0.1 * Exp(0.00252 * Tpx - 0.12 * tmp.sio2 + 5.027);
% PLeex = (log(species.Si) - b0 + b3 * species.FeII + b4 * species.Ca.^ 2) ...
%     ./ (b1 ./ (Tpx + 273.15) + b2 * (Tpx + 273.15).^ 0.5 - b5 * species.H);
% 
% 
% % Lee T-independent barometer
% c0 = 4.05
% c1 = -1.49
% c2 = 0.042
% c3 = 0.012
% 
% PLeez = (Log(species.Si) - c0 + c3 * species.FeII) / (c1 / (species.Mg ^ (1 / 3)) + c2 * species.Mg ^ (1 / 2));
% 
% 
% % Sugawara, Lee, Albarede
% if K = 0
%      ActiveCell.Offset(x, 83) = 0;
%  else
%      ActiveCell.Offset(x, 95).GoalSeek Goal:=0, ChangingCell:=ActiveCell.Offset(x, 83)
%  end
%  Tsug1 = ActiveCell.Offset(x, 83)
% 
%  if Tsug1 < 0 | Tsug1 > 2000
%      warning('');
%  end
% 
%  if K = 0
%      PLeeSug = 0;
%  else
%      PLeeSug = (log(species.Si) - b0 + b3 * species.FeII + b4 * species.Ca.^ 2) ./ (b1 ./ (Tsug1 + 273.15) + b2 * (Tsug1 + 273.15).^ 0.5 - b5 * species.H);
%  end
%  if K = 0
%      ActiveCell.Offset(x, 84) = 0;
%  else
%      ActiveCell.Offset(x, 84) = PLeeSug;
%  end
% 
%  if K = 0
%      ActiveCell.Offset(x, 86) = 0;
%  else
%      ActiveCell.Offset(x, 96).GoalSeek Goal:=0, ChangingCell:=ActiveCell.Offset(x, 86);
%  end
%  Tsug2 = ActiveCell.Offset(x, 86)
% 
% 
%  if K = 0
%      PAlbSug = 0;
%  else
%      PAlbSug = 0.1 * Exp(0.00252 * Tsug2 - 0.12 * tmp.sio2 + 5.027);
%  end
%  if K = 0;
%      ActiveCell.Offset(x, 87) = 0;
%  else
%      ActiveCell.Offset(x, 87) = PAlbSug;
%  end