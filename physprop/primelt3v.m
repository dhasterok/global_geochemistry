function [afm,bm] = primelt3v(data,P,varargin)
% method by Herzberg and Asimow, G^3, 2015, doi:10.1002/2014GC005631
% doesn't work yet

global CATIONS OXIDES NCAT MW R

% ideal gas constant (From CODATA)
R = 8.314462618;

%CATIONS = {'Si','Ti','Al','Cr','FeIII','FeII','Mn','Mg','Ca','Na','K','H'};
CATIONS = {'Si','Ti','Al','Cr','FeIII','FeII','Mn','Mg','Ca','Na','K','Ni','P'};
%OXIDES = {'SiO2','TiO2','Al2O3','Cr2O3','Fe2O3','FeO','MnO','MgO','CaO','Na2O','K2O','NiO','P2O5','H2O'};
OXIDES = {'SiO2','TiO2','Al2O3','Cr2O3','Fe2O3','FeO','MnO','MgO','CaO','Na2O','K2O','NiO','P2O5'};

NCAT = [1 1 2 2 2 1 1 1 1 2 2 1 2];

% molecular weights used by Herzberg and Asimow
%MW = [60.085 79.899 101.962 151.99 159.962 71.846 70.937 40.311 56.079 61.979 94.203 74.71 141.945 18.0152];
%MW = [60.085 79.899 101.962 151.99 159.962 71.846 70.937 40.311 56.079 61.979 94.203 74.71 141.945];

% molecular weights using atomic weights table
MW = molecularwt(OXIDES);

% Corrected errors in Herzberg and Asimow's original Excel code:
%   1. I think they transposed a couple of digits in the molecular weight
%   of Fe2O3 (159.962 rather than 159.692.  It is only a ~0.2% error so it
%   shouldn't affect much at all.

% olivine step size
ol_step = 0.01;

% Source peridotite composition
src.feo = 8.02;
src.mgo = 38.12;

% add H2O column if needed
%if ~any(strcmp(data.Properties.VariableNames,'h2o'))
%    data.h2o = zeros(height(data),1);
%end

% molar Fe2+/(Fe2+ + Fe3+)
data.fe_ratio = repmat(0.9,height(data),1);

% Fe2O3/TiO2 ratio
% What is a reasonable value for the default here? I have no idea.  It is
% only needed if fix_feti == true
data.fe2o3_tio2 = repmat(0.89547,height(data),1);

% Fix Fe2O3/TiO2 ratio
fix_feti = false(height(data),1);

if nargin > 2
    opt = 1;
    while opt + 1 <= nargin - 2
        switch lower(varargin{opt})
            case ''
                fe_ratio = varargin{opt+1};
                opt = opt + 2;
            case 'feti'
                fe2o3_tio2 = varargin{opt+1};
                if length(fe2o3_tio2) == 1
                    data.fe2o3_tio2 = repmat(fe_ratio,height(data),1);
                elseif length(fe2o3_tio2) == height(data)
                    data.fe2o3_tio2 = fe2o3_tio2(:);
                else
                    error('fe2o3_tio2 has incorrect dimensions.');
                end
                fix_feti = true(height(data),1);
                fix_feti(isnan(data.fe2o3_tio2)) = false;
                opt = opt + 2;
            case 'fe_ratio'
                fe_ratio = varargin{opt+1};
                if length(fe_ratio) == 1
                    data.fe_ratio = repmat(fe_ratio,height(data),1);
                elseif length(fe_ratio) == height(data)
                    data.fe_ratio = fe_ratio(:);
                else
                    error('fe_ratio has incorrect dimensions.');
                end
                opt = opt + 2;
            otherwise
                error(['Unknown option (',varargin{1},').']);
        end
    end
end

% compute Fe2O3 and FeO from FeO_tot using a fixed molar Fe2+/(Fe2+ + Fe3+)
% ratio or use a fixed Fe2O3/TiO2 ratio to comppute molar
% Fe2+/(Fe2+ + Fe3+) ratio

% fix_feti = true
if any(strcmp(data.Properties.VariableNames,'feo_tot'))
    if any(strcmp(data.Properties.VariableNames,'feo'))
        data.feo = [];
    end
    if any(strcmp(data.Properties.VariableNames,'fe2o3'))
        data.fe2o3 = [];
    end
    
    data = addvars(data,nan(height(data),1),nan(height(data),1),'NewVariableNames',{'fe2o3','feo'},'Before','feo_tot');
else
    error('All iron should be expressed as FeOT');
end
data(fix_feti,:) = fixed_feti(data(fix_feti,:));

% fix_feti = false
data(~fix_feti,:) = fe_conversion(data(~fix_feti,:),data.fe_ratio(~fix_feti));
data.fe2o3_tio2(~fix_feti) = data.fe2o3(~fix_feti)./data.tio2(~fix_feti);

for i = 1:length(OXIDES)
    if ~any(strcmpi(data.Properties.VariableNames,OXIDES{i}))
        data{:,lower(OXIDES{i})} = zeros(height(data),1);
    end
end

if length(P) == 1
    data.pressure = repmat(P,height(data),1);
elseif length(P) == height(data)
    data.pressure = P(:);
else
    error('pressure has incorrect dimensions.');
end
data.temperature = nan(height(data),1);

% normalize chemistry
data = oxide_normalize(data);

% normalize cations
data = cation_normalize(data,1);

% compute primary melt for accumulated fractional melting (AFM) and batch
% melting (BM)
[afm,bm] = fractionation(data,src,ol_step);

return


% fix Fe2O3/TiO2 ratio
function t = fixed_feti(t)

t.fe2o3 = t.tio2.*t.fe2o3_tio2;
t.fe_ratio = 1 - t.fe2o3./t.feo_tot*2*molecularwt('FeO')/molecularwt('Fe2O3');
t.feo = t.feo_tot.*t.fe_ratio;

ind = t.fe_ratio > 1;
if any(ind)
    warning('Fe ratio out of bounds, setting to 1');
end
t(ind,:) = fe_conversion(t(ind,:),1);
t.fe_ratio(ind) = 1;
t.fe2o3_tio2(ind) = t.fe2o3(ind,:)./t.tio2(ind,:);

ind = t.fe_ratio < 0;
if any(ind)
    warning('Fe ratio out of bounds, setting to 0');
end
t(ind,:) = fe_conversion(t(ind,:),0);
t.fe_ratio(ind) = 0;
t.fe2o3_tio2(ind) = t.fe2o3(ind,:)./t.tio2(ind,:);

return


% Primary melt composition
% -------------------------------------------------------------------------
function t = kd_iterator(t)

global R

[SiO2p,psi] = SiO2pound(t);

T = t.temperature + 273.15;
P = t.pressure;

Fo = repmat(0.9,height(t),1);
KD = nan(height(t),1);
oldKD = repmat(0.3,height(t),1);

maxit = 30;
tol = 1e-4;

ind = true(height(t),1);
for c = 1:maxit
    % equilibrium constant
    KD(ind) = exp( (-6766 ./ (R * T(ind)) - 7.34 / R) + ...
        log(0.036 * SiO2p(ind) - 0.22) + ...
        3000 .* (1 - 2 * Fo(ind)) ./ (R * T(ind)) + ...
        0.035 * (P(ind) - 1) ./ (R * T(ind)) );

    % Forsterite number
    Fo(ind) = 1 ./ (1 + KD(ind) .* t.FeII(ind) ./ t.Mg(ind));
    
    % ending conditions
    ind = abs(KD - oldKD) >= tol;
    if ~any(ind)
        break
    end
    
    oldKD(ind) = KD(ind);
end

t.KD = KD;
t.Fo = Fo;
t.SiO2p = SiO2p;
t.psi = psi;

return


function [SiO2p,psi] = SiO2pound(t)

global OXIDES MW

% initialize mole fraction table
molefrac = t(:,lower(OXIDES));

molefrac{:,lower(OXIDES)} = molefrac{:,lower(OXIDES)} ./ repmat(MW,height(molefrac),1);
molefrac{:,lower(OXIDES)} = 100 * molefrac{:,lower(OXIDES)} ./ ...
    repmat(sum(molefrac{:,lower(OXIDES)},2),1,width(molefrac));

psi = nan(height(t),1);
ind = molefrac.sio2 < 60;
psi(ind) = (0.46 * (100 ./ (100 - molefrac.sio2(ind))) - 0.93) .* ...
    (molefrac.na2o(ind) + molefrac.k2o(ind)) + ...
    (-5.33 * (100 ./ (100 - molefrac.sio2(ind))) + 9.69);
psi(~ind) = (11 - 5.5 * (100 ./ (100 - molefrac.sio2(~ind)))) .* ...
    exp(-0.13 * (molefrac.na2o(~ind) + molefrac.k2o(~ind)));

SiO2p = molefrac.sio2 + psi .* (molefrac.na2o + molefrac.k2o);

% Don't know why this adjustment is included in the spreadsheet supplement 
% to the original paper since H2O is not included in the calculations in 
% any other capacity.
%SiO2p = SiO2p - 0.8 * t.h2o;

return


% melt olivine line
function [afm,bm] = fractionation(t,src,ol_step)

global OXIDES

oxides = lower(OXIDES);

olivine_added = zeros(height(t),1);

% initialize KD, Fo, and SiO2p to NaN's
t.KD = nan(height(t),1);
t.Fo = nan(height(t),1);
t.SiO2p = nan(height(t),1);
t.psi = nan(height(t),1);

% start with observed melt composition
melt_init = t;
melt_init.feo_tot = []; % don't need feo_tot anymore
melt_init = addvars(melt_init,zeros(height(t),1),'NewVariableNames','olivine_added','Before',1);

% compute melt temperature
melt_init = temperature_estimator(melt_init);

% compute KD, Fo and SiO2p for the observed melt composition
melt_init = kd_iterator(melt_init);

% compute olivine composition in equilibrium with observed composition
olivine_init = equilibrium_olivine(melt_init);

% compute melt ternary composition
melt_init = mineralfrac(melt_init);

% compute melt fractions
melt_init = meltfrac(melt_init, src);

% set condition for backward fractionation (up-temperature, +ol_step) or
% forward fractionation (down-temperature, -ol_step)
ol_step = repmat(ol_step,height(t),1); % backward fractionation
ind = melt_init.Fproj - melt_init.Fafm < 0;
ol_step(ind) = -ol_step(ind); % forward fractionation

% initialize loop variables
fields = {'sio2','tio2','al2o3','cr2o3','feo','fe2o3','mno','mgo', ...
    'cao','na2o','k2o','nio','p2o5', ...
    'temperature','KD','Fo','SiO2p','psi', ...
    'ol','an','qz','s_cpx','s_gtlhzhz','domain', ...
    'Fbm','Fafm','Fproj'};
melt_old = melt_init;

ind_afm = true(height(t),1);   % flag indicating AFM solution has been reached
ind_bm = true(height(t),1);    % flag indicating batch solution has been reached

% will stop olivine addition when the difference between the melt
% projections change sign
sign_afm = sign(melt_init.Fproj - melt_init.Fafm);  % sign parameter for AFM solution
sign_bm = sign(melt_init.Fproj - melt_init.Fbm);   % sign parameter for batch solution

maxit = 60;
melt_new = melt_init;
olivine_old = olivine_init;
olivine_new = olivine_init;

afm = melt_init;
afm.olivine_added = nan(height(t),1);
afm.residual_lithology = repmat("",height(t),1);
afm.fe_ratio = nan(height(t),1);
afm.fe2o3_tio2 = nan(height(t),1);
afm.potential_temperature = nan(height(t),1);

bm = melt_init;
bm.olivine_added = nan(height(t),1);
bm.residual_lithology = repmat("",height(t),1);
bm.fe_ratio = nan(height(t),1);
bm.fe2o3_tio2 = nan(height(t),1);
bm.potential_temperature = nan(height(t),1);

for c = 1:maxit
    % logical for samples to continue computing
    %[ind_afm ind_bm]
    ind = ind_afm | ind_bm;
    if ~any(ind)
        break
    end
    
    % set old to previous new
    melt_old(ind,:) = melt_new(ind,:);
    olivine_old(ind,:) = olivine_new(ind,:);
    
    % increment olivine added
    olivine_added(ind) = olivine_added(ind) + ol_step(ind) * 100;
    melt_new.olivine_added(ind,:) = olivine_added(ind,:);
    
    % compute new melt composition
    melt_new{ind,oxides} = (melt_old{ind,oxides} + ...
        olivine_old{ind,oxides} .* repmat(ol_step(ind),1,length(oxides))) ./ ...
        (1 + repmat(ol_step(ind),1,length(oxides)));

    % normalize cations
    melt_new(ind,:) = cation_normalize(melt_new(ind,:),1);

    % compute temperature of melt
    melt_new(ind,:) = temperature_estimator(melt_new(ind,:));

    % compute KD, Fo and SiO2p from new melt composition
    melt_new(ind,:) = kd_iterator(melt_new(ind,:));

    % compute olivine composition in equilibrium with new melt composition
    olivine_new(ind,:) = equilibrium_olivine(melt_new(ind,:));

    % compute melt ternary composition
    melt_new(ind,:) = mineralfrac(melt_new(ind,:));
    
    % compute melt fractions
    melt_new(ind,:) = meltfrac(melt_new(ind,:), src);
    
    %[sign_afm sign(melt_new.Fproj - melt_new.Fafm) sign_bm sign(melt_new.Fproj - melt_new.Fbm)]
    %if phi_afm_new >= phi_afm_old & ~afm_flag
    ind2 = sign_afm ~= sign(melt_new.Fproj - melt_new.Fafm) & ind_afm;
    
    % initialize afm table
    afm(ind2,melt_old.Properties.VariableNames) = melt_old(ind2,:);

    % project olivine added to Fafm - Fproj = 0
    afm.olivine_added(ind2) = -(melt_new.olivine_added(ind2) - melt_old.olivine_added(ind2)) .* ...
        (melt_new.Fafm(ind2) - melt_new.Fproj(ind2)) ./ ...
        ( (melt_new.Fafm(ind2) - melt_new.Fproj(ind2)) - ...
        (melt_old.Fafm(ind2) - melt_old.Fproj(ind2)) ) + ...
        melt_new.olivine_added(ind2);

    % average other fields
    afm{ind2,fields} = avg_fields(melt_old{ind2,fields}, melt_new{ind2,fields}, ...
        repmat(melt_old.olivine_added(ind2),1,length(fields)), ...
        repmat(melt_new.olivine_added(ind2),1,length(fields)), ...
        repmat(afm.olivine_added(ind2),1,length(fields)));

    % set flag (AFM iterations are done!)
    ind_afm(ind2) = false;

    
    %if phi_bm_new >= phi_bm_old & ~bm_flag
    ind3 = sign_bm ~= sign(melt_new.Fproj - melt_new.Fbm) & ind_bm;
    % initialize bm table
    bm(ind3,melt_old.Properties.VariableNames) = melt_old(ind3,:);

    % project olivine added to Fbm - Fproj = 0
    bm.olivine_added(ind3) = -(melt_new.olivine_added(ind3) - melt_old.olivine_added(ind3)) .* ...
        (melt_new.Fbm(ind3) - melt_new.Fproj(ind3)) ./ ...
        ( (melt_new.Fbm(ind3) - melt_new.Fproj(ind3)) - ...
        (melt_old.Fbm(ind3) - melt_old.Fproj(ind3)) ) + ...
        melt_new.olivine_added(ind3);

    % average other fields
    bm{ind3,fields} = avg_fields(melt_old{ind3,fields}, melt_new{ind3,fields}, ...
        repmat(melt_old.olivine_added(ind3),1,length(fields)), ...
        repmat(melt_new.olivine_added(ind3),1,length(fields)), ...
        repmat(bm.olivine_added(ind3),1,length(fields)));

    % set flag (BM iterations are done!)
    ind_bm(ind3) = false;
    
    % if prim_afm == 0
    %     warning('No AFM solution found.');
    % end
end

% lookup lithology (based on domains)
% set all to domain 2 first, then correct
afm.residual_lithology = repmat("spinel peridotite",height(t),1);
bm.residual_lithology = repmat("spinel peridotite",height(t),1);
% domain 1
afm.residual_lithology(afm.domain == 1) = "garnet peridotite";
bm.residual_lithology(bm.domain == 1) = "garnet peridotite";
% domain 3
afm.residual_lithology(afm.domain == 3) = "harzburgite";
bm.residual_lithology(bm.domain == 3) = "harzburgite";

% correct cation ratios
bm = cation_normalize(bm);
afm = cation_normalize(afm);

% correct Fe2+/(Fe2+ + Fe3+) ratio 
bm.fe_ratio = bm.feo/molecularwt('FeO') ./ ...
    (bm.feo/molecularwt('FeO') + 2*bm.fe2o3/molecularwt('Fe2O3'));
afm.fe_ratio = afm.feo/molecularwt('FeO') ./ ...
    (afm.feo/molecularwt('FeO') + 2*afm.fe2o3/molecularwt('Fe2O3'));

% correct Fe2O3/TiO2 ratio
bm.fe2o3_tio2 = bm.fe2o3./bm.tio2;
afm.fe2o3_tio2 = afm.fe2o3./afm.tio2;

% pressure function
p_fcn = pressurefunction;

% compute potential temperatures and pressure range of melting
afm = potential_temperature(afm);
[afm.pressure_initial,afm.pressure_final] = estimate_pressures(afm,p_fcn);

bm = potential_temperature(bm);
[bm.pressure_initial,bm.pressure_final] = estimate_pressures(bm,p_fcn);

%t = t(ip)*(prim_afm - ol_percent(ip-1)) + t(ip-1)*(ol_percent(ip) - prim_afm) / (ol_percent(ip) - ol_percent(ip-1));
% if abs(t.Fafm - proj_f1) < 0.0000001
%     fractional.residual = 'Harzburgite';
% elseif ABS(Fafm - proj_f2) < 0.000001
%     fractional.residual = 'Spinel Peridotite';
% else
%     fractional.residual = 'Garnet Peridotite';
% end

afm.error1(t.cao < 13.81 - 0.274 * t.mgo) = "Pyroxenite Source, no solution";
afm.error1(t.cao > 2.318 * t.sio2 - 93.626) = "Volatile Peridotite Source, no solution";

afm.error2(afm.mgo <= 20.6 & afm.cao < 11.436 - 0.104*afm.mgo | ...
    afm.mgo > 20.6 & afm.cao < -23.209 + 3.643*afm.mgo - 0.1*afm.mgo.^2 ) = ...
    "Augite fractionation warning";
afm.error2(afm.cao > 1.095 + 0.154*afm.mgo + 116.58./afm.mgo) = "Augite accumulation warning";

afm.error3(afm.feo < 0.577*afm.mgo - 0.00026*afm.mgo.^3 + 5.48./afm.mgo) = "FeO/MgO forbidden, no solution";
afm.error4(afm.qz < 0) = "Negative normative Quartz";

return


% Computes olivine composition in equilibrium with the melt
function olivine = equilibrium_olivine(t)

global CATIONS OXIDES NCAT MW

% Partition coefficient constants
c = [0       0.0035 0.067 0      0.263 0.214 1  0.0071 0      0       3.346 0;
     0.0300 -0.031  0.183 0.0001 0.196 0.118 0 -0.019  0.0001 0.0001 -3.665 0.0001;
     0       0.093  0     0      0     0     0  0.063  0      0       0     0];

oxides = lower(OXIDES);

% Olivine-melt partition coefficients
D = t(:,oxides);   % initialize D
D{:,oxides} = NaN;

D{:,oxides(2:end)} = 100 * (2/3 * repmat(c(1,:),height(t),1) .* t.Fo ./ t.Mg + ...
    repmat(c(2,:),height(t),1) + ...
    repmat(c(3,:),height(t),1) .* t.Mg ./ (2 * t.Fo / 3)) .* t{:,CATIONS(2:end)};

sum1 = 0.5*(D.al2o3 - D.cr2o3 - D.na2o - D.k2o - D.p2o5) + ...
    D.cr2o3 + D.fe2o3 + D.feo + D.mno + D.mgo + D.cao + ...
    2*(D.na2o + D.k2o + D.p2o5) + D.nio;

% cation fraction olivine?
olivine = t(:,oxides);   % initialize olivine_eq
olivine{:,oxides} = NaN;

easy = {'tio2','al2o3','cr2o3','fe2o3','mno','cao','na2o','k2o','p2o5'};
olivine{:,easy}  = D{:,easy} * 200/3 ./ sum1;

sum2 = 0.5*(olivine.al2o3 - olivine.cr2o3 - olivine.na2o - ...
    olivine.k2o - olivine.p2o5) + olivine.cr2o3 + ...
    olivine.fe2o3 + olivine.mno + olivine.cao + ...
    2*(olivine.na2o + olivine.k2o + olivine.p2o5) + D.nio;

olivine.sio2 = 100/3 - olivine.tio2 - 0.5*(olivine.al2o3 - ...
    olivine.cr2o3 - olivine.na2o - olivine.k2o - ...
    olivine.p2o5) - olivine.cr2o3;

olivine.mgo = t.Fo .* (200/3 - sum2);
olivine.feo = (1 - t.Fo) .* (200/3 - sum2);
olivine.nio = (3.346 * (olivine.mgo ./ t.Mg) / 100 - 3.665) .* t.Ni * 100;

% compute mass fraction
nr = height(t);
olivine{:,:} = 100 * olivine{:,:} .* repmat(MW,nr,1) ./ repmat(NCAT,nr,1);
olivine{:,:} = 100 * olivine{:,:} ./ sum(olivine{:,:},2);

return


% Estimate the melt temperature at a given pressure
function t = temperature_estimator(t)

global CATIONS R

% WTF are these values? The come out of nowhere in the spreadsheet.
c = [0    0      0      0      0      0.279  0.259 1 0.0056 0      0       3.346 0;
     0.82 0.0001 0.0001 0.0001 0.0001 0.031 -0.049 0 0.0135 0.0001 0.0001 -3.665 0.0001];

sumA = sum(t{:,CATIONS(2:end)}.*c(1,2:end),2);
sumB = sum(t{:,CATIONS(2:end)}.*c(2,2:end),2);
dmgo2 = (2/3 - sumB) ./ sumA;

% compute temperature
t.temperature = (113100/R + t.Si * 0.00000411 / R) ./ (52.05 / R + 2 * log(dmgo2) + ...
    2 * log(1.5 * (t.FeII + t.Mn + t.Mg + t.Ca + t.Ni)) + ...
    2*log(3 * t.Si) - (3.5 * log(1 - t.Al) + 7 * log(1 - t.Ti))) - 273.15;

return


% Composition functions
% -------------------------------------------------------------------------
function t = meltfrac(t, src)

tol = 1e-4;

feo = t.feo;
mgo = t.mgo;

KD2 = 0.3813 - 0.7896 ./ mgo + 1.0389 ./ mgo.^2;

% initialize table variable
% not set to quite zero to handle future numerical instabilities
%t.Fbm = repmat(1e-5,height(t),1);

% initialize loop variables
fb = (mgo .* (src.feo ./ feo) - KD2 .* src.mgo) ./ ...
    (mgo .* (src.feo ./ feo) - mgo .* KD2);
fa = zeros(size(fb));
fc = zeros(size(fb));

ind = feo ~= 0;
fa(~ind) = 1e-5;

maxit = 50;
for c = 1:maxit
    fa(ind) = ( mgo(ind) .* (( (src.feo - fb(ind) .* feo(ind)) ./ ...
        (1 - fb(ind)) ) ./ feo(ind)) - KD2(ind) * src.mgo ) ...
        ./ ( mgo(ind) .* (((src.feo - fb(ind) .* feo(ind)) ./ ...
        (1 - fb(ind))) ./ feo(ind)) - mgo(ind) .* KD2(ind) );

    ind2 = fc ~= 0 & (fa - fb) ./ (fc - fb) > 0.9;
    fa(ind & ind2) = (fa(ind & ind2) + fb(ind & ind2)) / 2;
    
    fc(ind) = fb(ind);
    fb(ind) = fa(ind);
    
    % keep going?
    ind = abs(fc - fb) >= tol;

    if ~any(ind)
        break
    end
end
t.Fbm = fa;

ind = t.Fbm >= 0;
t.Fafm(ind) = 1./(0.98 ./ t.Fbm(ind) - 0.9 * log(t.Fbm(ind)).^2 + ...
    0.07 * log(t.Fbm(ind)).^4);
t.Fafm(~ind) = t.Fbm(~ind);

% F1
ind = t.domain == 1;
t.Fproj(ind) = 6.2819 * t.an(ind).^2 - 14.7789 * t.an(ind).^3 + ...
    0.00825 * (1 ./ t.an(ind)).^2;

% F2
ind = t.domain == 2;
t.Fproj(ind) = ((-1.994 + 2.25 * t.qz(ind) + 0.041 ./ t.qz(ind)) + ...
    (-1.183 - 3.005 * t.qz(ind) + 13.774 * t.qz(ind).^2 - 12.615 * t.qz(ind).^3)) / 2 + ...
    exp(0.931 + 1.623.*t.qz(ind)) .* t.qz(ind).^0.245 .* t.ol(ind) + ...
    exp(0.769 - 7.514 * t.qz(ind)) .* t.qz(ind).^0.577 ./ t.ol(ind);
t.Fproj(ind & t.qz <= 0) = 0;

ind2 = ind & t.Fproj > 0 & t.qz < -0.1773 + 0.1404 ./ t.ol - 0.008434 ./ t.ol.^2;
t.Fproj(ind2) = -t.Fproj(ind2);

% F3
ind = t.domain == 3;
t.Fproj(ind) = -2.5345 + 5.329 * (t.qz(ind) + t.ol(ind)*0.348) + 0.3012 ./ ...
    (t.qz(ind) + 0.348 * t.ol(ind));

ind2 = ind & t.Fproj > 0 & t.qz < -0.0896 + 0.02002 ./ t.ol + 0.02989 ./ t.ol.^2;
t.Fproj(ind2) = -t.Fproj(ind2);

return


function t = mineralfrac(t)

global OXIDES MW

molefrac = t(:,lower(OXIDES));

molefrac{:,lower(OXIDES)} = molefrac{:,lower(OXIDES)} ./ repmat(MW,height(molefrac),1);
molefrac{:,lower(OXIDES)} = 100 * molefrac{:,lower(OXIDES)} ./ ...
    repmat(sum(molefrac{:,lower(OXIDES)},2),1,width(molefrac));
% olivine
ol = -1.5*molefrac.tio2 + 0.5*(molefrac.al2o3 + molefrac.cr2o3 + ...
    molefrac.feo + molefrac.mno + molefrac.mgo + molefrac.nio) - ...
    0.5*molefrac.cao - 0.5*molefrac.na2o + 3*molefrac.k2o;
% anorthosite
an = molefrac.tio2 + molefrac.al2o3 + molefrac.cr2o3;
% quartz
qz = molefrac.sio2 - 0.5*(molefrac.al2o3 + molefrac.cr2o3 + ...
    molefrac.feo + molefrac.mno + molefrac.mgo + molefrac.nio) - ...
    1.5*molefrac.cao - 3*(molefrac.na2o + molefrac.k2o);

% normalize
sum_min = (ol + an + qz);
t.ol = ol./sum_min;
t.an = an./sum_min;
t.qz = qz./sum_min;

% clinopyroxene
t.s_cpx = -0.074 + 0.1713./t.ol - 0.0135./t.ol.^2;

t.s_gtlhzhz = zeros(height(t),1);
ind = t.an < 700;
t.s_gtlhzhz(ind) = 1./(16.843 + 28.733*t.an(ind) - 14.183*exp(t.an(ind)));

% determine compositional domain
t.domain = repmat(2,height(t),1);
t.domain(t.ol > 0.5 & t.qz < t.s_gtlhzhz) = 3;
t.domain(t.domain ~= 3 & t.qz > t.s_cpx) = 1;

return


% mantle potential temperature
function t = potential_temperature(t)

t.potential_temperature = 1025 + 28.6*t.mgo - 0.084*t.mgo.^2;

return


% averaging function for table fields all at once
function t_avg = avg_fields(t1,t2,ol1,ol2,ol0)

% linear interpolation
t_avg = (t2 - t1) ./ (ol2 - ol1) .* (ol0 - ol1) + t1;

return


function p_fcn = pressurefunction
% following the procedure of Hole and Millett (J. Petrology 2016)
% doi:10.1093/petrology/egw014

mgo = [12:0.01:24];

feo = [-0.0132*mgo.^2 + 0.7500*mgo + 0.094;                     % 1 GPa
    -0.0158*mgo.^2 + 0.8234*mgo + 0.0326;                       % 2 GPa
    -0.0142*mgo.^2 + 0.7410*mgo + 1.2802;                       % 3 GPa
    -0.00663379*mgo.^3 + 0.444349*mgo.^2 - 9.7573*mgo + 81.082246];   % 4 GPa
% in the supplment the equation for the 4 GPa is given as, but it
% yields values of feo that are too high for a given mgo.
% feo = -0.0069*mgo.^3 + 0.4643*mgo.^2 - 10.182*mgo + 84.089;

p = repmat([1:4]',1,length(mgo));

mgo = repmat(mgo,4,1);

feo = feo(:);
mgo = mgo(:);
p = p(:);

ind = ~(p == 4 & mgo < 18);
feo = feo(ind);
mgo = mgo(ind);
p = p(ind);

p_fcn = scatteredInterpolant(feo,mgo,p);

return


function [pi,pf] = estimate_pressures(t,p_fcn)
% following the procedure of Hole and Millett (J. Petrology 2016)
% doi:10.1093/petrology/egw014

% compute initial pressure
pi = 11.248*t.mgo - 13700*(1./t.mgo).^3 - 8.13*(log(t.mgo)).^3;

% estimate final pressure
pf = nan(height(t),1);

ind = t.Fafm >= 0.1;
pf(ind) = p_fcn(t.feo(ind),t.mgo(ind));

ind = t.Fafm < 0.1;
pf(ind) = pi(ind) - 4.659*t.Fafm(ind) + 10.240*t.Fafm(ind).^2;

return
