function [afm,bm] = primelts3(data,P,varargin)
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
if ~any(strcmp(data.Properties.VariableNames,'feo'))
    data = addvars(data,nan(height(data),1),nan(height(data),1),'NewVariableNames',{'fe2o3','feo'},'Before','feo_tot');
end
data(fix_feti,:) = fixed_feti(data(fix_feti,:));

% fix_feti = false
data(~fix_feti,:) = fe_conversion(data(~fix_feti,:),data.fe_ratio(~fix_feti));
data.fe2o3_tio2(~fix_feti) = data.fe2o3(~fix_feti)./data.tio2(~fix_feti);


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
for i = 1:height(data)
    [afm(i,:),bm(i,:)] = fractionation(data(i,:),src,ol_step);
end

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

R = 8.3143;
[SiO2p,psi] = SiO2pound(t);

T = t.temperature + 273.15;

if T == 0
    KD = 0.3;
else
    Fo = 0.9;
    oldKD = 0.3;
    iter = 1;
    maxit = 30;
    tol = 1e-4;
    while 1
        % equilibrium constant
        KD = exp( (-6766 ./ (R * T) - 7.34 / R) + ...
            log(0.036 * SiO2p - 0.22) + ...
            3000 .* (1 - 2 * Fo) ./ (R * T) + ...
            0.035 * (t.pressure - 1) ./ (R * T) );
        
        % Forsterite number
        Fo = 1 ./ (1 + KD .* t.FeII ./ t.Mg);

        % ending conditions
        if abs(KD - oldKD) < tol
            break
        end
        oldKD = KD;
        
        iter = iter + 1;
        
        if iter > maxit
            %warning('kd_iterator: maximum number of iterations exceeded.');
            break
        end
    end
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

molefrac{:,lower(OXIDES)} = t{:,lower(OXIDES)} ./ MW;
molefrac{:,lower(OXIDES)} = 100 * molefrac{:,lower(OXIDES)} ./ ...
    repmat(sum(molefrac{:,lower(OXIDES)},2),1,length(MW));

if molefrac.sio2 < 60
    psi = (0.46 * (100 / (100 - molefrac.sio2)) - 0.93) .* ...
        (molefrac.na2o + molefrac.k2o) + ...
        (-5.33 * (100 ./ (100 - molefrac.sio2)) + 9.69);
else
    psi = (11 - 5.5 * (100 ./ (100 - molefrac.sio2))) .* ...
        exp(-0.13 * (molefrac.na2o + molefrac.k2o));
end

SiO2p = molefrac.sio2 + psi * (molefrac.na2o + molefrac.k2o);

% Don't know why this adjustment is included in the code.  H2O is not included in the
% calculations in any other capacity.
%SiO2p = SiO2p - 0.8 * t.h2o;

return


% melt olivine line
function [afm,bm] = fractionation(t,src,ol_step)

global OXIDES

oxides = lower(OXIDES);

olivine_added = 0;

% initialize KD, Fo, and SiO2p to NaN's
t.KD = nan(height(t),1);
t.Fo = nan(height(t),1);
t.SiO2p = nan(height(t),1);
t.psi = nan(height(t),1);

% start with observed melt composition
melt_init = t;
melt_init.feo_tot = []; % don't need feo_tot anymore
melt_init = addvars(melt_init,0,'NewVariableNames','olivine_added','Before',1);

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
if melt_init.Fproj - melt_init.Fafm < 0
    ol_step = -ol_step; % forward fractionation
end

% initialize loop variables
fields = {'sio2','tio2','al2o3','cr2o3','feo','fe2o3','mno','mgo', ...
    'cao','na2o','k2o','nio','p2o5', ...
    'temperature','KD','Fo','SiO2p','psi', ...
    'ol','an','qz','s_cpx','s_gtlhzhz','domain', ...
    'Fbm','Fafm','Fproj'};
melt_old = melt_init;

afm_flag = false;   % flag indicating AFM solution has been reached
bm_flag = false;    % flag indicating batch solution has been reached

% will stop olivine addition when the difference between the melt
% projections change sign
sign_afm = sign(melt_init.Fproj - melt_init.Fafm);  % sign parameter for AFM solution
sign_bm = sign(melt_init.Fproj - melt_init.Fbm);   % sign parameter for batch solution

maxit = 60;
melt_new = melt_init;
olivine_old = olivine_init;
olivine_new = olivine_init;
iter = 1;
while 1
    % set old to previous new
    melt_old = melt_new;
    olivine_old = olivine_new;
    
    % increment olivine added
    olivine_added = olivine_added + ol_step * 100;
    melt_new.olivine_added = olivine_added;
    
    % compute new melt composition
    melt_new{1,oxides} = (melt_old{1,oxides} + ...
        olivine_old{1,oxides} * ol_step) / (1 + ol_step);

    % normalize cations
    melt_new = cation_normalize(melt_new,1);

    % compute temperature of melt
    melt_new = temperature_estimator(melt_new);

    % compute KD, Fo and SiO2p from new melt composition
    melt_new = kd_iterator(melt_new);

    % compute olivine composition in equilibrium with new melt composition
    olivine_new = equilibrium_olivine(melt_new);

    % compute melt ternary composition
    melt_new = mineralfrac(melt_new);
    
    % compute melt fractions
    melt_new = meltfrac(melt_new, src);
    
    %if phi_afm_new >= phi_afm_old && ~afm_flag
    if sign_afm ~= sign(melt_new.Fproj - melt_new.Fafm) && ~afm_flag
        % initialize afm table
        afm = melt_old;

        % project olivine added to Fafm - Fproj = 0
        afm.olivine_added = -(melt_new.olivine_added - melt_old.olivine_added) * ...
            (melt_new.Fafm - melt_new.Fproj) / ...
            ((melt_new.Fafm - melt_new.Fproj) - (melt_old.Fafm - melt_old.Fproj)) + ...
            melt_new.olivine_added;

        % average other fields
        afm{:,fields} = avg_fields(melt_old{:,fields}, melt_new{:,fields}, ...
            melt_old.olivine_added, melt_new.olivine_added, afm.olivine_added);

        % lookup lithology (based on domains)
        if (afm.ol > 0.5 && afm.qz < afm.s_gtlhzhz) % domain 1
            afm.residual_lithology = {'garnet peridotite'};
        elseif afm.qz > afm.s_cpx % domain 3
            afm.residual_lithology = {'harzburgite'};
        else % domain 2
            afm.residual_lithology = {'spinel peridotite'};
        end
        
        % set flag (AFM iterations are done!)
        afm_flag = true;
    end
    
    %if phi_bm_new >= phi_bm_old && ~bm_flag
    if sign_bm ~= sign(melt_new.Fproj - melt_new.Fbm) && ~bm_flag
        % initialize bm table
        bm = melt_old;

        % project olivine added to Fbm - Fproj = 0
        bm.olivine_added = -(melt_new.olivine_added - melt_old.olivine_added) * ...
            (melt_new.Fbm - melt_new.Fproj) / ...
            ((melt_new.Fbm - melt_new.Fproj) - (melt_old.Fbm - melt_old.Fproj)) + ...
            melt_new.olivine_added;

        % average other fields
        bm{:,fields} = avg_fields(melt_old{:,fields}, melt_new{:,fields}, ...
            melt_old.olivine_added, melt_new.olivine_added, bm.olivine_added);

        % lookup lithology (based on domains)
        if (bm.ol > 0.5 && bm.qz < bm.s_gtlhzhz) % domain 1
            bm.residual_lithology = {'garnet peridotite'};
        elseif bm.qz > bm.s_cpx % domain 3
            bm.residual_lithology = {'harzburgite'};
        else % domain 2
            bm.residual_lithology = {'spinel peridotite'};
        end

        % set flag (BM iterations are done!)
        bm_flag = true;
    end
    % if prim_afm == 0
    %     warning('No AFM solution found.');
    % end

    if afm_flag && bm_flag
        break
    end
    
    if iter > maxit
        if ~exist('afm','var')
            afm = melt_old;
            afm.olivine_added = NaN;
            afm.residual_lithology = {''};
            afm.fe_ratio = NaN;
            afm.fe2o3_tio2 = NaN;
            afm.potential_temperature = NaN;
        end
        if ~exist('bm','var')
            bm = melt_old;
            bm.olivine_added = NaN;
            bm.residual_lithology = {''};            
        end
        afm.fe_ratio = NaN;
        afm.fe2o3_tio2 = NaN;
        afm.potential_temperature = NaN;
            
        bm.fe_ratio = NaN;
        bm.fe2o3_tio2 = NaN;
        bm.potential_temperature = NaN;
        return
    end
    
    iter = iter + 1;
end

% correct cation ratios
bm = cation_normalize(bm);
afm = cation_normalize(afm);

% correct Fe2+/(Fe2+ + Fe3+) ratio 
bm.fe_ratio = bm.feo/molecularwt('FeO') / ...
    (bm.feo/molecularwt('FeO') + 2*bm.fe2o3/molecularwt('Fe2O3'));
afm.fe_ratio = afm.feo/molecularwt('FeO') / ...
    (afm.feo/molecularwt('FeO') + 2*afm.fe2o3/molecularwt('Fe2O3'));

% correct Fe2O3/TiO2 ratio
bm.fe2o3_tio2 = bm.fe2o3/bm.tio2;
afm.fe2o3_tio2 = afm.fe2o3/afm.tio2;

afm = potential_temperature(afm);
bm = potential_temperature(bm);

%t = t(ip)*(prim_afm - ol_percent(ip-1)) + t(ip-1)*(ol_percent(ip) - prim_afm) / (ol_percent(ip) - ol_percent(ip-1));
% if abs(t.Fafm - proj_f1) < 0.0000001
%     fractional.residual = 'Harzburgite';
% elseif ABS(Fafm - proj_f2) < 0.000001
%     fractional.residual = 'Spinel Peridotite';
% else
%     fractional.residual = 'Garnet Peridotite';
% end

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

D{:,oxides(2:end)} = 100 * (c(1,:) * 2/3 * t.Fo/t.Mg + ...
    c(2,:) + ...
    c(3,:) .* t.Mg / (2 * t.Fo / 3)) .* t{:,CATIONS(2:end)};

sum1 = 0.5*(D.al2o3 - D.cr2o3 - D.na2o - D.k2o - D.p2o5) + ...
    D.cr2o3 + D.fe2o3 + D.feo + D.mno + D.mgo + D.cao + ...
    2*(D.na2o + D.k2o + D.p2o5) + D.nio;

% cation fraction olivine?
olivine = t(:,oxides);   % initialize olivine_eq
olivine{:,oxides} = NaN;

easy = {'tio2','al2o3','cr2o3','fe2o3','mno','cao','na2o','k2o','p2o5'};
olivine{:,easy}  = D{:,easy} * 200/3 / sum1;

sum2 = 0.5*(olivine.al2o3 - olivine.cr2o3 - olivine.na2o - ...
    olivine.k2o - olivine.p2o5) + olivine.cr2o3 + ...
    olivine.fe2o3 + olivine.mno + olivine.cao + ...
    2*(olivine.na2o + olivine.k2o + olivine.p2o5) + D.nio;

olivine.sio2 = 100/3 - olivine.tio2 - 0.5*(olivine.al2o3 - ...
    olivine.cr2o3 - olivine.na2o - olivine.k2o - ...
    olivine.p2o5) - olivine.cr2o3;

olivine.mgo = t.Fo * (200/3 - sum2);
olivine.feo = (1 - t.Fo) * (200/3 - sum2);
olivine.nio = (3.346 * (olivine.mgo / t.Mg) / 100 - 3.665) * t.Ni * 100;

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

KD2 = 0.3813 - 0.7896./t.mgo + 1.0389./t.mgo.^2;

tol = 1e-4;
maxit = 50;

t.Fbm = 0;
if t.feo ~= 0
    fb = (t.mgo * (src.feo / t.feo) - KD2 * src.mgo) / (t.mgo * (src.feo / t.feo) - t.mgo * KD2);
    fa = 0;
    fc = 0;
    iter = 1;
    while abs(fc - fb) > tol  
        fa = (t.mgo * (((src.feo - fb * t.feo) / (1 - fb)) / t.feo) - KD2 * src.mgo) ...
            / (t.mgo * (((src.feo - fb * t.feo) / (1 - fb)) / t.feo) - t.mgo * KD2);

        if fc ~= 0 && (fa - fb) / (fc - fb) > 0.9
            fa = (fa + fb) / 2;
        end
        fc = fb;
        fb = fa;

        iter = iter + 1;
        if iter > maxit
            warning('Fbm: maximum number of iterations exceeded.');
            break
        end
    end
    t.Fbm = fa;
end

if t.Fbm == 0
    t.Fbm = 1E-05;
end

if t.Fbm < 0
    t.Fafm = t.Fbm;
else
    t.Fafm = 1/(0.98 / t.Fbm - 0.9 * log(t.Fbm)^2 + 0.07 * log(t.Fbm)^4);
end

if t.domain == 1
    % F1
    t.Fproj = 6.2819 * t.an.^2 - 14.7789 * t.an.^3 + 0.00825 * (1 ./ t.an).^2;
elseif t.domain == 2
    % F2
    t.Fproj = ((-1.994 + 2.25 * t.qz + 0.041 / t.qz) + ...
        (-1.183 - 3.005 * t.qz + 13.774 * t.qz.^2 - 12.615 * t.qz.^3)) / 2 + ...
        exp(0.931 + 1.623*t.qz) * t.qz^0.245 * t.ol + ...
        exp(0.769 - 7.514 * t.qz) * t.qz^0.577 / t.ol;
    t.Fproj(t.qz <= 0) = 0;

    if t.Fproj > 0 && t.qz < -0.1773 + 0.1404 / t.ol - 0.008434 / t.ol^2
        t.Fproj = -t.Fproj;
    end
else
    % F3
    t.Fproj = -2.5345 + 5.329 * (t.qz + t.ol*0.348) + 0.3012 / (t.qz + 0.348 * t.ol);
    if t.Fproj > 0 && t.qz < -0.0896 + 0.02002 / t.ol + 0.02989 / t.ol^2
        t.Fproj = -t.Fproj;
    end
end

return


function t = mineralfrac(t)

global OXIDES MW

molefrac = t(:,lower(OXIDES));

molefrac{:,lower(OXIDES)} = t{:,lower(OXIDES)} ./ MW;
molefrac{:,lower(OXIDES)} = 100 * molefrac{:,lower(OXIDES)} ./ ...
    repmat(sum(molefrac{:,lower(OXIDES)},2),1,length(MW));

ol = -1.5*molefrac.tio2 + 0.5*(molefrac.al2o3 + molefrac.cr2o3 + ...
    molefrac.feo + molefrac.mno + molefrac.mgo + molefrac.nio) - ...
    0.5*molefrac.cao - 0.5*molefrac.na2o + 3*molefrac.k2o;
an = molefrac.tio2 + molefrac.al2o3 + molefrac.cr2o3;
qz = molefrac.sio2 - 0.5*(molefrac.al2o3 + molefrac.cr2o3 + ...
    molefrac.feo + molefrac.mno + molefrac.mgo + molefrac.nio) - ...
    1.5*molefrac.cao - 3*(molefrac.na2o + molefrac.k2o);

sum_min = (ol + an + qz);
t.ol = ol./sum_min;
t.an = an./sum_min;
t.qz = qz./sum_min;

t.s_cpx = -0.074 + 0.1713./t.ol - 0.0135./t.ol.^2;

if t.an < 700
    t.s_gtlhzhz = 1./(16.843 + 28.733*t.an - 14.183*exp(t.an));
else
    t.s_gtlhzhz = 0;
end

% determine compositional domain
if t.ol > 0.5 & t.qz < t.s_gtlhzhz
    t.domain = 3;
elseif t.qz > t.s_cpx
    t.domain = 1;
else
    t.domain = 2;
end

return


% mantle potential temperature
function t = potential_temperature(t);

t.potential_temperature = 1025 + 28.6*t.mgo - 0.084*t.mgo.^2;

return


% averaging function for table fields all at once
function t_avg = avg_fields(t1,t2,ol1,ol2,ol0)

t_avg = (t2 - t1)/(ol2 - ol1)*(ol0 - ol1) + t1;

return
