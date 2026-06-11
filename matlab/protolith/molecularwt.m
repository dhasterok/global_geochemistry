function wt = molecularwt(formula);
% MOLECULARWT - Computes the molecular weight for a given chemical
% formula.
%
%   WT = MOLECULARWT(F) where F is a string giving a chemical formula
%   unit (e.g. Fe2O3).  Correct capitalization is essential to correctly
%   interpreting the formula.  The molecular weight, WT, is returned in
%   g/mol.
%
%   If F is a cell array, each element of F is a separate formula
%   string.
%

if isa(formula,'cell');
    wt = zeros(size(formula));
    formula = formula(:);
    N = length(formula);
    for i = 1:N
        wt(i) = computemass(formula{i});
    end
else
    wt = computemass(formula);
end

return

    
function wt = computemass(formula);

wt = 0;
ind = 1;
len = length(formula);
while ind <= len;
    [mass, ind] = unitmass(formula,ind);
    wt = wt + mass;
end

return


function [wt, rind] = unitmass(formula,lind)

len = length(formula);

if strcmp(formula(lind),'(');   % ion
    lind = lind + 1;
    rind = lind + 1;
    while 1         % find ion formula
        if rind > len
            error(['ERROR: Error in chemical formula. Did not find closing parenthesis.']);
        end
        if strcmp(formula(rind),')');
            break;
        end
        rind = rind + 1;
    end
    rind = rind + 1;
    if lind > rind - 2
        error('ERROR: Empty parentheses in formula.');
    end
    ion = formula(lind:rind-2);

    % can find the mass of the ion by a recursive call
    mass = molecularwt(ion);
else                            % element
    rind = lind + 1;
    if rind <= len  % find element formula
        if strcmp(formula(rind),upper(formula(rind)))
            sym = formula(lind);
        else
            sym = formula(lind:rind);
            rind = rind + 1;
        end
    else
        sym = formula(lind);
    end
    % lookup mass of element
    mass = element_lookup(sym);
end

% how many moles of an elements or ion is there in the formula unit?
n = 0;
while 1
    if rind > len
        break;
    end
    if isnum(formula(rind));
        n = 10*n + str2num(formula(rind));
    else
        break;
    end
    rind = rind + 1;
end
% there must be at least one.
if n == 0
    n = 1;
end

wt = n*mass;

return


function mass = element_lookup(sym);

%      No. Name            Symbol Atomic weight
element = {1 'hydrogen'       'H'    1.00794;        2 'helium'         'He'   4.002602;
           3 'lithium'        'Li'   6.941;          4 'beryllium'      'Be'   9.012182;
           5 'boron'          'B'   10.811;          6 'carbon'         'C'   12.0107;
           7 'nitrogen'       'N'   14.0067;         8 'oxygen'         'O'   15.9994;
           9 'fluorine'       'F'   18.9984032;     10 'neon'           'Ne'  20.1797;
          11 'sodium'         'Na'  22.98976928;    12 'magnesium'      'Mg'  24.3050;
          13 'aluminum'       'Al'  26.9815386;     14 'silicon'        'Si'  28.0855;
          15 'phosphorus'     'P'   30.973762;      16 'sulfur'         'S'   32.065;
          17 'chlorine'       'Cl'  35.453;         18 'argon'          'Ar'  39.948;
          19 'potassium'      'K'   39.0983;        20 'calcium'        'Ca'  40.078;
          21 'scandium'       'Sc'  44.955912;      22 'titanium'       'Ti'  47.867;
          23 'vanadium'       'V'   50.9415;        24 'chromium'       'Cr'  51.9961;
          25 'manganese'      'Mn'  54.938045;      26 'iron'           'Fe'  55.845;
          27 'cobalt'         'Co'  58.933195;      28 'nickel'         'Ni'  58.6934;
          29 'copper'         'Cu'  63.546;         30 'zinc'           'Zn'  65.409;
          31 'gallium'        'Ga'  69.723;         32 'germanium'      'Ge'  72.64;
          33 'arsenic'        'As'  74.92160;       34 'selenium'       'Se'  78.96;
          35 'bromine'        'Br'  79.904;         36 'krypton'        'Kr'  83.798;
          37 'rubidium'       'Rb'  85.4678;        38 'strontium'      'Sr'  87.62;
          39 'yttrium'        'Y'   88.90585;       40 'zirconium'      'Zr'  91.224;
          41 'niobium'        'Nb'  92.90638;       42 'molybdenum'     'Mo'  95.94;
          44 'ruthenium'      'Ru' 101.07;          45 'rhodium'        'Rh' 102.90550;
          46 'palladium'      'Pd' 106.42;          47 'silver'         'Ag' 107.8682;
          48 'cadmium'        'Cd' 112.411;         49 'indium'         'In' 114.818;
          50 'tin'            'Sn' 118.710;         51 'antimony'       'Sb' 121.760;
          52 'tellurium'      'Te' 127.60;          53 'iodine'         'I'  126.90447;
          54 'xenon'          'Xe' 131.293;         55 'cesium'         'Cs' 132.9054519;
          56 'barium'         'Ba' 137.327;         57 'lanthanum'      'La' 138.90547;
          58 'cerium'         'Ce' 140.116;         59 'praseodymium'   'Pr' 140.90765;
          60 'neodymium'      'Nd' 144.242;         62 'samarium'       'Sm' 150.36;
          63 'europium'       'Eu' 151.964;         64 'gadolinium'     'Gd' 157.25;
          65 'terbium'        'Tb' 158.92535;       66 'dysprosium'     'Dy' 162.500;
          67 'holmium'        'Ho' 164.93032;       68 'erbium'         'Er' 167.259;
          69 'thulium'        'Tm' 168.93421;       70 'ytterbium'      'Yb' 173.04;
          71 'lutetium'       'Lu' 174.967;         72 'hafnium'        'Hf' 178.49;
          73 'tantalum'       'Ta' 180.94788;       74 'tungsten'       'W'  183.84;
          75 'rhenium'        'Re' 186.207;         76 'osmium'         'Os' 190.23;
          77 'iridium'        'Ir' 192.217;         78 'platinum'       'Pt' 195.084;
          79 'gold'           'Au' 196.966569;      80 'mercury'        'Hg' 200.59;
          81 'thallium'       'Tl' 204.3833;        82 'lead'           'Pb' 207.2;
          83 'bismuth'        'Bi' 208.98040;       90 'thorium'        'Th' 232.03806;
          91 'protactinium'   'Pa' 231.03588;       92 'uranium'        'U'  238.02891};

Ne = length(element);

i = 1;
for i = 1:Ne
    if strcmp(sym,element{i,3})
        mass = element{i,4};
        return;
    end
end
error(['ERROR: Element not found.  Either no natural isotopes ' ...
       'exist or there is an error in the formula.           ']);

return


function flag = isnum(str)

flag = 0;
num = ['0123456789'];
for i = 1:10
    if strcmp(str,num(i));
        flag = 1;
        return;
    end
end

return
