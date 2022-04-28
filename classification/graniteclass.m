function [sia_scheme,frost_class1,frost_class2,frost_class3] = graniteclass(data);
% GRANITECLASS - chemical classification systems for granites
%
%   [sia_scheme,frost_class1,frost_class2,frost_class3] = graniteclass(data)
%   classifies granites according to the Frost et al. (J. Petrol., 2001)
%   classification scheme and the SIA scheme of Chappell and White (Trans.
%   R. Soc. Edinburgh, 1992) and Whalen et al. (Contrib. Mineral. Petrol.,
%   1987).
%
%   The sia_scheme will return a cell array of 'S','I', or 'A'.  The frost
%   classes are divided into frost_class1 ('ferroan' or 'magnesian') as
%   determined from the iron-number, frost_class2 ('alkalic',
%   'alkali-calcic', 'calc-alkalic', or 'calcic') as determined by the
%   modified alkali-lime index (MALI), and frost_class3 ('peraluminous',
%   'peralkaline', or 'metaluminous') as determined by the alumina
%   saturation index (ASI).

% Last Modified: 27 April 2021
%   Changed definition of A-type granites.  Previous definition: below 70 
%   wt.% SiO2 used Frost ferroan class to define A-type.  Above 70 wt.%
%   used several indicators by Whalen et al. (1987) which all had to be
%   met.  These definitions may have been too restrictive.  Now use two
%   conditions: 
% by D. Hasterok, University of Adelaide

% Granite classification schemes
% --------------------------------------
% initalize SIA-scheme
sia_scheme = cell(size(data.sio2));
sia_scheme(:) = {''};

% initialize temporary for defining Frost classification
tmpfe = cell(size(data.sio2));
tmpfe(:) = {''};
tmpca = cell(size(data.sio2));
tmpca(:) = {''};
tmpasi = cell(size(data.sio2));
tmpasi(:) = {''};


% magnesian (iron poor) vs. ferroan (iron rich)
ind1 = (data.Fe_number < 0.486 + 0.0046*data.sio2);
tmpfe(ind1) = {'magnesian'};

ind2 = (data.Fe_number >= 0.486 + 0.0046*data.sio2);
tmpfe(ind2) = {'ferroan'};

% A-type rocks
d1 = data.sio2 < 70;
d2 = data.f_ppm < 1500;
d3 = data.AI < 1/0.85;  % agpaitic index, note defined in geochem_index.m
                        % as Al/(Na + K)
d4 = data.zr_ppm + data.nb_ppm + data.ce_ppm + data.y_ppm > 350;
d5 = data.ga_ppm./(data.al2o3*2*molecularwt('Al')/molecularwt('Al2O3')) > 2.6;

% post April 27, 2021
sia_scheme(ind2 & d1 & d3) = {'A'};
sia_scheme(~d1 & d3) = {'A'};

% Special conditions must apply above 70% SiO2 because of the significant
% overlap between I- and A- types
% pre April 27, 2021

% for SiO2 < 70 wt.%
%sia_scheme(ind2 & ~d1) = {'A'};

% for SiO2 >= 70 wt.%

%sia_scheme(ind2 & d1 & d2 & d3 & d4 & d5) = {'A'};
%sia_scheme(ind2 & d1 & isnan(data.f_ppm) & d3 & d4 & d5) = {'A'};


% alkalic vs. alkali-calcic vs. calc-alkalic vs. calcic
ind1 = (data.MALI > -41.86 + 1.12*data.sio2 - 0.00572*data.sio2.^2);
tmpca(ind1) = {'alkalic'};

ind2 = (data.MALI > -44.72 + 1.094*data.sio2 - 0.00527*data.sio2.^2);
tmpca(ind2 & ~ind1) = {'alkali-calcic'};

ind3 = (data.MALI > -45.36 + 1.0043*data.sio2 - 0.00427*data.sio2.^2);
tmpca(ind3 & ~ind2) = {'calc-alkalic'};

ind4 = (data.MALI <= -45.36 + 1.0043*data.sio2 - 0.00427*data.sio2.^2);
tmpca(ind4) = {'calcic'};

% peraluminous and peralkaline/metaluminous
ind = (data.ASI >= 1.0);
tmpasi(ind) = {'peraluminous'};

ind = (data.ASI < 1.0);
Al = 2*data.al2o3/molecularwt('Al2O3');

NaK = nansum(2*[data.na2o/molecularwt('Na2O'),data.k2o/molecularwt('K2O')],2);
indp = ind1 & (NaK > Al);
tmpasi(ind & indp) = {'peralkaline'};
tmpasi(ind & ~indp) = {'metaluminous'};

%frost_class = strcat(tmpfe,'>',tmpca,'>',tmpasi);
frost_class1 = tmpfe;
frost_class2 = tmpca;
frost_class3 = tmpasi;

ind = (data.ASI >= 1.1);
% S- and I-types
sia_scheme(ind & strcmp('',sia_scheme)) = {'S'};
sia_scheme(~ind & strcmp('',sia_scheme)) = {'I'};

return
