function [sia_scheme,frost_class1,frost_class2,frost_class3] = graniteclass_sw(data);
% data = graniteclass(data);

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

% A-type rocks are ferroan
% Special conditions must apply above 70% SiO2 because of the significant
% overlap between I- and A- types
d1 = data.sio2 >= 70;

% for SiO2 < 70 wt.%
sia_scheme(ind2 & ~d1) = {'A'};

% for SiO2 >= 70 wt.%
%d2 = data.f_ppm < 1500;
d3 = data.k2o > data.na2o & data.k2o > 3;
%d4 = data.zr_ppm + data.nb_ppm + data.ce_ppm + data.y_ppm > 350;
%d5 = data.ga_ppm./(data.al2o3*2*molecularwt('Al')/molecularwt('Al2O3')) > 2.6;

sia_scheme(ind2 & d1 & d3) = {'A'};
%sia_scheme(ind2 & d1 & isnan(data.f_ppm) & d3) = {'A'};

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
