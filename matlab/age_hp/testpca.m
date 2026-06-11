function testpca(data)


datamatrix = [data.sio2 data.tio2 data.al2o3 data.cr2o3 data.feo_tot data.mgo data.cao ...
    data.mno data.nio data.k2o data.na2o data.sro data.p2o5 data.h2o_plus ...
    data.co2 data.so3 data.bao];

y = data.heat_production;
y(any(isnan(datamatrix), 2), :) = [];

datamatrix(any(isnan(datamatrix), 2), :) = [];

datamatrix = [ones(size(datamatrix,1),1) datamatrix];

%[coeff,score,latent,tsquared,explained,mu] = pca(datamatrix);

% majors = {'sio2';'tio2';'al2o3';'feo_tot';'mgo';'cao';'na2o';'k2o';'p2o5'};
% 
% [coeff, score, latent, Tsquared, explained, el] = hp_pca(data,majors)

%NIR = 60*401 (60 models, of 401 points)
%OCTANE + 60x1 (60 models, one value)


b = regress(y,datamatrix)


end