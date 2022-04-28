function data = lambdaree(data,varargin)
% LAMBDAREE - Orthonormal polynomial approximation of REE concentrations
%
%   data = lambdaree(data) produces a set of orthonormal polynomial
%   constants to fit REE patterns (ignoring Eu).  The result are 5 lambda
%   constants lambda0 to lambda4 and a goodness of fit estimate
%   (MSWD value) lambda_mswd.
%
%   This method is based on the paper by H. O'Neill (J. Petrol., 2016).
%   The constants can be used to reveal processes and/or source mineralogy
%   from these patterns.

% to fit REE concentrations (in parts per million by weight) to a
% fourth-order orthogonal polynomial
% Original algorithm: Hugh O'Neill December 2013 in Excel
% Translated to Matlab by D. Hasterok 2020

nREE = 14;
nlambda = 5;
normalize = true;

ree = {'la_ppm','ce_ppm','pr_ppm','nd_ppm','sm_ppm','eu_ppm','gd_ppm', ...
    'tb_ppm','dy_ppm','ho_ppm','er_ppm','tm_ppm','yb_ppm','lu_ppm'};

indeu = strcmp(ree,'eu_ppm');

% REE ionic radii in Å from Shannon (1976)
ionic_radius = [1.16, 1.143, 1.126, 1.109, 1.079, 1.066, 1.053, ...
    1.04, 1.027, 1.015, 1.004, 0.994, 0.985, 0.977];

% CI Chondrite normalization values O'Neill (J. Petrol. 2016)
CI = [0.2472, 0.6308, 0.095, 0.4793, 0.15419, 0.0592, 0.2059, ...
	0.0375, 0.254, 0.0554, 0.1645, 0.0258, 0.1684, 0.0251];

% initialize fields
data.lambda0 = nan(height(data),1);
data.lambda1 = nan(height(data),1);
data.lambda2 = nan(height(data),1);
data.lambda3 = nan(height(data),1);
data.lambda4 = nan(height(data),1);
data.nree = nan(height(data),1);
data.chi2_lambda = nan(height(data),1);

% extract REE data and store in a matrix for faster sampling
sample = data{:,ree};
sample(sample <= 0) = nan;

% normalize REE concentrations?
if normalize
    sample = log(sample ./ repmat(CI,height(data),1));
else
    sample = log(sample);
end

% constant arrays used in a couple of different places
B1 = (ionic_radius - 1.054769);
B2 = (ionic_radius - 1.005327429) .* (ionic_radius - 1.128236038);
B3 = (ionic_radius - 1.060548105) .* (ionic_radius - 1.145519887) .* ...
    (ionic_radius - 0.991412204);
B4 = (ionic_radius - 1.104414973) .* (ionic_radius - 1.153426708) .* ...
    (ionic_radius - 0.984820219) .* (ionic_radius - 1.030518142);

% used for building the operator matrix
A1 = ionic_radius'.^[0:4]';
A2 = A1 .* repmat(B1,5,1);
A3 = A1 .* repmat(B2,5,1);
A4 = A1 .* repmat(B3,5,1);
A5 = A1 .* repmat(B4,5,1);

% initialize storage variabiles
lambda = nan(height(data),nlambda);
Nree = nan(height(data),1);
chi2_lambda = nan(height(data),1);
for k = 1:height(data)
    % extract sample REE data
    X = sample(k,:);
    % must have a minimum of 7 REE in order to do this modeling
    if sum(~isnan(X)) < 7
        continue
    end

    % number of REE used, exclude Eu
    ind = ~isnan(X) & ~indeu;
    Nree(k) = sum(ind);
    
    % construct operator matrix
    A = [sum(A1(:,ind),2) sum(A2(:,ind),2) sum(A3(:,ind),2) sum(A4(:,ind),2) sum(A5(:,ind),2)];

    % construct data matrix
    d = sum(repmat(X(ind),nlambda,1) .* A1(:,ind),2);
    
    % solve for lambda values
    % A \ d is more stable than inv(A)*d
    Q = A \ d;
    lambda(k,:) = Q';
    
    % calculate model REE
    Xm = sum(repmat(Q,1,nREE) .* [ones(1,nREE); B1; B2; B3; B4],1);
    
    % reduced chi-sq or MSWD misfit extimates
    chi2_lambda(k) = sum((X(ind) - Xm(ind)).^2) / (0.01^2 * (sum(ind) - 4));
end

% save results to data table
data{:,{'lambda0','lambda1','lambda2','lambda3','lambda4'}} = lambda;
data.nree = Nree;
data.lambda_MSWD = chi2_lambda;

return
