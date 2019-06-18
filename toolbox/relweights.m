function epsilon = relweights(A,d,rsq)
% RELWEIGHTS - Relative weights.
%
%   Computes the relative weights of parameters used for multiple linear
%   regression.
%
%   [beta,Lambda,epsilon] = relweights(A,d)

% Tonidandel, S., LeBreton, J.M., 2009. Determining the Relative Importance
% of Predictors in Logistic Regression: An Extension of Relative Weight
% Analysis. Organizational Research Methods, Organizational Research
% Methods 13, 767?781. https://doi.org/10.1177/1094428109341993

sum(A(:,1))
if sum(A(:,1)) == length(d)
    A = A(:,2:end);
    disp('first')
elseif sum(A(:,end)) == length(d)
    A = A(:,2:end);
    disp('last')
end

[U,S,V] = svd(A);

%Z = U*V';
Z = zscore(A)*V*diag(diag(1./sqrt(S)))*V';

T = inv(Z'*Z)*Z';
beta = T*d;

%Lambda = T*A;
Lambda = V*diag(diag(sqrt(S)))*V'

epsilon = Lambda^2*beta.^2;
scalefactor = rsq/sum(epsilon);
epsilon = epsilon*scalefactor;

return