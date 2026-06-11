function [K,Th,U]= hpe_age(tau,K,Th,U);

lambda_U235 = log(2)/7.04e2;
lambda_U238 = log(2)/4.468e3;
lambda_Th = log(2)/1.40e4;
lambda_K40 = log(2)/1.248e3;

U235 = 0.007204*U.*exp(lambda_U235*tau);
U238 = 0.992742*U.*exp(lambda_U238*tau);

U = U235 + U238;

Knr = (1 - 0.000117)*K;
K40 = 0.000117*K.*exp(lambda_K40*tau);

Ktotal = Knr + K40;

Th = Th.*exp(lambda_Th*tau);

return