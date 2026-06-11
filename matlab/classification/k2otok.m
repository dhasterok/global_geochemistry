function K = k2otok(K2O)
% K2OTOK - Converts wt% K2O to ppm K

mK = 39.102;
mO = 15.9994;

K = K2O/1e-4*2*mK/(mK*2 + mO);

return
