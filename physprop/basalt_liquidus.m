function T = basalt_liquidus(data)
% Estimating liquidus temperature for basaltic lavas using the empirical
% relation in Niu, Y., T. Gilmore, S. Mackie, A. Greig & W. Bach, Mineral
% chemistry, whole-rock compositions and petrogenesis of ODP Leg 176
% gabbros: Data and discussion, Proc. Ocean Drill. Prog. Sci. Results, 176,
% 1--60, 2002.
%
% This is only reasonable for SiO2 <60% wt.%



T = 1026 * exp(0.01894 * data.mgo);

% Method by Beattie
% Herzberg et al. G^3 2007
%T = 935 + 33*data.mgo - 0.37*data.mgo.^2;

% Herzberg & O'Hara, 2002 (potential temperature
%T = 1463 + 12.74*data.mgo - 2924./data.mgo;

return