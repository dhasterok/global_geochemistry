function y = linInterp(x,a0,a1,b0,b1)

y = (x - a0).*(b1 - b0)./(a1 - a0) + b0;

return