function A = heron(u,v,w)
% HERON - Computes area of non-right triangle using Heron's formula.
%
%   A = heron(a,b,c) where a, b, and c are the lengths of the triangle
%   sides.
%
%   A = heron(u,v,w) where u, v, and w are cartesian positions of the 
%   triangle vertices.  The vertices can come from a array with two columns
%   for (x,y) positions and each row representing a vertex from the
%   triangle.


if length(u(1,:)) == 1 && length(v(1,:)) == 1 && length(w(1,:)) == 1
    a = u;
    b = v;
    c = w;
else
    a = sqrt(sum((u - v).^2,2));
    b = sqrt(sum((v - w).^2,2));
    c = sqrt(sum((w - u).^2,2));
end

s = 0.5*(a + b + c);

A = sqrt(s.*abs(s - a).*abs(s - b).*abs(s - c));

return

