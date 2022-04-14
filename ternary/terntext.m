function t = terntext(a,b,c,txt)
% TERNTEXT - places text within ternary axes.
%
%    t = terntext(a,b,c,txt)
%
%    The order of the axes are as follows:
%
%                    A
%                   / \
%                  /   \
%                 B --- C
%
%    See also TERNARY

[x,y] = tern2xy(a,b,c);

t = text(x,y,txt);

return
