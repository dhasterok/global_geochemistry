function map = palatte(left,center,right,n,varargin);
% PALATTE - Creates a RGB colormap.
%   MAP = PALATTE(LEFT,CENTER,RIGHT,N) creates a colormap grading
%   from the color LEFT through CENTER to RIGHT in N steps.  LEFT,
%   CENTER, and RIGHT must be color RGB triplets from 0 to 1.
%
%   An optional output file may be given using the command
%   MAP = PALATTE(LEFT,CENTER,RIGHT,N,FILE).
%
%   Last Modified: 09-Mar 2004 by D. Hasterok

if mod(n,2) == 0
    n = n + 1;
end

rl = left(1);
gl = left(2);
bl = left(3);

rc = center(1);
gc = center(2);
bc = center(3);

rr = right(1);
gr = right(2);
br = right(3);

R = colorgrad(rl,rc,rr,n);
G = colorgrad(gl,gc,gr,n);
B = colorgrad(bl,bc,br,n);

map = [R G B];

if nargin == 5
    filename = varargin{1}
else
    filename = [];
end

if ~isempty(filename)
    fid = fopen(filename,'w');
    for i = 1:length(map(:,1))
        fprintf(fid,'%6.4f  %6.4f  %6.4f\n',map(i,:));
    end
    fclose(fid);
end

return


function C = colorgrad(cl,cc,cr,n)
C = zeros(n,1);

len = (n - 1)/2;

C([1:len])   = cl + (cc - cl)*[0:len-1]'/len;
C(len+1)     = cc;
C([len+2:n]) = cr - (cr - cc)*[len-1:-1:0]/len;

return