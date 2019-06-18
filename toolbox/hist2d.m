function [n,outliers,cx,cy] = hist2d(x,y,ex,ey);
% HIST2D - Computes two dimensional histograms.
%    N = HIST2D(X,Y,EX,EY) bins the elements of X and Y into
%    rectangular cells with borders defined by the edges
%    defined in EX and EY.  LENGTH(X) must equal LENGTH(Y);
%    however, LENGTH(EX) need not equal LENGTH(EY) but both
%    must be monotonically increasing.
%
%    [N,OUT] = HIST2D(X,Y,EX,EY) returns the number of
%    outliers not found in any of the requested bins as the
%    variable OUT.
%
%    [N,OUT,CX,CY] = HIST2D(X,Y,EX,EY) returns the center of
%    the bin locations as CX and CY, a row vector and column
%    vector respectively for use in MESHGRID or other matlab
%    plotting routines.

% Modified: 08-JUL-2003
% (c) D. Hasterok - University of Utah
%     dhasterok@mines.utah.edu

if length(x) ~= length(y)
    error('ERROR: Length of X and Y must be equal');
end

lx  = length(x);
ly  = length(y);

lex = length(ex);
if lex == 1;
    lex = ex;
    ex = [min(x),max(x),lex];
end

ley = length(ey);
if ley == 1;
    ley = ey;
    ey = [min(y),max(y),ley];
end

%h = waitbar(0,'Histogram Status...');
n = zeros(length(ey)-1,length(ex)-1);
for i = 1:ley-1
    for j = 1:lex-1        
        % Find x's and y's within rectangular bins
        indx = find(ex(j) < x & x <= ex(j+1));
        indy = find(ey(i) < y & y <= ey(i+1));
        
        n(i,j) = length(intersect(indx,indy));

        % Status Bar
        %waitbar((j + (i - 1)*j)/((lex - 1)*(ley - 1)),h);
    end
end
%close(h);

% Coordinates of bin centers
cx(1:lex-1) = [ex(2:lex)+ex(1:lex-1)]/2;
cx = cx(:)';
cy(1:ley-1) = [ey(2:ley)+ey(1:ley-1)]/2;
cy = cy(:);
        
if sum(sum(n)) ~= length(x)
    outliers = length(x) - sum(sum(n));
else
    outliers = 0;
end

return
