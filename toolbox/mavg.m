function r = mavg(x,w,varargin);
% MAVG - moving average
%
%   r = mavg(x,w) produces a moving average (r) of x with a window size of
%   w points.
if mod(w,2) == 0
    w = w + 1;
end

m = (w - 1)/2;

% convert to average along column dimension and convert back to row later
if any(size(x) == 1)% vector
    row = 0;
    if size(x,1) < size(x,2)
        row = 1;
        x = x(:);
    end
else                % matrix
    dim = 1;
    if nargin == 3
        dim = varargin{1};
    end
    
    if dim == 2
        x = x';
    end    
end

r = zeros(size(x));
n = size(x,1);

for i = 1:m
    r(i,:) = sum(x(1:m+i,:),1)/(w - m + i - 1);
    r(end-i+1,:) = sum(x(end-m-i+1:end,:),1)/(w - m + i -1);
end
for i = (m+1):(n-m)
    r(i,:) = sum(x(i-m:i+m,:),1) / w;
end

% if average requested along dimension
if any(size(x) == 1)
    if row
        r = r';
    end
else
    if dim == 2
        r = r';
    end
end

return
