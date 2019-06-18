function varargout = isect(A,B)
% ISECT - returns all indicies of A in B / B in A
%
% IA = isect(A,B) returns the indicies of all items in A which match items
% in B.
%
% To return both, [IA,IB] = isect(A,B)

[tf,ib] = ismember(A,B);

% vector of all indicies into A
id = [1:numel(A)];

% indicies of A that match values in B
ia = accumarray(nonzeros(ib), id(tf), [], @(x) {x} );
ia = cell2mat(ia);

varargout{1} = ia;

if nargout == 2
    % do it again! but the other way around
    % indicies of B that match values in A
    ib = isect(B,A);

    varargout{2} = ib;
end

return