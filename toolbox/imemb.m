function ib = imemb(A,B)

[tf,ia] = ismember(B,A);

% vector of all indicies into B
id = [1:numel(B)];

% indicies of B that match values in A
ib = accumarray(nonzeros(ia), id(tf), [], @(x){x});

return