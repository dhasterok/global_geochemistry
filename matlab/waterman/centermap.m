function [map,cb] = centermap(map,i0,v0,clim)
% CENTERMAP - splits colormap about a value
%
% [map,cb] = centermap(map,i0,v0,clim);

du = clim(2) - v0;
dl = v0 - clim(1);

n = size(map,1);
nu = n - i0;
nl = i0;

% check to see if colormap is already split at v0
if abs(du/dl - nu/nl) < 1e-2
    return
end

if du/dl > nu/nl    % nu is too small
    nu_new = round(nl*du/dl)
    nt = nl + nu_new
    nl
    newmap = zeros(nt,3);
    newmap(1:nl) = map(1:nl,:);
    for k = 1:3
        newmap(nl+1:nt,k) = interp1([nl+1:n]',map(nl+1:n,k),linspace(nl+1,n,nt-nl)','linear');
    end
else    % nl is too small
    nl_new = round(nu*dl/du);
    nt = nl_new + nu;

    newmap = zeros(nt,3);
    newmap(nt-nu:nt,:) = map(n-nu:n,:);
    for k = 1:3
        newmap(1:nt-nu-1,k) = interp1([1:nl],map(1:nl,k),linspace(1,nl,nt-nu-1)','linear');
    end
end

colormap(newmap);
caxis(clim);
cb = colorbar;

return