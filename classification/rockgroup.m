function ind = rockgroup(data,type)
% ROCKGROUP - Finds data for desired rock groups.
%
%   ind = ROCKGROUP(data,type) returns a vector of logicals containing
%   samples in DATA from the desired group in TYPE.
%
%   Rock types:
%       'igneous'
%       'all igneous' (igneous and metaigneous)
%       'igneous protolith' (igneous, metaigneous, and
%           predicted metaigneous)
%       'volcanic'
%       'all volcanic' (volcanic and metavolcanic)
%       'plutonic'
%       'all plutonic' (plutonic and metaplutonic)
%       'sedimentary'
%       'all seds' (sedimentary and metasedimentary)
%       'sedimentary protolith' (sedimentary, metasedimentary,
%           and predicted metasedimentary)
%       'metamorphic'
%       'metaigneous' (metaigneous, metaplutonic, and metavolcanic)
%       'metaplutonic'
%       'metavolcanic'
%       'metasedimentary'


switch type
    case 'igneous'
        ind = strcmpi('igneous',data.rock_group);
    case 'all igneous'
        ind = strcmpi('igneous',data.rock_group) ...
            | strcmpi('metaigneous',data.rock_origin) ...
            | strcmpi('metaplutonic',data.rock_origin) ...
            | strcmpi('metavolcanic',data.rock_origin);
    case 'igneous protolith'
        if any(strcmpi('protolith_est',data.Properties.VariableNames))
            ind = strcmpi('igneous',data.rock_group) ...
                | strcmpi('metaigneous',data.rock_origin) ...
                | strcmpi('metaplutonic',data.rock_origin) ...
                | strcmpi('metavolcanic',data.rock_origin) ...
                | (strcmpi('metamorphic',data.rock_group) ...
                & strcmpi('metaigneous',data.protolith_est) ...
                & ~strcmpi('metasedimentary',data.rock_origin));
        else
            warning('Asked for estimated protoliths, but they have not yet been determined.');
            ind = strcmpi('igneous',data.rock_group) ...
                | strcmpi('metaigneous',data.rock_origin) ...
                | strcmpi('metaplutonic',data.rock_origin) ...
                | strcmpi('metavolcanic',data.rock_origin);
        end
    case 'volcanic'
        ind = strcmpi('volcanic',data.rock_origin);
    case 'metavolcanic'
        ind = strcmpi('metavolcanic',data.rock_origin);
    case 'all volcanic'
        ind = strcmpi('volcanic',data.rock_origin) | strcmpi('metavolcanic',data.rock_origin);    
    case 'plutonic'
        ind = strcmpi('plutonic',data.rock_origin);
    case 'metaplutonic'
        ind = strcmpi('metaplutonic',data.rock_origin);
    case 'all plutonic'
        ind = strcmpi('plutonic',data.rock_origin) | strcmpi('metaplutonic',data.rock_origin);    
    case 'sedimentary'
        ind = strcmpi('sedimentary',data.rock_group);
    case 'all seds'
        ind = strcmpi('sedimentary',data.rock_group) ...
            | strcmpi('metasedimentary',data.rock_origin);
    case 'sedimentary protolith'
        if any(strcmpi('protolith_est',data.Properties.VariableNames))
            ind = strcmpi('sedimentary',data.rock_group) ...
                | strcmpi('metasedimentary',data.rock_origin) ...
                | (strcmpi('metamorphic',data.rock_group) ...
                & strcmpi('metasedimentary',data.protolith_est) ...
                & ~strcmpi('metaigneous',data.rock_origin) ...
                & ~strcmpi('metaplutonic',data.rock_origin) ...
                & ~strcmpi('metavolcanic',data.rock_origin));
        else
            warning('Asked for estimated protoliths, but they have not yet been determined.');
            ind = strcmpi('sedimentary',data.rock_group) ...
                | strcmpi('metasedimentary',data.rock_origin);
        end 
    case 'metamorphic'
        ind = strcmpi('metamorphic',data.rock_group);
    case 'metaigneous'
        ind = strcmpi('metaigneous',data.rock_origin) ...
            | strcmpi('metaplutonic',data.rock_origin) ...
            | strcmpi('metavolcanic',data.rock_origin);
    case 'metasedimentary'
        ind = strcmpi('metasedimentary',data.rock_origin);
    case 'mineral'
        if any(strcmp(data.Properties.VariableNames,'mineral'))
            ind = strcmpi('mineral',data.material);
        else
            ind = logical(zeros([height(data) 1]));
        end
    otherwise
        error('Unknown type.');
end
%unique(data.rock_origin(ind))

return