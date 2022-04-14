function refchem = loadrefchem(varargin);

% load table of reference compositions
eref = readtable('earthref.xlsx');

if mod(nargin,2) == 1
    error('Invalid number of input arguments.');
end

%set options
mflag = 0;
feflag = 0;
physflag = 0;
for i = 1:nargin/2
    switch lower(varargin{2*i-1})
        case 'model' % select reference lines
            models = varargin{2*i};
            if isstr(models)
                refchem = eref(strcmp(eref.model,models),:);
            elseif iscell(models)
                ind = zeros([height(eref) 1]);
                for i = 1:length(models)
                    ind = ind & strcmp(eref.model,models{i});
                end
                refchem = eref(ind,:);
            end
            mflag = 1;
        case 'fefix' % flag to convert all FeO to FeOT
            feflag = varargin{2*i};
        case 'physprop' % compute physical properties.
            physflag = varargin{2*i};
    end
end

if ~mflag
    refchem = eref;
end

if feflag
    disp('Recomputing Fe to FeOT...');
    refchem = fefix(refchem);
end

if physflag & refchem.sigma == 0
    disp('Computing rock properties...');
    % geochemical indicies
    refchem = geochem_index(refchem);
    
    % estimate seismic velocities
    refchem = vpest(refchem);
    
    % estimate density
    refchem = densest2(refchem);

    % estimate heat production
    refchem = hpest(refchem);

    ind = (refchem.sigma ~= 0);

    refchem.p_velocity(ind) = NaN;
    refchem.heat_production_mass(ind) = NaN;
    refchem.heat_production(ind) = NaN;
    refchem.density_bk(ind) = NaN;
end

return
