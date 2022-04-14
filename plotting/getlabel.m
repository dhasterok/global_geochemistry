function lbl = getlabel(el)
% GETLABEL - Creates plot label from geochemistry text
%
%   lbl = getlabel(el) where el is a string.  Units are added.  For ratios,
%   el should be a cell array where ratio = el{1}/el{2}

if iscell(el)
    lbl = [el{1},'/',el{2}];
else
    switch lower(el)
        case {'sio2','tio2','al2o3','feo_tot','mgo','cao','mno','na2o','k2o','p2o5','loi'}
            lbl = [el,' (wt.%)'];
        case {'heat production','heat_production'}
            lbl = 'Heat Production (µW m^{-3})';
        case {'density model','density_model'}
            lbl = 'Density (kg m^{-3})';
        case {'vmoho','moho','crustal_thickness','crustal thickness'}
            lbl = 'Crustal Thickness (km)';
        case {'p_velocity','s_velocity'}
            lbl = [upper(el(1)),' Velocity (km/s)']; 
        case {'mg_number'}
            lbl = 'Mg Number';
        otherwise
            lbl = [el,' (ppm)'];
    end
end

return