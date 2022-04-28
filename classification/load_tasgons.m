function [tasgons,uhmgons] = load_tasgons

% (SiO2, Na2O+K2O) verticies for polygon TAS field names
% first name is volcanic name, second is plutonic
tasgons = {[41 0; 41 3; 33 3; 33 0], ...
    '', 'peridotite';
    [41 0; 41 3; 45 3; 45 0; 41 0], ...
    'picrobasalt', 'peridotgabbro';
    [45 2; 45 5; 52 5; 45 2], ...
    'alkalic basalt', 'alkalic gabbro';
    [45 0; 45 2; 52 5; 52 0; 45 0], ...
    'subalkalic basalt', 'subalkalic gabbro';
    [52 0; 52 5; 57 5.9; 57 0; 52 0], ...
    'basaltic andesite', 'gabbroic diorite';
    [57 0; 57 5.9; 63 7; 63 0; 57 0], ...
    'andesite', 'diorite';
    [63 0; 63 7; 69 8; 77.3 0; 63 0], ...
    'dacite', 'granodiorite';
    [77.3 0; 69 8; 71.8 13.5; 85.9  6.8; 87.5  4.7; 77.3 0], ...
    'rhyolite', 'granite';
    [45 5; 49.4 7.3; 52 5; 45 5], ...
    'trachybasalt', 'monzogabbro';
    [52 5; 49.4 7.3; 53 9.3; 57 5.9; 52 5], ...
    'basaltic trachyandesite', 'monzodiorite';
    [57 5.9; 53 9.3; 57.6 11.7; 63 7; 57 5.9], ...
    'trachyandesite', 'monzonite';
    [63 7; 61.1 8.65; 71.8 13.5; 69 8; 63 7], ...
    'trachydacite', 'quartz monzonite';
    [61.1 8.65; 57.6 11.7; 61 13.5; 63 16.2; 71.8 13.5; 61.1 8.65], ...
    'trachyte', 'syenite';
    [41 3; 41 7; 45 9.4; 49.4 7.3; 45 5; 45 3; 41 3], ...
    'tephrite', 'foid gabbro';
    [49.4 7.3; 45 9.4; 48.4 11.5; 53 9.3; 49.4 7.3], ...
    'phonotephrite', 'foid monzodiorite';
    [53 9.3; 48.4 11.5; 52.5 14; 57.6 11.7; 53 9.3], ...
    'tephriphonolite', 'foid monzosyenite';
    [57.6 11.7; 52.5 14; 52.5 18; 57 18; 63 16.2; 61 13.5; 57.6 11.7], ...
    'phonolite', 'foid syenite';
    [33 10; 33 3; 41 3; 41 7; 33 10], ...
    'ultramafic foidite','ultramafic foidolite';
    [41 7; 33 10; 37 14; 45 9.4; 41 7], ...
    'mafic foidite', 'mafic foidolite';
    [37 14; 52.5 18; 52.5 14; 48.4 11.5; 45 9.4; 37 14], ...
    'intermediate foidite', 'intermediate foidolite';
    [77.3 0; 87.5 4.7; 85.9 6.8; 93.2 6.8; 100 0; 77.3 0], ...
    'silexite', 'quartzolite';
    [93.2 6.8; 85.9 6.8; 71.8 13.5; 63 16.2; 57 18; 52.5 18; 37 14; 37 63; 93.2 6.8], ...
    'ultra-high alkali volcanic','ultra-high alkali plutonic'};
    %[33 3; 33 7; 37 14; 52.5 18; 52.5 14; 48.4 11.5; 45 9.4; 41 3; 33 3], ...
    %'foidite', 'foidolite';

% TAS polygons adapted for high-Mg lavas
uhmgons = {[33 0; 33 3; 52 3; 52 0; 33 0], 'picrite', '', 'intrusive picrite';
    [33 3; 33 6; 52 6; 52 3; 33 3], 'alkali picrite', '', 'intrusive alkali picrite';
    [33 0; 33 1; 45 1; 45 0; 33 0], 'komatiite', 'meimechite', 'intrusive komatiite';
    [45 0; 45 1; 52 1; 52 0; 45 0], 'basaltic komatiite', 'meimechite', 'intrusive komatiite';
    [52 0; 52 6; 65 6; 65 0; 52 0], 'boninite', '', 'sanukitoid'};

return