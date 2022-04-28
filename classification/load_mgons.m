function mgons = load_mgons

mgons = { ...
    [100   0   0; 90 10  0; 90  0 10; 100   0   0], 'anorthosite';
    [  0 100   0; 10 90  0;  0 90 10;   0 100   0], 'orthopyroxenite';
    [  0   0 100; 10  0 90;  0 10 90;   0   0 100], 'clinopyroxenite';
    [ 10  5 85; 5 5 90; 0 10 90; 0 90 10; 5 90 5; 10 85 5; 10 5 85], 'plagioclase-bearing websterite';
    [5 90 5; 10 90 0; 90 10 0; 90 5 5; 5 90 5], 'norite';
    [5 5 90; 10 0 90; 90 0 10; 90 5 5; 5 5 90], 'gabbro';
    [10 45 45; 90 5 5; 10 85 5; 10 45 45], 'orthopyroxene gabbro';
    [10 45 45; 90 5 5; 10 5 85; 10 45 45], 'clinopyroxene norite'};

return