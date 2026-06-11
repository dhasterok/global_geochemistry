function qapgons = load_qapgons

qapgons = { ...
    [90 10     0    0; 100   0     0    0; 90  0    10    0; 90 10     0    0],                   'quartzolite', 'silexite';
    [90 10     0    0;  90   0    10    0; 60  0    40    0; 60 40     0    0; 90 11     0    0], 'quartz-rich granitoid', 'quartz-rich rhyolite';
    [60 40     0    0;  60  36     4    0; 20 72     8    0; 20 80     0    0; 60 40     0    0], 'alkali feldspar granite', 'alkali feldspar rhyolite';
    [60  4    36    0;  20   8    72    0; 20 28    52    0; 60 14    26    0; 60  4    36    0], 'granodiorite', 'dacite';
    [60 26    14    0;  20  52    28    0; 20 72     8    0; 60 36     4    0; 60 26    14    0], 'syenogranite', 'rhyolite';
    [60 14    26    0;  20  28    52    0; 20 52    28    0; 60 26    14    0; 60 14    26    0], 'monzogranite', 'rhyolite';
    [60  0    40    0;  60   4    36    0; 20  8    72    0; 20  0    80    0; 60  0    40    0], 'tonalite', 'dacite';
    [20 80     0    0;   5  95     0    0;  5 85.5   9.5  0; 20 72     8    0],                   'alkali feldspar quartz syenite', 'alkali feldspar quartz trachyte';
    [20 72     8    0;   5  85.5   9.5  0;  5 61.75 33.25 0; 20 52    28    0; 20 72     8    0], 'quartz syenite', 'quartz trachyte';
    [20 52    28    0;   5  61.75 33.25 0;  5 33.25 61.75 0; 20 28    52    0; 20 52    28    0], 'quartz monzonite', 'quartz latite';
    [20 28    52    0;   5  33.25 61.75 0;  5  9.5  85.5  0; 20  8    72    0; 20 28    52    0], 'quartz monzodiorite', 'andesite';
    [20  8    72    0;   5   9.5  85.5  0;  5  0    95    0; 20  0    80    0; 20  8    72    0], 'quartz diorite', 'andesite'; % /gabbro/anorthosite
    [ 5 95     0    0;   0 100     0    0;  0 90    10    0;  5 85.5   9.5  0;  5 95     0    0], 'alkali feldspar syenite', 'alkali feldspar trachyte';
    [ 5 85.5   9.5  0;   0  90    10    0;  0 65    35    0;  5 61.75 33.25 0;  5 85.5   9.5  0], 'syenite', 'trachyte';
    [ 5 61.75 33.25 0;   0  65    35    0;  0 35    65    0;  5 33.25 61.75 0;  5 61.75 33.25 0], 'monzonite', 'latite';
    [ 5 33.25 61.75 0;   0  35    65    0;  0 10    90    0;  5  9.5  85.5  0;  5 33.25 61.75 0], 'monzodiorite', 'andesite'; % /monzogabbro
    [ 5  9.5  85.5  0;   0  10    90    0;  0  0   100    0;  5  0    95    0;  5  9.5  85.5  0], 'diorite', 'andesite'; % gabbro
    [ 0  0     0  100;   0  10     0   90;  0  0    10   90;  0  0     0  100],                   'foidolite', 'foidite';
    [ 0  10    0   90;   0  40     0   60;  0 20    20   60;  0  5     5   90;  0 10     0   90], 'foidolite', 'phonolitic foidite';
    [ 0  0    10   90;   0   0    40   60;  0 20    20   60;  0  5     5   90;  0  0    10   90], 'foidolite', 'tephritic foidite';
    [ 0 40     0   60;   0  36     4   60;  0 81     9   10;  0 90     0   10;  0 40     0   60], 'foid syenite', 'phonolite';
    [ 0 100    0    0;   0  90    10    0;  0 81     9   10;  0 90     0   10;  0 100    0    0], 'foid-bearing alkali feldspar syenite', 'foid-bearing alkali feldspar trachyte';
    [ 0 90    10    0;   0  81     9   10;  0 58.5  31.5 10;  0 65    35    0;  0 90    10    0], 'foid-bearing syenite', 'foid-bearing trachyte';
    [ 0 65    35    0;   0  35    65    0;  0 31.5  58.5 10;  0 58.5  31.5 10;  0 65    35    0], 'foid-bearing monzonite', 'foid-bearing latite';
    [ 0 10    90    0;   0   9    81   10;  0 31.5  58.5 10;  0 35    65    0;  0 10    90    0], 'foid-bearing monzodiorite', 'andesite'; % /monzogabbro
    [ 0  0   100    0;   0  10    90    0;  0  9    81   10;  0  0    90   10;  0  0   100    0], 'foid-bearing diorite', 'andesite'; % /gabbro/anorthosite
    [ 0  0    40   60;   0   4    36   60;  0  9    81   10;  0  0    90   10;  0  0    40   60], 'foid diorite', 'tephrite'; % /gabbro
    [ 0 45    45   10;   0   9    81   10;  0  4    36   60;  0 20    20   60;  0 45    45   10], 'foid monzodiorite', 'phonolitic tephrite'; % /monzogabbro
    [ 0 45    45   10;   0  81     9   10;  0 36     4   60;  0 20    20   60;  0 45    45   10], 'foid monzosyenite', 'tephritic phonolite'};

return
