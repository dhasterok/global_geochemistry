function ttggons = load_ttggons

addpath ternary

[a1,b1,c1] = ternintersect([0 100],[70 0],[30 0],[25 0],[75 0],[0 100])
[a2,b2,c2] = ternintersect([0 80],[80 0],[20 20],[30 a1],[70 b1],[0 c1]);
[a3,b3,c3] = ternintersect([0 65],[65 0],[35 35],[0 a1],[0 b1],[100 c1]);
[a4,b4,c4] = ternintersect([0 50],[50 0],[50 50],[0 a1],[0 b1],[100 c1]);

% Normative An Ab Or
ttggons = {[0 100 0; 30 70 0; a1 b1 c1; 0 70 30; 0 100 0],'tronjhemite'; ...
    [30 70 0; 100 0 0; 80 0 20; a2 b2 c2; 30 70 0],'tonalite'; ...
    [80 0 20; a2 b2 c2; a1 b1 c1; a3 b3 c3; 65 0 35; 80 0 20],'granodiorite'; ...
    [0 70 30; 0 0 100; a1 b1 c1; 0 70 30],'leucogranite'; ...
    [65 0 35; a3 b3 c3; a4 b4 c4; 50 0 50; 65 0 35],'quartz monzonite'};

return