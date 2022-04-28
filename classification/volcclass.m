function volcgons = volcclass

% A new cation plot for classifying sub-alkaline volcanic rocks
% LS Jenson - Ontario Div. Min., Misc. Paper, 1976
addpath ternary

[a1,b1,c1] = ternintersect([10 55],[90 0],[0 45],[0 20],[80 80],[20 0]);
[a2,b2,c2] = ternintersect([10 55],[90 0],[0 45],[0 30],[70 70],[30 0]);
[a3,b3,c3] = ternintersect([10 55],[90 0],[0 45],[0 40],[60 60],[40 0]);
[a4,b4,c4] = ternintersect([10 55],[90 0],[0 45],[0 49],[51 51],[49 0]);
[a5,b5,c5] = ternintersect([0 50],[100 0],[0 50],[50 0],[50 50],[0 50]);

% Normative An Ab Or
volcgons = {[0 40 60; 40 0 60; 0 0 100; 0 40 60],'komatiite'; ...
    [70 15 15; 0 50 50; 0 40 60; 40 0 60; 85 0 15; 70 15 15],'basaltic komatiite'; ...
    [70 15 15; 85 0 15; 100 0 0; 85 15 0; 70 15 15], 'high Ti-Fe other'; ...
    [10 90 0; a1 b1 c1; 0 80 20; 0 100 0; 10 90 0],'calc-alkaline rhyolite'; ...
    [10 90 0; a2 b2 c2; 30 70 0; 10 90 0],'tholeiitic rhyolite'; ...
    [a1 b1 c1; 0 80 20; 0 70 30; a2 b2 c2; a1 b1 c1],'calc-alkaline dacite'; ...
    [a2 b2 c2; 30 70 0; 40 60 0; a3 b3 c3; a2 b2 c2],'tholeiitic dacite'; ...
    [a2 b2 c2; 0 70 30; 0 60 40; a3 b3 c3; a2 b2 c2],'calc-alkaline andesite'; ...
    [a3 b3 c3; 40 60 0; 50 50 0; 32 50 18; a4 b4 c4; a3 b3 c3],'tholeiitic andesite'; ...
    [0 60 40; a3 b3 c3; a4 b4 c4; 28 50 22; a5 b5 c5; 13 51 36; 0 57.5 42.5; 0 60 40],'calc-alkaline basalt'; ...
    [50 50 0; 32 50 18; a4 b4 c4; 28 50 22; a5 b5 c5; 100/3 100/3 100/3; 70 15 15; 85 15 0; 50 50 0],'high Fe-tholeiitic basalt'; ...
    [0 50 50; 0 57.5 42.5; 13 51 36; a5 b5 c5; 100/3 100/3 100/3; 0 50 50],'high Mg-tholeiitic basalt'};

return