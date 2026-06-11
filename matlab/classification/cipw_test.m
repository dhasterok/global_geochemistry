close all;
clear all;

s = [58.88; % SiO2
      0.79; % TiO2
     15.29; % Al2O3
      0.87; % Fe2O3
      5.47; % FeO
      0.07; % MnO
      3.79; % MgO
      5.10; % CaO
      3.26; % Na2O
      3.01; % K2O
      0.50; % P2O5
      0   ; % CO2
      0   ; % SO3
      0   ; % S
      0   ; % F
      0   ; % Cl
      0   ; % Sr
      3200; % Ba
      0   ; % Ni
      0   ; % Cr
      0  ]; % Zr
  
minerals = cipw(s);
fields = fieldnames(minerals);
for i = 1:length(fields)
    temp = getfield(minerals,fields{i});
    
    if temp.moles ~= 0 & ~isnan(temp.moles)
        fprintf('%-15s %3.4f\n',fields{i},temp.mode);
    end
end