function [Atot,summary] = radheat2(density,K,Rb,Sm,Th,U,age,varargin);
% check decay rates for Sm-147 and fix Sm for two radiogenic isotopes

opt = 0;
if nargin == 7
    opt = varargin{1};
end

if opt
    % convert wt.% K to ppm K
    K = K*1e4;
else
    % Convert from wt.% K2O to ppm K
    K = 2*molecularwt('K')/molecularwt('K2O') * K * 1e4;
end

% order of elements in summary
element = [K,Rb,Sm,Th,U];

% Constants and conversion factors
alpha_mass = 4.00260325413;  % amu | g/mol mass of alpha particle
alpha_uncertainty = 6e-11; % amu IAEA
% 1 eV = 1.60217653(14)e-19 J  i.e. charge of an electron
Ce = 1.6021766208e-19;
Ce_uncertainty = 98e-29;

eV2kg = 1.782661907e-36; % NIST (2014 CODATA)
eV2kg_uncertainty = 11e-45;
amu2kg = 1.660539040e-27; % NIST (2014 CODATA)
amu2kg_uncertainty = 20e-36;

% Avogadro's Number
Na = 6.022140857e23; % mol^-1 NIST (2014 CODATA)
Na_uncertainty = 74e-32;

% Seconds per year (365.2422 days/yr)
secperyr = 31556926.08; % s/yr

% convert time in Ma to seconds
time = age * 1e6 * secperyr;


% --------------------------------
% Radiogenic decay constant (in s)
% Decay constants and energies from
% UN International Atomic Energy Agency,
% Nuclear Data Services
% https://www-nds.iaea.org/relnsd/vcharthtml/VChartHTML.html
% 14 June 2017
% --------------------------------
% Data on atomic masses and isotopic abundance from
% US Physical measurement laboratory, NIST,
% http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
% 14 June 2017

% ratio of each mode of decay (only matters for K)
% 40-K->Ca 40-K->Ar 232-Th 235-U 238-U 87-Rb
summary(1).element = 'K';
summary(1).isotope = {'K-40'};
summary(1).abundance = 0.000117;
summary(1).abundance_uncertainty = 0.000001;
summary(1).daughter = {'Ca-40','Ar-40'};
summary(1).P_decay = [0.8928, 0.1072];
summary(1).P_decay_uncertainty = [0.0013, 0.0013];
summary(1).half_life = 1.248e9; % a
summary(1).half_life_uncertainty = 3e6; % a
summary(1).decay_energy = [0.56018, 1.460822]; % MeV
summary(1).decay_uncertainty = [0.00005, 0.000006]; % MeV Hasterok & Gard 2017
summary(1).decay_rate = 1.1559e-14;
summary(1).decay_rate_uncertainty = 0.0050e-14;
% summary(1).decay_energy = [0.590, 1.461]; % MeV Fiorentini et al., 2007
% summary(1).decay_energy = [0.60, 1.460]; % MeV Rybach, 1988
% summary(1).atomic_mass = 39.963998166; % NIST amu | g/mol
summary(1).atomic_mass = 39.96399817; % IAEA amu | g/mol
summary(1).atomic_mass_uncertainty = 6e-8; % IAEA amu | g/mol
% summary(1).daughter_mass = [39.962590863, 39.9623831237]; % NIST amu | g/mol
summary(1).daughter_mass = [39.962590865, 39.9623831238]; % IAEA amu | g/mol
summary(1).number_alpha = [0, 0];
summary(1).neutrino_energy = [0.75089, 0.043578]; % MeV Hasterok & Gard 2017
summary(1).neutrino_uncertainty = [0.00012, 0.000060]; % MeV Hasterok & Gard 2017

summary(2).element = 'Rb';
summary(2).isotope = {'Rb-87'};
summary(2).abundance = 0.2783;
summary(2).abundance_uncertainty = 0.0002;
summary(2).daughter = {'Sr-87'};
summary(2).P_decay = 1;
summary(2).P_decay_uncertainty = 0;
summary(2).half_life = 4.97e10; % a
summary(2).half_life_uncertainty = 3e8; % a
summary(2).decay_energy = 0.08167; % MeV Hasterok & Gard 2017
summary(2).decay_uncertainty = 0.00036; % MeV Hasterok & Gard 2017
summary(2).decay_rate = 3.609e-17;
summary(2).decay_rate_uncertainty = 0.027e-17;
% summary(2).decay_energy = 0.122; % MeV Fiorentini et al., 2007
% summary(2).atomic_mass = 86.9091805310; % NIST amu | g/mol
summary(2).atomic_mass = 86.909180531; % IAEA amu | g/mol
summary(2).atomic_mass_uncertainty = 6e-9; % IAEA amu | g/mol
% summary(2).daughter_mass = 86.9088775; % NIST amu | g/mol
summary(2).daughter_mass = 86.908877496; % IAEA amu | g/mol
summary(2).number_alpha = 0;
summary(2).neutrino_energy = 0.20061; % MeV Hasterok & Gard 2017
summary(2).neutrino_uncertainty = 0.00036; % MeV Hasterok & Gard 2017

summary(3).element = 'Sm';
summary(3).isotope = {'Sm-148'};
summary(3).abundance = 0.1124;
summary(3).abundance_uncertainty = 0.0010;
summary(3).daughter = {'Nd-144'};
summary(3).P_decay = 1;
summary(3).P_decay_uncertainty = 0;
summary(3).half_life = 7e15; % a
summary(3).half_life_uncertainty = 3e15; % a
summary(3).decay_energy = 1.9868; % MeV Hasterok & Gard 2017
summary(3).decay_uncertainty = 0.00004; % MeV Hasterok & Gard 2017
summary(3).decay_rate = 1e-16;
summary(3).decay_rate_uncertainty = 4e-17;
% summary(3).decay_energy = 0.122; % MeV Fiorentini et al., 2007
% summary(3).atomic_mass = 86.9091805310; % NIST amu | g/mol
summary(3).atomic_mass = 147.9148290; % IAEA amu | g/mol
summary(3).atomic_mass_uncertainty = 1.5e-6; % IAEA amu | g/mol
% summary(3).daughter_mass = 86.9088775; % NIST amu | g/mol
summary(3).daughter_mass = 143.9100929; % IAEA amu | g/mol
summary(3).number_alpha = 1;
summary(3).neutrino_energy = 0; % MeV Hasterok & Gard 2017
summary(3).neutrino_uncertainty = 0; % MeV Hasterok & Gard 2017

summary(3).element = 'Sm';
summary(3).isotope = {'Sm-147'};
summary(3).abundance = 0.1499;
summary(3).abundance_uncertainty = 0.0018;
summary(3).daughter = {'Nd-143'};
summary(3).P_decay = 1;
summary(3).P_decay_uncertainty = 0;
summary(3).half_life = 1.060e11; % a
summary(3).half_life_uncertainty = 0.011e11; % a
summary(3).decay_energy = 2.31096; % MeV Hasterok & Gard 2017
summary(3).decay_uncertainty = 0.00035; % MeV Hasterok & Gard 2017
summary(3).decay_rate = 1e-16;
summary(3).decay_rate_uncertainty = 4e-17;
% summary(3).decay_energy = 0.122; % MeV Fiorentini et al., 2007
% summary(3).atomic_mass = 86.9091805310; % NIST amu | g/mol
summary(3).atomic_mass = 146.9149044; % IAEA amu | g/mol
summary(3).atomic_mass_uncertainty = 1.9e-6; % IAEA amu | g/mol
% summary(3).daughter_mass = 86.9088775; % NIST amu | g/mol
summary(3).daughter_mass = 142.9098200; % IAEA amu | g/mol
summary(3).number_alpha = 1;
summary(3).neutrino_energy = 0; % MeV Hasterok & Gard 2017
summary(3).neutrino_uncertainty = 0; % MeV Hasterok & Gard 2017

summary(4).element = 'Th';
summary(4).isotope = {'Th-232'};
summary(4).abundance = 1;
summary(4).abundance_uncertainty = 0;
summary(4).daughter = {'Pb-208'};
summary(4).P_decay = 1;
summary(4).P_decay_uncertainty = 0;
summary(4).half_life = 1.40e10; % a
summary(4).half_life_uncertainty = 1e8; % a
summary(4).decay_energy = 38.3652; % MeV Hasterok & Gard 2017
summary(4).decay_uncertainty = 0.0071; % MeV Hasterok & Gard 2017
summary(4).decay_rate = 6.019e-14;
summary(4).decay_rate_uncertainty = 0.016e-14;
%summary(4).decay_energy = 40.4; % MeV Fiorentini et al., 2007
%summary(4).decay_energy = 38.99; % MeV Rybach, 1988
%summary(4).atomic_mass = 232.0380558; % NIST amu | g/mol
summary(4).atomic_mass = 232.0380537; % IAEA amu | g/mol
summary(4).atomic_mass_uncertainty = 1.5e-6; % IAEA amu | g/mol
%summary(4).daughter_mass = 207.9766525; % NIST amu | g/mol
summary(4).daughter_mass = 207.9766519; % IAEA amu | g/mol
summary(4).number_alpha = 6;
summary(4).neutrino_energy = 4.2805; % MeV Hasterok & Gard 2017
summary(4).neutrino_uncertainty = 0.0047; % MeV Hasterok & Gard 2017

summary(5).element = 'U';
summary(5).isotope = {'U-235','U-238'};
summary(5).abundance = [0.007204; 0.992742];
summary(5).abundance_uncertainty = [0.000006; 0.000010];
summary(5).daughter = {'Pb-207','Pb-206'};
summary(5).P_decay = 1;
summary(5).P_decay_uncertainty = 0;
summary(5).half_life = [7.04e8; 4.468e9]; % a
summary(5).half_life_uncertainty = [1e6; 6e6]; % a
summary(5).decay_energy = [44.2543; 46.596]; % MeV Hasterok & Gard 2017
summary(5).decay_uncertainty = [0.0038; 0.019]; % MeV Hasterok & Gard 2017
summary(5).decay_rate = [1.38068e-12; 2.2906e-13];
summary(5).decay_rate_uncertainty = [0.00073e-12; 0.0014e-13];
%summary(5).decay_energy = [44; 47.7]; % MeV Fiorentini et al., 2007
%summary(5).decay_energy = [45.26; 46.34]; % MeV Rybach, 1988
%summary(5).atomic_mass = [235.0439301; 238.0507884]; % NIST amu | g/mol
summary(5).atomic_mass = [235.0439282; 238.0507870]; % IAEA amu | g/mol
summary(5).atomic_mass_uncertainty = [1.2e-6; 1.6e-6]; % IAEA amu | g/mol
%summary(5).daughter_mass = [206.9758973; 205.9744657]; % NIST amu | g/mol
summary(5).daughter_mass = [206.9758967; 205.9744651]; % IAEA amu | g/mol
summary(5).number_alpha = [7; 8];
summary(5).neutrino_energy = [2.1424; 5.098]; % MeV Hasterok & Gard 2017
summary(5).neutrino_uncertainty = [0.0079; 0.017]; % MeV Hasterok & Gard 2017


Atot = 0;
for i = 1:length(summary)
    % compute decay constant (s^-1)
    summary(i).decay_constant = log(2)./summary(i).half_life / secperyr;
    summary(i).decay_constant_uncertainty = log(2) / secperyr * ...
        summary(i).half_life_uncertainty ./ summary(i).half_life.^2;
        summary(i).P_decay .* ...
        summary(i).decay_constant .* ...
        summary(i).decay_energy * 1e6 * Ce;
    
    summary(i).decay_rate = Na ./ summary(i).atomic_mass .* ...
        summary(i).decay_rate * 1e6 * Ce;
    summary(i).decay_rate = sum(summary(i).decay_rate,2);
     summary(i).decay_rate_uncertainty = 1e6 * ...
        sqrt( ...
        Na_uncertainty^2 * ...
            (Ce .* ...
            summary(i).decay_rate ./ summary(i).atomic_mass).^2 + ...
        Ce_uncertainty^2 * ...
            (Na * summary(i).decay_rate ./ ...
            summary(i).atomic_mass).^2 + ...
        summary(i).decay_rate_uncertainty.^2 .* ...
            (Na * ...
            Ce ./ summary(i).atomic_mass).^2 + ...
        summary(i).atomic_mass_uncertainty.^2 .* ...
            (Na * Ce .* ...
            summary(i).decay_rate ./ summary(i).atomic_mass.^2).^2 ...
        );
    
    summary(i).isotope_heat_uncertainty = 1e9 * ...
        sqrt( ...
        Na_uncertainty^2 * ...
            (summary(i).P_decay .* summary(i).decay_energy .* Ce .* ...
            summary(i).decay_constant ./ summary(i).atomic_mass).^2 + ...
        summary(i).P_decay_uncertainty.^2 .* ...
            (Na * summary(i).P_decay .* summary(i).decay_energy .* Ce .* ...
            summary(i).decay_constant ./ summary(i).atomic_mass).^2 + ...
        summary(i).decay_uncertainty.^2 .* ...
            (Na * summary(i).P_decay .* Ce .* ...
            summary(i).decay_constant ./ summary(i).atomic_mass).^2 + ...
        Ce_uncertainty^2 * ...
            (Na * summary(i).P_decay .* summary(i).decay_energy .* ...
            summary(i).decay_constant ./ summary(i).atomic_mass).^2 + ...
        summary(i).decay_constant_uncertainty.^2 .* ...
            (Na * summary(i).P_decay .* summary(i).decay_energy .* ...
            Ce ./ summary(i).atomic_mass).^2 + ...
        summary(i).atomic_mass_uncertainty.^2 .* ...
            (Na * summary(i).P_decay .* summary(i).decay_energy .* Ce .* ...
            summary(i).decay_constant ./ summary(i).atomic_mass.^2).^2 ...
        );
    
    summary(i).isotope_heat = sum(summary(i).isotope_heat,2);
    
    summary(i).isotope_heat_uncertainty = sqrt(sum(summary(i).isotope_heat_uncertainty.^2,2));
    % concentration of isotope in rock sample
    summary(i).concentration = element(i) * summary(i).abundance .* exp(summary(i).decay_constant * time);
    summary(i).concentration_uncertainty = element(i) .* summary(i).abundance_uncertainty .*  exp(summary(i).decay_constant * time);
    
    % total heat produced
    summary(i).heat = density * summary(i).concentration .* summary(i).isotope_heat;
    summary(i).heat_uncertainty = density * sqrt(sum(summary(i).concentration_uncertainty.^2 .* summary(i).isotope_heat.^2 + summary(i).concentration.^2 .* summary(i).isotope_heat_uncertainty.^2));
    
    summary(i).heat_rate = density * summary(i).concentration .* summary(i).decay_rate;
    summary(i).heat_rate_uncertainty = density * sqrt(sum(summary(i).concentration_uncertainty.^2 .* summary(i).decay_rate.^2 + summary(i).concentration.^2 .* summary(i).decay_rate_uncertainty.^2));
    summary(i).heat_rate_uncertainty = density * sqrt(sum(summary(i).concentration.^2 .* summary(i).decay_rate_uncertainty.^2));
    summary(i).mass_diff = amu2kg * ( ...
        (summary(i).atomic_mass - summary(i).daughter_mass) - ...
        summary(i).number_alpha * alpha_mass );
    
    summary(i).mass_energy = summary(i).mass_diff / eV2kg / 1e6 - summary(i).neutrino_energy;
    
    summary(i).isotope_heat2 = 1e3*Na./summary(i).atomic_mass .* ...
        summary(i).P_decay .* ...
        summary(i).decay_constant .* ...
        summary(i).mass_energy * 1e6 * Ce;
    summary(i).isotope_heat2 = sum(summary(i).isotope_heat2,2);
    
    summary(i).heat2 = density * summary(i).concentration .* summary(i).isotope_heat2;
    
    %summary(i)
    
    Atot = Atot + sum(summary(i).heat_rate);
end

return


% figure;
% subplot(231);
% hold on;
% plot(age,summary.Ck40 / (K * X(1)),'k');
% plot(age,summary.Cth232 / (Th * X(2)),'r');
% plot(age,summary.Cu235 / (U * X(3)),'b--');
% plot(age,summary.Cu238 / (U * X(4)),'b:');
% plot(age,summary.Cu / (U * X(3) + U * X(4)),'b-');
% title('Concentration of HPE''s');
% xlabel('Time before present [Ga]');
% ylabel('Cp/Cp_{t=0}');
% set(gca,'Box','on');
% ylim([0 10]);
% xlim([-4.55 4.55]);
% axis square;
% 
% subplot(232);
% hold on;
% plot(age,Ha(1,:)/(H(1) + H(2)),'k');
% plot(age,Ha(2,:)/H(3),'r--');
% plot(age,Ha(3,:)/(X(3)*H(4) + X(4)*H(5)),'b');
% title('Heat Production Per Unit Mass');
% xlabel('Time before present [Ga]');
% ylabel('H/H_{t=0}');
% set(gca,'Box','on');
% xlim([-4.55 4.55]);
% ylim([0 2.5]);
% axis square;
% 
% ind = find(t == 0);
% subplot(233);
% hold on;
% plot(age,Ak/Ak(ind),'k-');
% plot(age,Ath/Ath(ind),'r-');
% plot(age,Au/Au(ind),'b-');
% plot(age,(Atot)/(Atot(ind)),'g-');
% Atot/Atot(ind)
% title('Heat Production');
% xlabel('Time before present [Ga]');
% ylabel('A/A_{t=0}');
% set(gca,'Box','on');
% xlim([-4.55 4.55]);
% ylim([0 10]);
% axis square;
% 
% subplot(234);
% hold on;
% plot(age,Ak,'k-');
% plot(age,Ath,'r-');
% plot(age,Au,'b-');
% plot(age,Ak+Ath+Au,'g-');
% title('Heat Production');
% xlabel('Time before present [Ga]');
% ylabel('A [\muW/m^3]');
% set(gca,'Box','on');
% xlim([-4.55 4.55]);
% axis square;

return
