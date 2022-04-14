
eref = loaderef;

sym = {'^','d','v'};

elayer = unique(eref.layer);
figure;
for i = 1:length(elayer)
    ind = strcmp(elayer{i},eref.layer);

    subplot(121);
    hold on;
    plot(eref.p_velocity(ind),log10(eref.heat_production(ind)),sym{i});

    subplot(122);
    hold on;
    plot(eref.density_bk(ind),log10(eref.heat_production(ind)),sym{i});
end

ind = strcmp('Rudnick and Gao (2003)',eref.reference);
subplot(121);
hold on;
plot(eref.p_velocity(ind),log10(eref.heat_production(ind)),'p');
xlabel('Estimated P-Velocity [km s^{-1}]');
ylabel('Heat Production [\muW m^{-3}]');
hpax([-2 1]);
xlim([5.8 8.4]);
golden;

subplot(122);
hold on;
plot(eref.density_bk(ind),log10(eref.heat_production(ind)),'p');
xlabel('Estimated Density [kg m^{-3}]');
ylabel('Heat Production [\muW m^{-3}]');
hpax([-2 1]);
xlim([2500 3500]);
golden;


