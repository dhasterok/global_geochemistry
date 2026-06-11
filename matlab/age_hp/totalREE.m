function totalREE(data)

    %REE indices
    agebin.ind_REE = find(~isnan(data.heat_production)...
        & ~isnan(data.sc_ppm) & ~isnan(data.y_ppm) & ~isnan(data.la_ppm)...
        & ~isnan(data.ce_ppm) & ~isnan(data.pr_ppm) & ~isnan(data.nd_ppm)...
        & ~isnan(data.sm_ppm) & ~isnan(data.eu_ppm) & ~isnan(data.gd_ppm)...
        & ~isnan(data.tb_ppm) & ~isnan(data.dy_ppm) & ~isnan(data.ho_ppm)...
        & ~isnan(data.er_ppm) & ~isnan(data.tm_ppm) & ~isnan(data.yb_ppm)...
        & ~isnan(data.lu_ppm));
    
    t_REE = [abs(data.sc_ppm(agebin.ind_REE)) abs(data.y_ppm(agebin.ind_REE)) ...
        abs(data.la_ppm(agebin.ind_REE)) abs(data.ce_ppm(agebin.ind_REE)) ...
        abs(data.pr_ppm(agebin.ind_REE)) abs(data.nd_ppm(agebin.ind_REE)) ...
        abs(data.sm_ppm(agebin.ind_REE)) abs(data.eu_ppm(agebin.ind_REE)) ...
        abs(data.gd_ppm(agebin.ind_REE)) abs(data.tb_ppm(agebin.ind_REE)) ...
        abs(data.dy_ppm(agebin.ind_REE)) abs(data.ho_ppm(agebin.ind_REE)) ...
        abs(data.er_ppm(agebin.ind_REE)) abs(data.tm_ppm(agebin.ind_REE)) ...
        abs(data.yb_ppm(agebin.ind_REE)) abs(data.lu_ppm(agebin.ind_REE)) ];

    REE = sum(t_REE,2);
    hp_REE = data.heat_production(agebin.ind_REE);

%Plot in age or HP?
%[agebin.Qage_REE,agebin.QREE] = whisker(avg_age_REE,REE,'Color',[0.5 0.5 0.5],'Scale','log');


filt_granit = (~cellfun(@isempty,regexp(data.rock_type(agebin.ind_REE),'granit','ignorecase')));
filt_diorit = (~cellfun(@isempty,regexp(data.rock_type(agebin.ind_REE),'diorit','ignorecase')));
filt_gabbro = (~cellfun(@isempty,regexp(data.rock_type(agebin.ind_REE),'gabbro','ignorecase')));

% filt_granit = (~cellfun(@isempty,regexp(data.rock_name(agebin.ind_REE),'granit','ignorecase')));
% filt_diorit = (~cellfun(@isempty,regexp(data.rock_name(agebin.ind_REE),'diorit','ignorecase')));
% filt_gabbro = (~cellfun(@isempty,regexp(data.rock_name(agebin.ind_REE),'gabbro','ignorecase')));
% 
% filt_granit = (~cellfun(@isempty,regexp(data.rock_name(agebin.ind_REE),'granit','ignorecase')) | ~cellfun(@isempty,regexp(data.rock_type(agebin.ind_REE),'granit','ignorecase')));
% filt_diorit = (~cellfun(@isempty,regexp(data.rock_name(agebin.ind_REE),'diorit','ignorecase')) | ~cellfun(@isempty,regexp(data.rock_type(agebin.ind_REE),'diorit','ignorecase')));
% filt_gabbro = (~cellfun(@isempty,regexp(data.rock_name(agebin.ind_REE),'gabbro','ignorecase')) | ~cellfun(@isempty,regexp(data.rock_type(agebin.ind_REE),'gabbro','ignorecase')));


figure()
%REE
subplot(4,3,3)
plot(REE,hp_REE,'b.')
xlabel('Total REE (ppm)');
ylabel('Heat production');
str = sprintf('Total REE: %f',length(hp_REE));
title(str);
set(gca,'Box','on');
xlim([0 800])
ylim([0 10])
%pbaspect([1.618 1 1])

subplot(4,3,6)
plot(REE(filt_gabbro,:),hp_REE(filt_gabbro,:),'b.')
xlabel('Total REE (ppm)');
ylabel('Heat production');
str = sprintf('Total REE - Gabbroic: %f',length(REE(filt_gabbro,:)));
title(str);
set(gca,'Box','on');
xlim([0 800])
ylim([0 10])
%pbaspect([1.618 1 1])

subplot(4,3,9)
plot(REE(filt_diorit,:),hp_REE(filt_diorit,:),'b.')
xlabel('Total REE (ppm)');
ylabel('Heat production');

str = sprintf('Total REE - Dioritic: %f',length(REE(filt_diorit,:)));
title(str);

set(gca,'Box','on');
xlim([0 800])
ylim([0 10])
%pbaspect([1.618 1 1])

subplot(4,3,12)
plot(REE(filt_granit,:),hp_REE(filt_granit,:),'b.')
xlabel('Total REE (ppm)');
ylabel('Heat production');
title('Total REE - Granitic');

str = sprintf('Total REE - Dioritic: %f',length(REE(filt_granit,:)));
title(str);

set(gca,'Box','on');
xlim([0 800])
ylim([0 10])
%pbaspect([1.618 1 1])

%sio2
sio2 = data.sio2(agebin.ind_REE);
subplot(4,3,1)
plot(sio2,hp_REE,'b.')
xlabel('sio2');
ylabel('Heat production');

title('sio2');

set(gca,'Box','on');
xlim([40 100])
ylim([0 10])
%pbaspect([1.618 1 1])
strmin = ['n = ',num2str(length(REE))];
text(40+0.0333*60,10 - 0.1*10,strmin,'HorizontalAlignment','left');

subplot(4,3,4)
plot(sio2(filt_gabbro,:),hp_REE(filt_gabbro,:),'b.')
xlabel('sio2');
ylabel('Heat production');
title('Sio2 - Gabbroic');
set(gca,'Box','on');
xlim([40 100])
ylim([0 10])
%pbaspect([1.618 1 1])
strmin = ['n = ',num2str(length(REE(filt_gabbro,:)))];
text(40+0.0333*60,10 - 0.1*10,strmin,'HorizontalAlignment','left');

subplot(4,3,7)
plot(sio2(filt_diorit,:),hp_REE(filt_diorit,:),'b.')
xlabel('sio2');
ylabel('Heat production');
title('Sio2 - Dioritic');
set(gca,'Box','on');
xlim([40 100])
ylim([0 10])
%pbaspect([1.618 1 1])
strmin = ['n = ',num2str(length(REE(filt_diorit,:)))];
text(40+0.0333*60,10 - 0.1*10,strmin,'HorizontalAlignment','left');

subplot(4,3,10)
plot(sio2(filt_granit,:),hp_REE(filt_granit,:),'b.')
xlabel('sio2');
ylabel('Heat production');
title('Sio2 - Granitic');
set(gca,'Box','on');
xlim([40 100])
ylim([0 10])
%pbaspect([1.618 1 1])
strmin = ['n = ',num2str(length(REE(filt_granit,:)))];
text(40+0.0333*60,10 - 0.1*10,strmin,'HorizontalAlignment','left');

%feo_tot
feo_tot = data.feo_tot(agebin.ind_REE);
subplot(4,3,2)
plot(feo_tot,hp_REE,'b.')
xlabel('feo_tot');
ylabel('Heat production');
title('feo_tot');
set(gca,'Box','on');
xlim([0 25])
ylim([0 10])
%pbaspect([1.618 1 1])

subplot(4,3,5)
plot(feo_tot(filt_gabbro,:),hp_REE(filt_gabbro,:),'b.')
xlabel('feo_tot');
ylabel('Heat production');
title('feo_tot - Gabbroic');
set(gca,'Box','on');
xlim([0 25])
ylim([0 10])
%pbaspect([1.618 1 1])

subplot(4,3,8)
plot(feo_tot(filt_diorit,:),hp_REE(filt_diorit,:),'b.')
xlabel('feo_tot');
ylabel('Heat production');
title('feo_tot - dioritic');
set(gca,'Box','on');
xlim([0 25])
ylim([0 10])
%pbaspect([1.618 1 1])

subplot(4,3,11)
plot(feo_tot(filt_granit,:),hp_REE(filt_granit,:),'b.')
xlabel('feo_tot');
ylabel('Heat production');
title('feo_tot - granitic');
set(gca,'Box','on');
xlim([0 25])
ylim([0 10])
%pbaspect([1.618 1 1])


return