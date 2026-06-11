function [TZr,Pqtz,Pfsp] = pt_plot(data,age_div)

[data.TZr,data.TZr_sd] = zr_sat(data);
data.TZr = data.TZr - 273.15;
[data.Pqtz,data.Pfsp] = pressest(data);

figure;
ind = ~isnan(data.TZr) & ~isnan(data.Pqtz);
eT = [500:20:900];
eP = [50:10:500];
n = hist2d(data.TZr(ind),data.Pqtz(ind),eT,eP);
imagesc(eT,eP,n);
xlabel('Zircon Sat. Temperature [^\circC]');
ylabel('Pressure (Qtz) [MPa]');

figure;
subplot(311);
TZr = plot_p_or_t(data,'TZr',age_div,'T');
ylabel('Zircon Sat. Temperature [^\circC]');

subplot(312);
Pqtz = plot_p_or_t(data,'Pqtz',age_div,'P');
ylabel('Pressure (Qz) [MPa]');
ylim([0 400]);

subplot(313);
Pfsp = plot_p_or_t(data,'Pfsp',age_div,'P');
ylabel('Pressure (Ab + Or) [MPa]');
xlabel('Age [Ga]');
ylim([0 400]);

figure;
subplot(121);
plot(TZr.val(:,3),Pqtz.val(:,3),'o-');
for i = 1:length(age_div) - 1
    text(TZr.val(i,3)+1,Pqtz.val(i,3),num2str((age_div(i+1) + age_div(i))/2));
end
xlabel('Zircon Sat. Temperature [^\circC]');
ylabel('Pressure (qtz) [MPa]');

t.Zircon_sat_temp_degC = TZr.val;
t.Pressure_quartz_MPa = Pqtz.val;

t = struct2table(t);
writetable(t,'Granite_PT_estimates.csv');

subplot(122);
plot(TZr.val(:,3),Pfsp.val(:,3),'o-');
for i = 1:length(age_div) - 1
    text(TZr.val(i,3)+1,Pfsp.val(i,3),num2str((age_div(i+1) + age_div(i))/2));
end
xlabel('Zircon Sat. Temperature [^\circC]');
ylabel('Pressure (ab+or) [MPa]');

return

function [out,agebin] = plot_p_or_t(data,field,age_div,type)

for i = 1:length(age_div)-1
    ind = age_div(i) <= data.avg_age & ...
        data.avg_age < age_div(i+1) & ...
        data{:,field} > 0;
    agebin.n(i) = sum(ind);
    agebin.ind{i} = ind;
    agebin.avg_age{i} = data.avg_age(ind);
    agebin.val{i} = data{ind,field};
end

% heat production versus age
%---------------------------
% linear scale
if strcmp(type,'T')
    [out.age,out.val] = whisker(agebin.avg_age,agebin.val,'Color',[0.5 0.5 0.5]);
    cla;
    sqwavefill(out.val,out.age(:,3),age_div,[0,0,0]);
    ylim([500 900]);
elseif strcmp(type,'P')
    [out.age,out.val] = whisker(agebin.avg_age,agebin.val,'Color',[0.5 0.5 0.5]);
    cla;
    sqwavefill(out.val-3,out.age(:,3),age_div,[0,0,0]);
    %hpax([-2 1]);
    ylim([100 1000]);
end

%Color = [0,0,0]+alpha for greys, [0,0,0] is black, alpha up to 1


%a1 = pettitt(log10(agebin.Qhp(:,3)));
%a1

%if ~isempty(a1)
%    a2 = pettitt(log10(agebin.Qhp(a1(1)+1:end,3)));
%    a2
%end

xlim([age_div(1) age_div(end)]);

set(gca,'Box','on');
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);

set(gca,'XTick',[0:500:4000],'XTickLabel',[0:0.5:4]);

return