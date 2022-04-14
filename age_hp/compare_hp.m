function stats = compare_hp(data);

ind = data.u_ppm > 0 & ...
    data.th_ppm > 0 & ...
    data.k2o > 0 & ...
    data.density_model > 0 & ...
    data.rb_ppm > 0 & ...
    data.sm_ppm > 0;

sum((ind & rockgroup(data,'all igneous')))/sum(ind)

data = data(ind,:);

[A,A_iso] = radtime(data.density_model, ...
    data.k2o, ...
    data.rb_ppm, ...
    data.sm_ppm, ...
    data.th_ppm, ...
    data.u_ppm, ...
    'k2o','age',0,'formula','hg17');

figure;
subplot(121);
histogram(A_iso(:,:,2),'BinWidth',0.001);
xlim([0 0.1]);
golden

subplot(122);
histogram(100*reshape(A_iso(:,:,2),size(A))./A,'BinWidth',0.001);
xlim([0 0.3]);
golden
per = 100*reshape(A_iso(:,:,2),size(A))./A;
[mean(per) std(per)]
length(per)
sum(per > 1)
sum(per > 10)

age = [0:1:4]*1e3;

Ahg17 = radtime(data.density_model, ...
    data.k2o, ...
    data.rb_ppm, ...
    data.sm_ppm, ...
    data.th_ppm, ...
    data.u_ppm, ...
    'k2o','age',age,'formula','hg17');

Ahg17_norb = radtime(data.density_model, ...
    data.k2o, ...
    zeros(size(data.k2o)), ...
    zeros(size(data.k2o)), ...
    data.th_ppm, ...
    data.u_ppm, ...
    'k2o','age',age,'formula','hg17');
Ar88 = radtime(data.density_model, ...
    data.k2o, ...
    zeros(size(data.k2o)), ...
    zeros(size(data.k2o)), ...
    data.th_ppm, ...
    data.u_ppm, ...
    'k2o','age',age,'formula','r88');
Ats14 = radtime(data.density_model, ...
    data.k2o, ...
    zeros(size(data.k2o)), ...
    zeros(size(data.k2o)), ...
    data.th_ppm, ...
    data.u_ppm, ...
    'k2o','age',age,'formula','ts14');
Ad12 = radtime(data.density_model, ...
    data.k2o, ...
    zeros(size(data.k2o)), ...
    zeros(size(data.k2o)), ...
    data.th_ppm, ...
    data.u_ppm, ...
    'k2o','age',age,'formula','d12');

figure(2);
xax = [-7 4];
for i = 1:length(age)
    subplot(5,1,i);
    dr88 = (Ahg17(:,i)-Ar88(:,i))./Ar88(:,i)*100;
    dr88_norb = (Ahg17_norb(:,i)-Ar88(:,i))./Ar88(:,i)*100;
    dts14 = (Ahg17(:,i)-Ats14(:,i))./Ats14(:,i)*100;
    dts14_norb = (Ahg17_norb(:,i)-Ats14(:,i))./Ats14(:,i)*100;
    dd12 = (Ahg17(:,i)-Ad12(:,i))./Ad12(:,i)*100;
    dd12_norb = (Ahg17_norb(:,i)-Ad12(:,i))./Ad12(:,i)*100;
    
    stats.dr88(i,:) = [mean(dr88) std(dr88) quantile(dr88,[0.01 0.99]) min(dr88) max(dr88)];
    stats.dr88_norb(i,:) = [mean(dr88_norb) std(dr88_norb) quantile(dr88_norb,[0.01 0.99]) min(dr88_norb) max(dr88_norb)];
    stats.dts14(i,:) = [mean(dts14) std(dts14) quantile(dts14,[0.01 0.99]) min(dts14) max(dts14)];
    stats.dts14_norb(i,:) = [mean(dts14_norb) std(dts14_norb) quantile(dts14_norb,[0.01 0.99]) min(dts14_norb) max(dts14_norb)];
    stats.dd12(i,:) = [mean(dd12) std(dd12) quantile(dd12,[0.01 0.99]) min(dd12) max(dd12)];
    
    histogram(dr88,'BinWidth',0.05,'DisplayStyle','stairs');
    hold on;
    histogram(dr88_norb,'BinWidth',0.05,'DisplayStyle','stairs');
    histogram(dts14,'BinWidth',0.05,'DisplayStyle','stairs');
    histogram(dts14_norb,'BinWidth',0.05,'DisplayStyle','stairs');
    histogram(dd12,'BinWidth',0.05,'DisplayStyle','stairs');
    histogram(dd12_norb,'BinWidth',0.05,'DisplayStyle','stairs');
    xlim(xax);
    ylim([0 10000]);
    
    set(gca,'XTick',[xax(1):xax(2)]);
    
    if i == 1
        legend('R88','R88, no Rb','TS14','TS14 no Rb','D12','D12 no Rb');
    end
    
    text(-6,7000,[num2str(age(i)),' Ga']);
    
    if i ~= length(age)
        set(gca,'XTickLabel',{});
    end
    
    if i == length(age)
        xlabel('Difference (%) HG17 to R88 or TS14');
    end
end

ind = rockgroup(data,'all igneous') & strcmp(data.rock_type,'subalkalic basalt');

[median(data.density_model(ind)), ...
    10.^median(log10(data.k2o(ind))), ...
    10.^median(log10(data.rb_ppm(ind))), ...
    10.^median(log10(data.sm_ppm(ind))), ...
    10.^median(log10(data.th_ppm(ind))), ...
    10.^median(log10(data.u_ppm(ind)))]

age = [0:5:4000];
[A,A_iso] = radtime(median(data.density_model(ind)), ...
    10.^median(log10(data.k2o(ind))), ...
    10.^median(log10(data.rb_ppm(ind))), ...
    10.^median(log10(data.sm_ppm(ind))), ...
    10.^median(log10(data.th_ppm(ind))), ...
    10.^median(log10(data.u_ppm(ind))), ...
    'k2o','age',age,'formula','hg17');

A(1:1000/5:end)
s = size(A_iso);
A_iso = reshape(A_iso,[s(2) s(3)])';

figure;
subplot(121);
plot(age/1000,A_iso);
hold on;
plot(age/1000,A);
legend('^{40}K','^{87}Rb','^{232}Th','^{235}U','^{238}U','A Total');
xlabel('Time before present (Ga)');
ylabel('Heat Production [\muW m^{-2}]');
golden;
ylim([0 4]);

subplot(122);
plot(age/1000,A_iso./repmat(A,[5 1])*100);
xlabel('Time before present (Ga)');
ylabel('Contribution to Heat Production (%)');
golden;

return
