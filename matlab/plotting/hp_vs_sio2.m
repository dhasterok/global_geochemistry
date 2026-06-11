function [Qsio2,Qhp] = hp_vs_sio2(sio2,hp,sio2_div);

for i = 1:(length(sio2_div)-1)
    S{i} = sio2(sio2 >= sio2_div(i) & sio2 < sio2_div(i+1) ...
        & ~isnan(hp) & ~isnan(sio2),:);
    A{i} = hp(sio2 >= sio2_div(i) & sio2 < sio2_div(i+1) ...
        & ~isnan(hp) & ~isnan(sio2),:);
end

[Qsio2,Qhp] = whisker(S,A,'Color',[0.5 0.5 0.5],'Scale','log');
xlabel('SiO_2 (wt.%)');
hpax([-1 2]);
ylim([-1 log10(30)]);

return