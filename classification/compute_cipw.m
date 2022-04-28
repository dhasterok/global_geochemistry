function data = compute_cipw(data)

for i = 1:length(data.ROCK_GROUP);
    oxides(1) = data.SIO2(i);
    oxides(2) = data.TIO2(i);
    oxides(3) = data.AL2O3(i);
    oxides(4) = data.FE2O3(i);
    oxides(5) = data.FEO(i);
    oxides(6) = data.MNO(i);
    oxides(7) = data.MGO(i);
    oxides(8) = data.CAO(i);
    oxides(9) = data.NA2O(i);
    oxides(10) = data.K2O(i);
    oxides(11) = data.P2O5(i);
    oxides(12) = data.CO2(i);
    oxides(13) = data.SO3(i);
    oxides(14) = data.S_PPM(i)/10000;
    oxides(15) = data.F_PPM(i)/10000;
    oxides(16) = data.CL_PPM(i)/10000;
    oxides(17) = data.SR_PPM(i)/10000;
    oxides(18) = data.BA_PPM(i)/10000;
    oxides(19) = data.NI_PPM(i)/10000;
    oxides(20) = data.CR_PPM(i)/10000;
    oxides(21) = data.ZR_PPM(i)/10000;

    tmp = cipw(oxides);
    data.minerals(i) = tmp;
end

return
