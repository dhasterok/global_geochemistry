function mantle_temp(data,age_div)
%Ni(ppm)/Mg(ppm) x 10000 again of igneous rocks with SiO2 between 40 and 56
%wt%. This is effectively and indication of the extent of mantle melting
%and thus mantle temperature.

    %IGNORE IAB's
    %Calc DF1 and DF2
%     ind = data.tio2 >= 0 & data.sio2 >= 0 & data.al2o3 >=0 & data.feo_tot >= 0 ...
%         & ~isnan(data.fe_fe2_ratio) & data.mno >= 0 & data.mgo >= 0 & data.cao >= 0 ...
%         & data.na2o >= 0 & data.k2o >= 0 & data.p2o5 >= 0 ...
%         & (strcmpi(data.rock_type,'subalkalic basalt') | strcmpi(data.rock_type,'alkalic basalt'));
%             %& data.sio2 <= 52 & data.sio2 >=25;
%     
%     %ratio is feo/(feo+fe2o3)
%     data.fe2o3 = data.feo_tot.*(1-data.fe_fe2_ratio);
%     data.feo = data.feo_tot.*(data.fe_fe2_ratio);
%     
%     DF1 = -4.6761.*log(data.tio2./data.sio2)...
%         +2.5330.*log(data.al2o3./data.sio2)...
%         -0.3884.*log(data.fe2o3./data.sio2)...
%         +3.9688.*log(data.feo./data.sio2)...
%         +0.8980.*log(data.mno./data.sio2)...
%         -0.5832.*log(data.mgo./data.sio2)...
%         -0.2896.*log(data.cao./data.sio2)...
%         -0.2704.*log(data.na2o./data.sio2)...
%         +1.0810.*log(data.k2o./data.sio2)...
%         +0.1845.*log(data.p2o5./data.sio2)...
%         +1.5445;
%     
%     DF2 = 0.6751.*log(data.tio2./data.sio2)...
%         +4.5895.*log(data.al2o3./data.sio2)...
%         +2.0897.*log(data.fe2o3./data.sio2)...
%         +0.8514.*log(data.feo./data.sio2)...
%         -0.4334.*log(data.mno./data.sio2)...
%         +1.4832.*log(data.mgo./data.sio2)...
%         -2.3627.*log(data.cao./data.sio2)...
%         -1.6558.*log(data.na2o./data.sio2)...
%         +0.6757.*log(data.k2o./data.sio2)...
%         +0.4130.*log(data.p2o5./data.sio2)...
%         +13.1639;
%     
%     line_segs = [1.160,-0.333,5.912,8;%
%         -0.266,0.020,-4.190,8;
%         -8,-2.490,-0.266,0.020;
%         1.160,-0.333,3.431,-8;%
%         1.160,-0.333,-0.266,0.020];
%     x=[line_segs(:,1) line_segs(:,3)];
%     y=[line_segs(:,2) line_segs(:,4)];
    
    ind = (isnan(data.mg_ppm) | data.mg_ppm < 0) & data.mgo >= 0;
    data.mg_ppm(ind,:) = 10000.*(data.mgo(ind,:)./((24.305 + 15.994)./24.305));
    ind = data.sio2 >= 40 & data.sio2 <= 56 & data.ni_ppm >= 0 & data.mg_ppm >=0;
    data = data(ind,:);

    
    
    %Age div indices
    for i = 1:length(age_div)-1
        ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1));
        n(i) = length(ind);
        agebin.ind{i} = ind;
        avg_age{i} = data.avg_age(ind);
        ni_ppm{i} = data.ni_ppm(ind);
        mg_ppm{i} = data.mg_ppm(ind);
    end
    
    %Loop through ages and find ratio and magnitude ratio of TA to line
    
    number_above = [];
    number_below = [];
    ratio_above = [];
    
    for i = 1:length(agebin.ind)
        datapoints = data(agebin.ind{i},:);
        avg_age{i} = datapoints.avg_age;
        mantle_temp_median(i) = nanmedian((datapoints.ni_ppm./datapoints.mg_ppm).*10000);
        mantle_cycle(i) = nanmedian(datapoints.feo_tot./datapoints.mno);
    end
    
    
    figure()
    subplot(2,1,1)
    hold on
    plot((age_div(2:end)+age_div(1:end-1))./2,mantle_temp_median,'-ob')
    hold off
    xlim([0 4000])
    title('Ni/Mg ppm * 10000')
    ylabel('Ni/Mg ppm * 10000')
    xlabel('Age [Ma]')
    
    subplot(2,1,2)
    hold on
    plot((age_div(2:end)+age_div(1:end-1))./2,mantle_cycle,'-ob')
    hold off
    xlim([0 4000])
    title('Fe/Mn')
    ylabel('Fe/Mn')
    xlabel('Age [Ma]')

    
return