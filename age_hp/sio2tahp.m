function hp_corrected = sio2tahp(data)
%Calculates the best fit plane in 3d space, hp vs sio2 vs ta

    ind = data.heat_production > 0 & data.heat_production < 50 &...
        data.sio2 > 30 & data.sio2 < 85 &...
        data.k2o > 0 & data.na2o > 0;
    data = data(ind,:);


    hp = data.heat_production;
    sio2 = data.sio2;
    TA = data.k2o + data.na2o;
    
    X = [sio2 TA];
    y = hp;
    b = robustfit(X,y);
    
    hp_corrected = b(1)+b(2).*sio2+b(3).*TA;
    
    
    scatter3(sio2,TA,hp,5,'r','fill')
    xlabel('sio2')
    ylabel('TA')
    zlabel('HP')
    set(gca,'xlim',[30 90],'ylim',[0 30],'zlim',[0 20])

    nBins = 100;
    xBins = linspace(min(sio2),max(sio2),nBins);
    yBins = linspace(min(TA),max(TA),nBins);
    zBins = linspace(min(hp),max(hp),nBins);
    D = zeros(nBins,nBins,nBins);
    
    for i = 1:numel(sio2)
        xi = find((sio2(i)>xBins),1,'last');
        yi = find((TA(i)>yBins),1,'last');
        zi = find((hp(i)>zBins),1,'last');
        D(xi,yi,zi) = D(xi,yi,zi) + 1;
    end
    sliceomatic(D)
    
    hp_corrected = 1;
return