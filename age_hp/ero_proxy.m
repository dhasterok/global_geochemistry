function [mafic_ind,felsic_ind] = ero_proxy(data,age_div,div_value)
    figure()
    for i = 1:length(age_div)-1
        mafic_ind{i} = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) ...
            & data.sio2<=div_value & data.sio2>0);
    
        felsic_ind{i} = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) ...
            & data.sio2>div_value & data.sio2<=100);

        felsic_num(i) = length(felsic_ind{i});
        mafic_num(i) = length(mafic_ind{i});
    end
    
    h = bar(age_div(1:end-1),(felsic_num-mafic_num)./(felsic_num+mafic_num),'histc');
    set(h,'FaceColor',[0.5 0.5 0.5])
    ylabel('Composition bias [+ve Felsic, -ve Mafic]');
    xlabel('Age [Ma]');
    title('Felsic vs. Mafic bias');
    set(gca,'Box','on')
return