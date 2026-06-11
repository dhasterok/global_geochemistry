%
% Sample code to demonstrate a uislider
% to update an ode solution
%

function TA_age(data,age_div)

    % Build the GUI
    init_ui;
    
    % This funciton sets up the figure and slider.
    % note the CreateFcn and the Callback.
    % The CreatFcn runs the first time while the Callback
    % runs each time you move the slider.
    function init_ui()
        %f1 = figure;
        slider1 = uicontrol('Style', 'slider',...
                            'Min',1,'Max',size(age_div,2)-1,'Value',1,...
                            'Position', [100 5 350 15],...
                            'CreateFcn', @solve_and_plot,...
                            'Callback',  @solve_and_plot,...
                            'SliderStep',[1/(size(age_div,2)-2) , 10/(size(age_div,2)-2)],...
                            'Units','normalized'); 
    end

    function agebin = agehpbox(data,age_div,name)
        %Tests for age column:
        %Set max age variance
        max_dage = 200;

        %Check for NaN age but age_min, age_max exist
        data.avg_age = data.age;
        ind = find(isnan(data.avg_age) & ~isnan(data.age_min) ...
            & ~isnan(data.age_max));
        dage = (data.age_max - data.age_min);
        dage(dage > max_dage) = NaN;
        data.avg_age(ind) = data.age_min(ind) + dage(ind)/2;
        for i = 1:length(age_div)-1
            ind = find(~isnan(data.na2o) & ~isnan(data.k2o) & ~isnan(data.sio2) & age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) ...
                & ~isnan(data.heat_production) & data.sio2 > 25 & (~strcmp(data.country,'unknown')|~strcmp(data.country,'ocean')));
            n(i) = length(ind);
            agebin.ind{i} = ind;
            avg_age{i} = data.avg_age(ind);
            TA{i} = data.na2o(ind) + data.k2o(ind);
            S{i} = data.sio2(ind);
        end
        agebin.TA = TA;
        agebin.S = S;
        agebin.avg_age = avg_age;
    end

    
    function plot_data(arrayname,age_ind, c)
        arrayname = arrayname(age_ind,:);
        TA = arrayname.na2o + arrayname.k2o;
        S = arrayname.sio2;
        
%         eS = [20:0.5:100];
%         eTA = [-6:0.2:30];
        eS = [20:1:100];
        eTA = [-6:0.5:30];

        n = hist2d(S,TA,eS,eTA);
        imagesc(eS,eTA,log10(n));
        colorbar;
        caxis(c);
        axis xy;
        xlabel('sio2 [wt.%]');
        ylabel('k2o + na2o [wt.%]');
        xlim([25 100]);
        ylim([-6 20]);
        set(gca,'Box','on');
        golden;
        
        hold on
        x = [0:100];
        y = 0.1776.*x -5.1874;
        plot(x,y,'-g')
        hold off
        
        hold on
        b = robustfit(S,TA);
        plot(x,b(1)+b(2).*x,'-r')
        hold off
    end




    % This funciton it what is run when you first run the 
    % code and each time you move the slider.
    function solve_and_plot(src,event)
        % Get the current slider value
        Tmax = get(src,'Value');
        Tmax = floor(Tmax);

        h3 = subplot(3,1,3);
        agebin = agehpbox(data,[0:200:Tmax*200],'Full data');
        xlim([age_div(1) age_div(end)])
        set(gca,'Box','on');
        pbaspect([2 1 1])
        ylim([-2 2])
        
                
        h2 = subplot(3,1,1:3);
        cla(h2)
        c = length(agebin.ind{Tmax})/(length(agebin.ind{Tmax})/2.5);
        plot_data(data,agebin.ind{Tmax},[-0.1 c]);
        str = strcat(num2str(age_div(Tmax)),'-',num2str(age_div(Tmax+1)),{' '},'Mya. Contains',{' '},num2str(length(agebin.ind{Tmax})),{' '},'data points');
        title(str)
  
    end

end