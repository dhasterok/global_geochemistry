%
% Sample code to demonstrate a uislider
% to update an ode solution
%

function geo_age(data,age_div)
    figure()
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
            ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) ...
                & ~isnan(data.heat_production) & data.sio2 > 25 & (~strcmp(data.country,'unknown')|~strcmp(data.country,'ocean')));
            n(i) = length(ind);
            agebin.ind{i} = ind;
            avg_age{i} = data.avg_age(ind);
            heat_production{i} = data.heat_production(ind);
        end
        [agebin.Qage,agebin.Qhp] = whisker(avg_age,heat_production,'Color',[0.5 0.5 0.5],'Scale','log');
        ylabel('Heat Production [\muW m^{-3}]');
        title(name);
        set(gca,'Box','on');
    end

    
    function plot_data(arrayname,age_ind)
        arrayname = arrayname(age_ind,:);
        
        id = arrayname.sample_id;
        lat_array = arrayname.latitude;
        lon_array = arrayname.longitude;
        cntry = arrayname.country;

        % check locations
        ind = find(lat_array < -90 | lat_array > 90 ...
            | lon_array < -180 | lon_array > 180);
        if ~isempty(ind)
            warning('LATITUDE or LONGITUDE out of bounds.');
        end
        ind_ocean = find(ismember(cntry, {'ocean','unknown'}));
        ind_country = find(~ismember(cntry, {'ocean','unknown'}));
        plot(lon_array(ind_ocean),lat_array(ind_ocean),'b.');
        plot(lon_array(ind_country),lat_array(ind_country),'g.');
    end




    % This function it what is run when you first run the 
    % code and each time you move the slider.
    function solve_and_plot(src,event)
        % Get the current slider value
        Tmax = get(src,'Value');
        Tmax = floor(Tmax);

        h3 = subplot(3,1,3);
        cla(h3)
        agebin = agehpbox(data,[0:200:Tmax*200],'Full data');
        xlim([age_div(1) age_div(end)])
        set(gca,'Box','on');
        pbaspect([2 1 1])
        ylim([-2 2])
        
                
        h2 = subplot(3,1,1:2);
        cla(h2)
        hold on
        plotcoast;
        plot_data(data,agebin.ind{Tmax});
        str = strcat(num2str(age_div(Tmax)),'-',num2str(age_div(Tmax+1)),{' '},'Mya. Contains',{' '},num2str(length(agebin.ind{Tmax})),{' '},'data points');
        title(str)
        hold off
  
    end

end
