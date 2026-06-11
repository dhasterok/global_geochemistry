%--------------------------------------------------------------------------
%                           Age adjustment
%--------------------------------------------------------------------------

%If no range is given, enter the age value into all of age_min/age/age_max

function [mean_age,weighted_avg,no_useful_data] = age_overlap(age,age_min,age_max,range_min,range_max,hp)
age_x = range_min:1:range_max;
weight = zeros(1,length(age));

for i=1:length(age)
    if age_min(i) == age(i) || age_max(i) == age(i) || age_min(i) == age_max(i)
        if age(i)>=range_min && age(i)<=range_max
            weight(i) = 1;
        else
            weight(i) = 0;
        end
    elseif age_min(i) > age(i)
        temp = age(i);
        age(i) = age_min(i);
        age_min(i) = temp;
        
    elseif age(i) > age_max(i)
        temp = age_max(i);
        age_max(i) = age(i);
        age(i) = temp;
        
    elseif floor(age_max(i)-age(i))==floor(age(i)-age_min(i))
        %Symmetric gaussian curve, assume age limits are 95.45% confidence\

        % Define normal distribution and integrate over range interval
        %normFun = @(x) 1/(sqrt(2*pi)*s)*exp(-(x-age(i)).^2./(2*s^2)); 
        %weight(i) = integral(normFun,range_min,range_max);
        
        %Also can auto produce plots using normspec
        %weight(i) = normspec([range_min,range_max],m,s)
        
        
        %Mean and std dev.
        m = age(i);
        s = (age(i)-age_min(i))./2;
        %s = age(i)-age_min(i);
        
        %Gaussian function, integrate under curve for range desired
        normFun = @(x)normpdf(x,m,s);
        weight(i) = integral(normFun,range_min,range_max);
        if weight(i) < 0.05
            weight(i) = 0;
        end
        
        %Optional plots to confirm working
%         hold on
%         plot(age_x,normpdf(age_x,m,s))
%         hold off
    else
        %Something to do with asymmetric distribution
        weight(i) = 0; %temp for now
    end
    
    %Calculate weighted averages
    weighted_avg = (sum(weight.*hp))./(sum(weight));
    if isnan(weighted_avg)
        weighted_avg = 0;
    end
    
    no_useful_data = length(weight(weight>0.05));
    mean_age = (range_max+range_min)./2;
end