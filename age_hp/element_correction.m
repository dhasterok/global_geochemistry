function [data,varargout] = element_correction(data)
% data = element_correction(data)
% Negative values are switched to positive for u, th and k
% Major elements with > 100 wt% are divided by 10000 (ppm to wt% convert)
% u_ppm and th_ppm > 1000 are divided by 1000 (assumed ppb instead of ppm)
% Total of majors is also added: if greater than 110, NaN'd
% Check k2o and k_ppm are ~similar

% Add an output (optional) - ind of those removed completely.

fprintf('\n------------------\n')
fprintf('element_correction\n')
fprintf('------------------\n\n')


%--------------------------------------------------------------------------
% 1. Switching negatives
%--------------------------------------------------------------------------
ind = data.k_ppm < 0;
data.k_ppm(ind) = NaN;
fprintf('k_ppm < 0 removed: %d\n',length(find(ind)))

ind = data.u_ppm < 0;
data.u_ppm(ind) = NaN;
fprintf('u_ppm < 0 removed: %d\n',length(find(ind)))

ind = data.th_ppm < 0;
data.th_ppm(ind) = NaN;
fprintf('th_ppm < 0 removed: %d\n',length(find(ind)))

%--------------------------------------------------------------------------
% 2. Major element check
%--------------------------------------------------------------------------

major_elements = {'sio2','tio2','al2o3','cr2o3','feo_tot','mgo','cao',...
    'mno','nio','k2o','na2o','sro','p2o5'};
volatiles = {'h2o_tot','co2','so3','f_ppm','cl_ppm'};

sum_majors = zeros(size(data.sio2,1),1);
for i = 1:length(major_elements)
    ind = (data.(major_elements{i}) > 102);
    fprintf('%s > 102wt%% converted to ppm (not done anymore): %d\n',major_elements{i},length(find(ind)))
   
%     % Pause
%     fprintf('%s\n',major_elements{i})
%     data(ind,{'filename','sample_name',major_elements{1:end}})
%     x = input('\n');
    
    %data.(major_elements{i})(ind) = data.(major_elements{i})(ind)./10000;
    
    % See individual spreadsheets for each elemennt
%     data(ind,{'filename','sample_name','age_min','age_max','age',major_elements{:}})
%     x = input('\n');
    
    sum_majors = nansum([sum_majors data.(major_elements{i})],2);
end


% If the sum of majors is > 85, normalise to 100wt%
% Remove the rest?





ind = (sum_majors > 120);


% Pause
fprintf('Total greater than 120:\n')

totaltable = cell2table(num2cell(sum_majors),'VariableNames',{'Total'});

% [data(ind,{'filename','sample_name','age_min','age_max','age',major_elements{:},volatiles{:}}),totaltable(ind,:)]
% x = input('\n');


data(ind,:) = [];
removed_ind = ind;
fprintf('Major element sum > 120 removed: %d\n',length(find(ind)))

%--------------------------------------------------------------------------
% 3. HP element check
%--------------------------------------------------------------------------

% Uranium and thorium > 1000 ppm
ppm_threshhold = 500;

ind = data.u_ppm >= ppm_threshhold;
fprintf('u_ppm >= %d identified (not changed): %d\n',ppm_threshhold,length(find(ind)))

ind = data.th_ppm >= ppm_threshhold;
fprintf('th_ppm >= %d identified (not changed): %d\n',ppm_threshhold,length(find(ind)))


%--------------------------------------------------------------------------
% 4. Check consistency between k2o and k_ppm
%--------------------------------------------------------------------------

ind = ~isnan(data.k2o) & ~isnan(data.k_ppm);
atomic_k = 39.0983;
atomic_o = 15.9994;
conv_fact = (2*atomic_k)./(2*atomic_k + atomic_o);

total_k = length(find(ind));

ind = ((data.k_ppm > ((data.k2o.*conv_fact.*10000)).*1.25) |  ...
    (data.k_ppm < ((data.k2o.*conv_fact.*10000)).*0.75));

fprintf('No. samples with inconsistent k_ppm and k2o (>25%% variance): %d\n',...
    length(find(ind)))

ind = data.k2o <= 1 & ((data.k_ppm > ((data.k2o.*conv_fact.*10000)).*1.25) |  ...
    (data.k_ppm < ((data.k2o.*conv_fact.*10000)).*0.75));
data.k2o(ind) = (data.k_ppm(ind)./conv_fact)./10000;
fprintf('Inconsistent K2O and k_ppm replaced with k_ppm (K2O <= 1 wt%%): %d\n',...
    length(find(ind)))

ind = data.k2o > 1 & ((data.k_ppm > ((data.k2o.*conv_fact.*10000)).*1.25) |  ...
    (data.k_ppm < ((data.k2o.*conv_fact.*10000)).*0.75));
data.k2o(ind) = (data.k_ppm(ind)./conv_fact)./10000;
fprintf('Removed k_ppm for inconsistent K2O > 1: %d\n',...
    length(find(ind)))
data.k_ppm(ind) = NaN;

%Convert k_ppm to k2o if not listed
ind = (data.k_ppm >= 0 & isnan(data.k2o));
data.k2o(ind,:) = ((2*atomic_k + atomic_o)/(2*atomic_k*1e4))*data.k_ppm(ind,:);
fprintf('k_ppm converted to missing K2O value: %d\n',length(find(ind)))


% Return removed indices:
if nargout == 2
    varargout{1} = removed_ind;
end

return





