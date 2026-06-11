function h = elratio_vs_age(data,age_div,scale,el1,varargin)
% ELRATIO_VS_AGE - Produces element and element ratio plots with time
%
%   elratio_vs_age(data,age_div,scale,el) will produce a age vs element
%   concentration plot using geochemical analyses found in the data table.
%   The element will plot on a 'log' (default) or 'linear' scaled axis.
%   The data will be binned by the bin edges supplied by age_div.
%
%   To create a plot of element ratios, use
%   elratio_vs_age(data,age_div,scale,el1,el2)
%
%   For K, Rb, Sm, U, and Th, the concentrations can be corrected for
%   radiogenic decay.  To include a decay correction,
%   elratio_vs_age(data,age_div,scale,el1,...,decaycorr) where decaycorr is
%   a 1 = yes, 0 = no.

% Last Modified: 29 Oct. 2020
% D. Hasterok (dhasterok@gmail.com), University of Adelaide

% keep element names as text for labels
txt1 = el1;

if ~strcmp(data.Properties.VariableNames,el1)
    el1 = lower(el1);
    if ~strcmp(data.Properties.VariableNames,el1)
        el1 = [el1,'_ppm'];
    end
end

rflag = 0;          % ratio flag
decaycorr = 0;      % flag for decay correction
colour = [0 0 0];
if nargin >= 5
    opt = 1;
    while opt + 4 <= nargin
        switch lower(varargin{opt})
            case 'decaycorrect'
                decaycorr = 1;
                opt = opt + 1;
            case 'color'
                colour = varargin{opt+1};
                opt = opt + 2;
            otherwise
                if ~ischar(varargin{1})
                    error('Field must be a char array.')
                end
                el2 = varargin{1};
                txt2 = el2;
                el2 = lower(el2);
                if ~strcmp(data.Properties.VariableNames,el2)
                    el2 = [el2,'_ppm'];
                end

                % flag indicating take ratio of el1 / el2
                rflag = 1;
                opt = opt + 1;
        end
    end
end

% ensure data are positive and have associated ages
if rflag
    ind = data{:,el1} > 0 & data{:,el2} > 0 & ...
        ~isnan(data.avg_age);
else
    ind = data{:,el1} > 0 & ...
        ~isnan(data.avg_age);
end
data = data(ind,:);

% decay correction
radel = {'k_ppm','rb_ppm','sm_ppm','th_ppm','u_ppm'};
dref = 'R88'; % decay constants reference (Rybach, 1988)

if decaycorr
    % element 1 if necessary
    if any(strcmp(el1,radel))
        data{:,el1} = decaycorrect(el1,data{:,el1},data.avg_age,dref);
    end
    
    % element 2 if necessary
    if exist(el2)
        if any(strcmp(el2,radel))
            data{:,el2} = decaycorrect(el2,data{:,el2},data.avg_age(ind),dref);
        end
    end
end

% set age divisions
for i = 1:length(age_div)-1
    ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1));
    agebin.n(i) = length(ind);
    agebin.ind{i} = ind;
    avg_age{i} = data.avg_age(ind);
    if rflag
        ratio{i} = data{ind,el1}./data{ind,el2};
    else
        ratio{i} = data{ind,el1};
    end
end


if strcmp(scale,'linear')
    % linear-scale
    [agebin.Qage,agebin.Qratio] = whisker(avg_age, ...
        ratio,'NoPlot');

    h = sqwavefill(agebin.Qratio,agebin.Qage(:,3),age_div,colour);
    
    ylim([0 ceil(max(agebin.Qratio(:,5)))]);
else
    % log-scale
    [agebin.Qage,agebin.Qratio] = whisker(avg_age, ...
        ratio,'Scale','log','NoPlot');

    h = sqwavefill(agebin.Qratio,agebin.Qage(:,3),age_div,colour);

    hpax([floor(min(agebin.Qratio(:,1))) ceil(max(agebin.Qratio(:,5)))]);
end 
set(gca,'Box','on', ...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',8);

xlim([age_div(1) age_div(end)]);
xlabel('Age [Ga]');
set(gca,'XTick',[0:500:4000],'XTickLabel',[0:0.5:4]);

a = pettitt(agebin.Qratio(:,3));
aflip = pettitt(flipud(agebin.Qratio(:,3)));
age_div_flip = fliplr(age_div);
if rflag
    try
    ylabel([txt1,'/',txt2],'FontSize',10);
    fprintf('  %-10s age: %.0f Ma,  p-value: %.4f\n     flipped age: %.0f Ma,  p-value: %.4f\n',[txt1,'/',txt2],age_div(a(1)+1),a(3),age_div_flip(aflip(1)+1),aflip(3));
    catch
    end
else
    switch lower(txt1)
        case {'eu_anomaly','sr_anomaly','pb_anomaly','nd_anomaly'}
            txt1 = [upper(txt1(1)),lower(txt1(2)),'*'];
            units = '';
        case {'sio2','tio2','feo','mgo','mno','cao','na2o','k2o'}
            units = ' (wt.%)';
        otherwise
            units = ' (ppm)';
    end
    ylabel([txt1,units],'FontSize',10);
    try
        fprintf('  %-10s age: %.0f Ma,  p-value: %.4f,   flipped age: %.0f Ma,  p-value: %.4f\n',[txt1],age_div(a(1)+1),a(3),age_div_flip(aflip(1)+1),aflip(3));
    catch
        % ignore error
    end
end


return