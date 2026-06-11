function harker(data,varargin)
% HARKER - Harker geochemical plots
%
%   harker(data) will produce Harker diagrams.
%
%   Options:
%       'PlotType'          select from 'scatter' (default) and 2D
%                           'histogram' (a.k.a. heatmap)
%
%       'Color'             provide color triplet as a 1x3 for all points,
%                           or height(data)x1 for each individual point or
%                           as a height(data)x3 color triplet for each
%                           point (scatter only)
%
%       'Symbol'            set symbol (see plot, scatter only)
%
%       'MarkerFaceAlpha'   set point transparency (scatter only), ranging
%                           from 0 for transparent to 1 for opaque
%
%       'NBins'             set number of bins for histogram (default, 25)
%
%       'BinScale'          set bin color scale to 'linear' or 'log' (default)

% D. Hasterok, University of Adelaide
% Last Modified: 5 May 2023

p = inputParser;

addRequired(p,'data');
addParameter(p,'PlotType','scatter',@ischar);
addParameter(p,'Color',[0.6 0.6 0.6],@isnumeric);
addParameter(p,'Symbol','o',@ischar);
addParameter(p,'MarkerFaceAlpha',1,@isnumeric);
addParameter(p,'NBins',25,@isnumeric);
addParameter(p,'BinScale','log',@isnumeric)

parse(p,data,varargin{:});

type = p.Results.PlotType;
C = p.Results.Color;
sym = p.Results.Symbol;
alpha = p.Results.MarkerFaceAlpha;
switch p.Results.BinScale
    case 'linear'
        logflag = false;
    case 'log'
        logflag = true;
    otherwise
        error('Unknown BinScale.');
end

sax = [floor(min(data.sio2)/5)*5 ceil(max(data.sio2)/5)*5];
ox = {'tio2','al2o3','feo_tot','mgo','cao','mno','na2o','k2o','p2o5'};
yl = {'TiO_2','Al_2O_3','FeO_T','MgO','CaO','MnO','Na_2O','K_2O','P_2O_5'};

switch type
    case 'histogram'
        % sio2 bins
        sbin = [sax(1):(sax(2)-sax(1))/25:sax(2)];
        for i = 1:length(ox)
            % determin oxide bins
            ymax = ceil(1.3*quantile(data.(ox{i}),0.95));
            yax = [0 ymax];

            dy = [yax(2) - yax(1)]/25;

            % compute oxide heatmap
            n = hist3([data.sio2,data.(ox{i})],{sbin [yax(1)+dy/2:dy:yax(2)+dy/2]});

            subplot(3,3,i);
            
            % log or linear scale?
            if logflag
                imagesc(sbin,[yax(1)+dy/2:dy:yax(2)-dy/2],log10(n(:,1:end-1))');
                clim([-0.1 ceil(max(log10(n),[],'all'))]);
                cb = colorbar;
                cb.Label.Text = 'log_{10}(N)';
            else
                imagesc(sbin,[yax(1)+dy/2:dy:yax(2)-dy/2],n(:,1:end-1)');
                cb = colorbar;
                cb.Label.String = 'N';
            end
            hold on;
            axis xy;
            % the fill works fine on its own, but throws a warning and will
            % not plot here... not sure why.
            %fill([0 100 100 0],[100 0 100 100],0.7*[1 1 1],'EdgeColor','none');
            xlabel('SiO_2 (wt.%)');
            ylabel([yl{i},' (wt.%)']);
            xlim(sax);
            ylim(yax);
            set(gca,'Box','on','TickDir','out');
            shading flat;
        end
    case 'scatter'
        for i = 1:length(ox)
            subplot(3,3,i); hold on;
            % no data can plot here
            fill([0 100 100 0],[100 0 100 100],[0.7 0.7 0.7],'EdgeColor','none');
        
            % plot data
            scatter(data.sio2,data.(ox{i}),[],C,sym,'filled','MarkerFaceAlpha',alpha);

            if height(C) == height(data)
                colorbar;
            end
        
            xlabel('SiO_2 (wt.%)');
            ylabel([yl{i},' (wt.%)']);
            xlim(sax);
            ymax = ceil(1.3*quantile(data.(ox{i}),0.95));
            ylim([0 ymax]);
            set(gca,'Box','on');
        end
    otherwise
        error('Unknown PlotType.');
end

return