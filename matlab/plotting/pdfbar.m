function varargout = pdfbar(v,p,varargin)
% PDFBAR - plots histograms or PDF data as a strip.
%
%   pdfbar(v,p) will produce a plot of PDFs as a stacked set of strips for
%   comparison, where p(v) is the probability.  The ordinate, v, is common
%   to each pdf and p can be a matrix of [N length(v)] where N is the
%   number of PDFs to display.
%
%   h = pdfbar(v,p) will return graphics handles for each plot in h.
%
%   Additional options can be supplied as option, value pairs, e.g.,
%   pdfbar(v,p,'Option',value).
%
%   Options:
%       'YTickLabel'        Each PDF labeled
%       'CMap'              Colormap name or set of color triplets.  Uses
%                           extended colormap2 colorset.
opt = 1;
cmap = 'red';
w = [];
plotmax = false(1);
if nargin > 2
    while opt + 2 <= nargin
        switch lower(varargin{opt})
            case 'width'
                w = varargin{opt+1};
                w = 100*[0; cumsum(w)]/sum(w);
                %figure;
                %plot(w);
                %ylabel('Cumulative Area (%)');
                %xlabel('Province');
                opt = opt + 2;
            case 'yticklabel'
                ylbl = varargin{opt+1};
                opt = opt + 2;
            case 'cmap'
                cmap = varargin{opt+1};
                opt = opt + 2;
            case 'max'
                plotmax = true(1);
                opt = opt + 1;
            otherwise
                error('Unknown option.');
        end
    end
end

% create a pseudocolor plot for each pdf strip
[ny,nv] = size(p);
for i = 1:ny
    if isempty(w)
        h(i) = pcolor([v; v],repmat([-0.325; 0.325]+i,1,nv),[p(i,:); p(i,:)]);
    else
        h(i) = pcolor([v; v],repmat([w(i); w(i+1)],1,nv),[p(i,:); p(i,:)]);

    end
    shading flat;
    if i == 1
        hold on;
    end
end

if plotmax
    [pmax,ind] = max(p,[],2);
    if isempty(w)
        c = plot(v(ind),i,'k.');
    else
        c = plot(v(ind),midpt(w),'k.');
    end
    %set(c,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',4);
end
%[v(ind),midpt(w)]
if isempty(w)
    axis([v(1) v(end) 0.25 ny+0.75]);
    % add y labels
    if exist('ylbl','var')
        set(gca,'YTick',[1:ny],'YTickLabel',ylbl);
    end
else
    axis([v(1) v(end) min(w) max(w)]);
    ylabel('Area (%)');
end

% change from default colormap
if ~isnumeric(cmap)
    colormap2(cmap);
else
    colormap(cmap);
end

cb = colorbar;
cb.Label.String = 'Probability';

set(gca,'Layer','top','TickDir','out','XGrid','on','FontName','Myriad Pro');

if nargout == 1
    varargout{1} = h;
end

return