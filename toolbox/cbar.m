function cb = cbar(varargin);
% CBAR - Makes a more compact colorbar.
%
%     Matlab's colobar is big and bulky.  So CBAR reduces the size
%     of the colorbar.  It also allow the user to set some properties
%     easily without the use of handles.
%
%     CBAR(TTXT) sets the title above the colorbar.
%
%     CBAR(TTXT,TICK,TICKLBL) sets the tick locations and labels.
%
% Last Modified: 22-Feb 2008 by D. Hasterok

% Create a colorbar
cb = colorbar;
% Shrink the size
pos = get(cb,'Position');
set(cb, ...
    'Position',[pos(1) pos(2)+pos(4)/2.75 pos(3)/1.5 pos(4)/4]);

if nargin > 0
    % Add a title
    if isa(varargin{1},'char')
        %set(cb,'XAxisLocation','top');
        xlbl = get(cb,'YLabel');
        set(xlbl,'String',varargin{1});
    else
        error('Error: Colorbar title must be a string.');
    end
    % Change the tick locations
    if nargin > 1
        if isa(varargin{2},'double')
            set(cb,'YTick',varargin{2});
        else
            error('Error: Colorbar ticks must be a vector of doubles.');
        end
        % Change the tick labels
        if nargin > 2
            if isa(varargin{3},'double') | isa(varargin{3},'char');
                set(cb,'YTickLabel',varargin{3});
            else
                error('Error: Colorbar Labels must be a vector of strings or doubles.');
            end
        else
            error('Error: Too many arguments.');
        end
    end
end

return
