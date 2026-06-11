function isocontrols(fig, onoff)
% Set up FIG to have an ISO surface controller on the bottom.
% ONOFF indicates if the controller is being turned ON or OFF

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

% Check variables
  error(nargchk(2,2,nargin))

  d = getappdata(fig, 'sliceomatic');
  
  if onoff

    lim=[min(min(min(d.data))) max(max(max(d.data)))];
      
    set(d.axiso,'handlevisibility','on');
    set(fig,'currentaxes',d.axiso);
    set(d.axiso, 'xlim',lim,...
                 'ylim',[1 5],...
                 'clim',lim);
    image('parent',d.axiso,'cdata',1:64,'cdatamapping','direct',...
          'xdata',lim,'ydata',[0 5],...
          'alphadata',.6, ...
          'hittest','off');
    activelabel('title','Iso Surface Controller');
    set(d.axiso,'handlevisibility','off');
    
  else
    % Turn off the controller
    
    delete(findobj(d.axis,'type','image'));

  end