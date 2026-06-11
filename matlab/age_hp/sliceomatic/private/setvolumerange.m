function setvolumerange
% Query for a new volume range based on the sliceomatic gui
% which should be GCF
  
% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  %d=getappdata(fig,'sliceomatic');
  
  p=get(fig,'position');
  np=[p(1)+20 p(2)+30 400 200];
  figure('position',np);
  
  uicontrol('units','norm','style','text','string','X Range',...
            'position',[0 .6 .3 .3]);
  uicontrol('units','norm','style','text','string','Y Range',...
            'position',[0 .3 .3 .3]);
  uicontrol('units','norm','style','text','string','Z Range',...
            'position',[0 0  .3 .3]);
  