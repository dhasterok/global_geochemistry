function showarrowtip (arrow)
% Display a tip for ARROW.
% Depends on tipdata being set on the handle to ARROW.
  
% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  d=getappdata(gcf,'sliceomatic');
  
  if isempty(arrow)
    tipdata = [];
  else
    
    ctrlarrow = getappdata(arrow,'controlarrow');
    
    if ~isempty(ctrlarrow)
      % In this case, the slice or isosurface passed in
      % most likely from the motion callback.  Lets redirect
      % our input argument and show the tip anyway.
      
      % See the "arrow" private function for the setting of this value;
      arrow = ctrlarrow(2);
    end
    
    tipdata = getappdata(arrow,'tipdata');
  end
  
  if ~isempty(tipdata)
  
    set(d.tip,'parent',tipdata.parentaxes, ...
              'string',sprintf('Value: %1.3f',tipdata.value),...
              'units','data', ...
              'position', tipdata.position, ...
              'verticalalignment', tipdata.verticalalign,...
              'horizontalalignment', tipdata.horizontalalign);
    set(d.tip,'units','pixels');
    set(d.tip,'visible','on');
    
  else
    
    set(d.tip,'visible','off');
    
  end
