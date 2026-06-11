function sliceomaticmotion(fig,action)
% Handle generic motion events for the figure window.

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  obj = hittest(fig);

  % Some objects get a special pointer when the mouse waves
  % over them.  Get it from appdata.
  if ~isempty(obj)
    t = getappdata(obj,'motionpointer');
    cc = get(fig,'pointer');
  
    if t
      newc = t;
    else
      newc = get(0,'defaultfigurepointer');
    end
  
    if isa(newc,'char') && isa(cc,'char') && ~strcmp(newc,cc)
      setpointer(fig, newc);
    end
  end
  
  d = getappdata(fig,'sliceomatic');

  % The motion meta slice is a line that is managed here.
  % There is only one.
  if isempty(d.motionmetaslice)
    d.motionmetaslice = animatedline('parent',d.axmain,... 
                    'vis','off',... 
                    'linestyle','--',... 
                    'marker','none',... 
                    'linewidth',2,... 
                    'clipping','off'); 
    setappdata(fig,'sliceomatic',d);
  end

  showarrowtip(obj);
  
  if isempty(obj) || (obj ~= d.axx && obj ~= d.axy && obj ~= d.axz)
    set(d.motionmetaslice,'visible','off');
    
    return
  end
  
  % OBJ can only be an Axes because of the previous IF statement.  
  aa = obj;
  apos=get(aa,'currentpoint');

  xl = d.xlim;
  yl = d.ylim;
  zl = d.zlim;
  
  if aa==d.axx || aa==d.axiso
    if aa==d.axiso
      % eh?
    else
      xdata = [ apos(1,1) apos(1,1) apos(1,1) apos(1,1) apos(1,1) ];
      ydata = [ yl(1) yl(2) yl(2) yl(1) yl(1) ];
      zdata = [ zl(2) zl(2) zl(1) zl(1) zl(2) ];
    end
  else
    % We are moving a Y or Z slice
    if aa==d.axy
      ydata = [ apos(1,2) apos(1,2) apos(1,2) apos(1,2) apos(1,2) ];
      xdata = [ xl(1) xl(2) xl(2) xl(1) xl(1) ];
      zdata = [ zl(2) zl(2) zl(1) zl(1) zl(2) ];
    else
      zdata = [ apos(1,2) apos(1,2) apos(1,2) apos(1,2) apos(1,2) ];
      ydata = [ yl(1) yl(2) yl(2) yl(1) yl(1) ];
      xdata = [ xl(2) xl(2) xl(1) xl(1) xl(2) ];
    end
  end

  addpoints(d.motionmetaslice,xdata,ydata,zdata); 
  drawnow