function [a, s]=getarrowslice
% Return the Arrow and Slice based on the GCO

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  if isempty(getappdata(gco,'controlarrow')) && ...
        isempty(getappdata(gco,'isosurface'))
    a = gco;
    s = getappdata(a,'arrowslice');
    if isempty(s)
      s=getappdata(a,'arrowiso');
    end
  else
    s = gco;
    if ~isempty(getappdata(s,'isosurface'))
      s=getappdata(s,'isosurface');
    end
    a = getappdata(s,'controlarrow');
  end

