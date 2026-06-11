function tf=propcheck(obj, prop, value)
% Check to see if PROP for OBJ has VALUE

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc
  
  try
    v=get(obj,prop);
  catch
    tf = 0;
    return
  end

  if isa(v,class(value))
    if isa(v,'char')
      tf=strcmp(v,value);
    else
      if v==value
        tf=1;
      else
        tf=0;
      end
    end
  else
    tf=0;
  end
