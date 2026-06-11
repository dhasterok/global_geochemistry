function pushset(handle,prop,value)
% PUSHSET - push new properties onto a value stack.
%
% PUSHSET(HANDLE, PROP, VALUE) will take the old value of PROP for
% HANDLE, and save it on a stack associated with HANDLE.  It then
% assigns VALUE as the new value.
%

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  proplist=fieldnames(get(handle(1)));
  prop=proplist{strcmpi(prop,proplist)};
  appstr = [prop '_hgstack'];

  for k=1:prod(size(handle))

    oldv = get(handle(k),prop);
    olds = getappdata(handle(k),appstr);
    set(handle(k),prop,value);
    setappdata(handle(k),appstr,{ oldv olds });
    
  end
