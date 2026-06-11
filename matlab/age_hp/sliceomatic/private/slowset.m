function slowset(handle, prop, value, increment)
% SLOWSET(H, 'PROPERTY', VALUE)
%
% Like SET, except that the property is set against the original
% value over several steps such that the value morphs from the
% starting value to the end value.
%
% H can be a vector of handles, but they must all accept PROPERTY.
%
% Only one property is supported at this time.
%
% Optional fourth argument INCREMENT specifies how many steps of
% animation to use.
  
% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  global INCREMENT;
  
  if nargin == 4
    INCREMENT = increment;
    if INCREMENT==0
      INCREMENT=1;
    end
  else
    INCREMENT=10;
  end
  
%  nargchk(3,3,'wrong number of arguments.');
  
  proplist=fieldnames(get(handle(1)));
  tprop={ proplist{strncmpi(prop,proplist,length(prop))} };
  prop=tprop{1};
  
  hp = [];
  
  for i = 1:length(handle)
    hp(i).handle = handle(i);
    hp(i).start = get(hp(i).handle,prop);
    hp(i).end = value;
    if isnumeric(hp(i).end) && isnumeric(hp(i).start)
      hp(i).values = VectorCalc(hp(i));
    else
      set(hp(i).handle,prop,value);
      hp(i).values = [];
    end
  end
  
  for inc = 1:INCREMENT
    for i = 1:length(handle)
      if ~isempty(hp(i).values)
        newval = reshape(hp(i).values(inc,:,:,:),...
                         size(hp(i).start,1),...
                         size(hp(i).start,2));
        
        set(hp(i).handle,prop,newval);
      end
    end
    pause(.05)
  end
  
function values =  VectorCalc(hp)
% Do nothing but go to end value.

  global INCREMENT;
  
  s = prod(size(hp.end));
  
  values = ones(INCREMENT,size(hp.end,1), size(hp.end,2),size(hp.end,3));
  
  for c = 1:s
    newval = linspace(hp.start(c),hp.end(c),INCREMENT);
    values(:,c) = newval';
  end
  
  values = reshape(values, INCREMENT, size(hp.end,1), size(hp.end,2), ...
                   size(hp.end,3));