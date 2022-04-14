function activelabel(label, string)
% ACTIVELABEL(LABEL, STRING) - Create a label on GCA which is
%     active.  LABEL is the property of GCA whose label you are
%     setting.  STRING is the initial text string for the label.

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  l = get(gca,label);
  
  set(l,'string',string);
  set(l,'buttondownfcn',@activelabelbuttondown);
  
function activelabelbuttondown(obj, action)
% Callback when one of our active labels is clicked on.
  
  set(obj,'edit','on');
  
    