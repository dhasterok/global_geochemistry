function levels = selectcontourlevels(data, min, max)
% L = SELECTCONTOURLEVELS(DATA) - select contour levels for DATA.
%              Returns a vector of levels.
% L = SELECTCONTOURLEVELS(DATA, min, max) - select contour levels for DATA.
%              MINinimum and MAXimum possible contour levels
%              specified.
  
% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  dlg_title = 'Sliceomatic contour select';
  
  if nargin == 1
    min = min(min(data));
    max = max(max(data));
  end

  levels = [];
  
  prompt = ['Enter contour values over '...
            num2str(min) ' to '...
            num2str(max) ':'];
  answer = inputdlg(prompt,dlg_title);
  if isempty(answer)~=1 & isequal(answer{1},'')~=1
    levels = str2num(answer{1});
    if isempty(levels)==1
      warndlg('Contour values must be numeric: data not accepted.');
    end
  else
    warndlg('Empty value: data not accepted.');
  end
  
% End

