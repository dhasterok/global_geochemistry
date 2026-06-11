function outd = figtoolbar(d)
% Set up the toolbar for Sliceomatic within structure D

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  set(gcf,'toolbar','none');
  
  if exist('uitoolfactory') == 2
    
    % Create a toolbar with just the elements useful on sliceomatic
    % on it.

    d.toolbar = uitoolbar('parent',gcf);
    uitoolfactory(d.toolbar, 'Annotation.InsertRectangle');
    uitoolfactory(d.toolbar, 'Annotation.InsertEllipse');
    uitoolfactory(d.toolbar, 'Annotation.InsertTextbox');
    uitoolfactory(d.toolbar, 'Annotation.InsertArrow');
    uitoolfactory(d.toolbar, 'Annotation.InsertLine');
    uitoolfactory(d.toolbar, 'Exploration.ZoomIn');
    uitoolfactory(d.toolbar, 'Exploration.ZoomOut');
    uitoolfactory(d.toolbar, 'Exploration.Pan');
    uitoolfactory(d.toolbar, 'Exploration.Rotate');
    
    cameratoolbar('show');
    cameratoolbar('togglescenelight');
    
  else

    % We are in R13 or earlier
    try
      cameratoolbar('show');
      cameratoolbar('togglescenelight');
      %cameratoolbar('setmode','orbit');
    catch
      disp('Could not display the camera toolbar.');
    end
    
  end
  
  outd = d;