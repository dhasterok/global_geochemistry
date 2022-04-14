function outd = figmenus(d)
% Set up sliceomatic's gui menus within structure D

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

%%% 2/18/05 RAB part in 3 parts %%% 

% Main Figure Menu
  set(gcf,'menubar','none');
  
  % File menu
  d.filemenu = uimenu(gcf,'label','File');
  d.fcopy = uimenu(d.filemenu, 'label', 'Copy figure','callback', 'sliceomatic copy');
  d.fprint  = uimenu(d.filemenu,'label','Print...','callback','sliceomatic print');

  %%% start patch 1of3 RAB 2/18/05 %%%
  d.fsaveprefs = uimenu(d.filemenu,'label','Save preferences','callback',@SavePrefs);
  %%% end patch 1of3 RAB 2/18/05 %%%

  % How do get these props onto the print figure?
  %d.fprints = uimenu(d.filemenu,'label','Print Setup...','callback','printdlg -setup');
  % ---
  d.fexit = uimenu(d.filemenu, 'label', 'Close','callback','closereq',...
                   'separator','on');

  % Controls Menu
  d.defcontrols = uimenu(gcf,'label','Controls', 'callback',@controlmenu);
  if exist('uitoolfactory') == 2
    d.anntoolbar = uimenu(d.defcontrols,'label','Annotations toolbar','callback', 'sliceomatic annotationtoolbar');
  end
  d.camtoolbar = uimenu(d.defcontrols,'label','Camera toolbar','callback', 'sliceomatic cameratoolbar');
  d.dcalpha = uimenu(d.defcontrols,'label','Controls Transparency');
  d.dcalpha1= uimenu(d.dcalpha,'label','1','callback','sliceomatic controlalpha 1');
  d.dcalpha8= uimenu(d.dcalpha,'label','.8','callback','sliceomatic controlalpha .8');
  d.dcalpha6= uimenu(d.dcalpha,'label','.6','callback','sliceomatic controlalpha .6');
  d.dcalpha5= uimenu(d.dcalpha,'label','.5','callback','sliceomatic controlalpha .5');
  d.dcalpha4= uimenu(d.dcalpha,'label','.4','callback','sliceomatic controlalpha .4');
  d.dcalpha2= uimenu(d.dcalpha,'label','.2','callback','sliceomatic controlalpha .2');
  d.dcalpha0= uimenu(d.dcalpha,'label','0','callback','sliceomatic controlalpha 0');
  d.dcanimstep = uimenu(d.defcontrols,'label','Animation','callback', 'sliceomatic toggleanimation');
  d.dclabels= uimenu(d.defcontrols','label','Tick Labels','callback','sliceomatic controllabels');
  d.dcvis   = uimenu(d.defcontrols','label','Visible','callback','sliceomatic controlvisible');
  % d.dsetrange= uimenu(d.defcontrols','label','Set Range','callback','@setvolumerange');
  % d.dcslice = uimenu(d.defcontrols,'label','Slice Controls','callback','sliceomatic useslicecontrols');
  % d.dciso   = uimenu(d.defcontrols,'label','Iso Surface Control','callback','sliceomatic useisocontrols','separator','on');

  % Remove this once we have more controls to enable and disable.
  %  set(d.defcontrols,'vis','off');
  
  % Default for new slices menu
  d.defmenu = uimenu(gcf,'label','Object_Defaults', 'callback', @defaultmenu);
  d.dfacet  = uimenu(d.defmenu,'label','Slice Color Faceted','callback','sliceomatic defaultfaceted');
  d.dflat   = uimenu(d.defmenu,'label','Slice Color Flat',   'callback','sliceomatic defaultflat');
  d.dinterp = uimenu(d.defmenu,'label','Slice Color Interp', 'callback','sliceomatic defaultinterp');
  d.dtex    = uimenu(d.defmenu,'label','Slice Color Texture','callback','sliceomatic defaulttexture');
  d.dcnone  = uimenu(d.defmenu,'label','Slice Color None','callback','sliceomatic defaultcolornone');
  d.dtnone  = uimenu(d.defmenu,'label','Slice Transparency None','callback','sliceomatic defaulttransnone','separator','on');
  d.dtflat  = uimenu(d.defmenu,'label','Slice Transparency Flat','callback','sliceomatic defaulttransflat');
  d.dtinterp= uimenu(d.defmenu,'label','Slice Transparency Interp','callback','sliceomatic defaulttransinterp');
  d.dttex   = uimenu(d.defmenu,'label','Slice Transparency Texture','callback','sliceomatic defaulttranstexture');
  d.dlflat  = uimenu(d.defmenu,'label','IsoSurface Lighting Flat','callback','sliceomatic defaultlightflat','separator','on');
  d.dlsmooth= uimenu(d.defmenu,'label','IsoSurface Lighting Smooth','callback','sliceomatic defaultlightsmooth');
  %d.dcsmooth= uimenu(d.defmenu,'label','Contour Line Smoothing','callback','sliceomatic defaultcontoursmooth');
  d.dcflat  = uimenu(d.defmenu,'label','Contour Color Flat',   'callback','sliceomatic defaultcontourflat','separator','on');
  d.dcinterp= uimenu(d.defmenu,'label','Contour Color Interp', 'callback','sliceomatic defaultcontourinterp');
  d.dcblack = uimenu(d.defmenu,'label','Contour Color Black',  'callback','sliceomatic defaultcontourblack');
  d.dcwhite = uimenu(d.defmenu,'label','Contour Color White',  'callback','sliceomatic defaultcontourwhite');
  d.dclinew = uimenu(d.defmenu,'label','Contour Line Width');
  d.dcl1    = uimenu(d.dclinew,'label','1','callback','sliceomatic defaultcontourlinewidth 1');
  d.dcl2    = uimenu(d.dclinew,'label','2','callback','sliceomatic defaultcontourlinewidth 2');
  d.dcl3    = uimenu(d.dclinew,'label','3','callback','sliceomatic defaultcontourlinewidth 3');
  d.dcl4    = uimenu(d.dclinew,'label','4','callback','sliceomatic defaultcontourlinewidth 4');
  d.dcl5    = uimenu(d.dclinew,'label','5','callback','sliceomatic defaultcontourlinewidth 5');
  d.dcl6    = uimenu(d.dclinew,'label','6','callback','sliceomatic defaultcontourlinewidth 6');
  
  d.defcolor='texture';
  d.defalpha='texture';
  d.deflight='smooth';
  d.defcontourcolor='black';
  d.defcontourlinewidth=1;
  % This exposes an unpleasant R14 bug
  d.defcontoursmooth='off';

  % investigate hardware opengl.
  inc = 0;
  try
    od = opengl('data');
    if isfield(od,'Software')
      % R14 version of MATLAB
      if ~od.Software
        inc = 10;
      end
    else
      % Older version of MATLAB
      if ~(strcmp(od.Renderer,'Mesa X11') || ...
           strcmp(od.Renderer, 'GDI Generic'))
        inc = 10;
      end
    end
  end
  
  d.animincrement=inc;
  
  %%% start patch 2of3 RAB 2/18/05 %%%
  d = OverrideStickyUserPreferences(d);
  %%% end patch 2of3 RAB 2/18/05 %%%
  
  % Set props for all slices menu
  d.allmenu = uimenu(gcf,'label','AllSlices');
  uimenu(d.allmenu,'label','Color Faceted','callback','sliceomatic allfacet');
  uimenu(d.allmenu,'label','Color Flat','callback','sliceomatic allflat');
  uimenu(d.allmenu,'label','Color Interp','callback','sliceomatic allinterp');
  uimenu(d.allmenu,'label','Color Texture','callback','sliceomatic alltex');
  uimenu(d.allmenu,'label','Color None','callback','sliceomatic allnone');
  uimenu(d.allmenu,'label','Transparency None','callback','sliceomatic alltnone','separator','on');
  uimenu(d.allmenu,'label','Transparency .5','callback','sliceomatic alltp5');
  uimenu(d.allmenu,'label','Transparency Flat','callback','sliceomatic alltflat');
  uimenu(d.allmenu,'label','Transparency Interp','callback','sliceomatic alltinterp');
  uimenu(d.allmenu,'label','Transparency Texture','callback','sliceomatic allttex');

  % Setup Help style options
  d.helpmenu = uimenu(gcf,'label','Help');
  uimenu(d.helpmenu,'label','Help','callback','doc sliceomatic/sliceomatic');
  uimenu(d.helpmenu,'label','Check for Updates','callback','web http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=764&objectType=FILE');
  uimenu(d.helpmenu,'label','About Author','callback','web http://www.mathworks.com/matlabcentral/fileexchange/loadAuthor.do?objectId=803709&objectType=author');
  
  % Context Menus
  % Slice Context Menu
  d.uic=uicontextmenu('callback', @slicecontextmenu);
  d.vistog = uimenu(d.uic,'label','Visible','callback','sliceomatic togglevisible');
  d.uicdelete = uimenu(d.uic,'label','Delete','callback','sliceomatic deleteslice');
  d.smcolorm  = uimenu(d.uic,'label','Color','separator','on');
  d.smfacet   = uimenu(d.smcolorm,'label','Color Faceted','callback','sliceomatic setfaceted');
  d.smflat    = uimenu(d.smcolorm,'label','Color Flat','callback','sliceomatic setflat');
  d.sminterp  = uimenu(d.smcolorm,'label','Color Interp','callback','sliceomatic setinterp');
  d.smtex     = uimenu(d.smcolorm,'label','Color Texture','callback','sliceomatic settexture');
  d.smnone    = uimenu(d.smcolorm,'label','Color None','callback','sliceomatic setnone');
  d.smtransm  = uimenu(d.uic,'label','Transparency');
  d.smtnone   = uimenu(d.smtransm,'label','Transparency None','callback','sliceomatic setalphanone');
  d.smtp5     = uimenu(d.smtransm,'label','Transparency .5','callback','sliceomatic setalphapoint5');
  d.smtflat   = uimenu(d.smtransm,'label','Transparency Flat','callback','sliceomatic setalphaflat');
  d.smtinterp = uimenu(d.smtransm,'label','Transparency Interp','callback','sliceomatic setalphainterp');
  d.smttex    = uimenu(d.smtransm,'label','Transparency Texture','callback','sliceomatic setalphatexture');
  d.smcontour = uimenu(d.uic,'label','Add Contour','separator','on');
  d.smcont0   = uimenu(d.smcontour,'label','Auto (Slice)','callback','sliceomatic slicecontour');
  d.smcont0v  = uimenu(d.smcontour,'label','Auto (Volume)','callback','sliceomatic slicecontourfullauto');
  d.smcont1   = uimenu(d.smcontour,'label','Select Levels','callback','sliceomatic slicecontour_select','separator','on');
  d.smcsetauto= uimenu(d.uic,'label','Set Auto Levels (Slice)','callback','sliceomatic slicecontour_setauto');
  d.smcsetav  = uimenu(d.uic,'label','Set Auto Levels (Volume)','callback','sliceomatic slicecontour_setfullauto');
  d.smclevels = uimenu(d.uic,'label','Set Levels','callback','sliceomatic slicecontour_setlevels');
  d.smrcontour= uimenu(d.uic,'label','Remove Contour','callback','sliceomatic deleteslicecontour');
  d.smccm     = uimenu(d.uic,'label','Contour Colors');
  d.smcflat   = uimenu(d.smccm,'label','Contour Flat','callback','sliceomatic slicecontourflat');
  d.smcinterp = uimenu(d.smccm,'label','Contour Interp','callback','sliceomatic slicecontourinterp');
  d.smcblack  = uimenu(d.smccm,'label','Contour Black','callback','sliceomatic slicecontourblack');
  d.smcwhite  = uimenu(d.smccm,'label','Contour White','callback','sliceomatic slicecontourwhite');
  d.smccolor  = uimenu(d.smccm,'label','Contour Color','callback','sliceomatic slicecontourcolor');
  d.smcsmooth = uimenu(d.uic,'visible','off','label','Smooth Contour Lines','callback','sliceomatic slicecontoursmooth');
  d.smclinew  = uimenu(d.uic,'label','Contour Line Width');
  d.smcl1     = uimenu(d.smclinew,'label','1','callback','sliceomatic slicecontourlinewidth 1');
  d.smcl2     = uimenu(d.smclinew,'label','2','callback','sliceomatic slicecontourlinewidth 2');
  d.smcl3     = uimenu(d.smclinew,'label','3','callback','sliceomatic slicecontourlinewidth 3');
  d.smcl4     = uimenu(d.smclinew,'label','4','callback','sliceomatic slicecontourlinewidth 4');
  d.smcl5     = uimenu(d.smclinew,'label','5','callback','sliceomatic slicecontourlinewidth 5');
  d.smcl6     = uimenu(d.smclinew,'label','6','callback','sliceomatic slicecontourlinewidth 6');
  
  % Isosurface Context Menu
  d.uiciso=uicontextmenu('callback',@isocontextmenu);
  d.vistogiso = uimenu(d.uiciso,'label','Visible','callback','sliceomatic isotogglevisible');
  d.isodelete = uimenu(d.uiciso,'label','Delete','callback','sliceomatic isodelete');
  d.isoflatlight=uimenu(d.uiciso,'label','Lighting Flat','callback','sliceomatic isoflatlight','separator','on');
  d.isosmoothlight=uimenu(d.uiciso,'label','Lighting Smooth','callback','sliceomatic isosmoothlight');
  d.isocolor = uimenu(d.uiciso,'label','Change Color','callback','sliceomatic isocolor','separator','on');
  d.isoalpha=uimenu(d.uiciso,'label','Change Transparency');
  uimenu(d.isoalpha,'label','.2','callback','sliceomatic isoalpha .2');
  uimenu(d.isoalpha,'label','.5','callback','sliceomatic isoalpha .5');
  uimenu(d.isoalpha,'label','.8','callback','sliceomatic isoalpha .8');
  uimenu(d.isoalpha,'label','1','callback','sliceomatic isoalpha 1');
  d.isocap=uimenu(d.uiciso,'label','Add IsoCaps','callback','sliceomatic isocaps','separator','on');
  
  outd = d;
  
function controlmenu(fig, action)  
% Handle doing things to the CONTROLS menu
  
  d=getappdata(gcf,'sliceomatic');

  if cameratoolbar('getvisible')
    set(d.camtoolbar,'checked','on');
  else
    set(d.camtoolbar,'checked','off');
  end
  
  if exist('uitoolfactory') == 2
    if propcheck(d.toolbar,'visible','on')
      set(d.anntoolbar,'checked','on');
    else
      set(d.anntoolbar,'checked','off');
    end
  end
  
  set([d.dcalpha1 d.dcalpha8 d.dcalpha6 d.dcalpha5 d.dcalpha6 d.dcalpha2 d.dcalpha0...
       d.dclabels d.dcvis ],...
      'checked','off');

  switch get(d.pxx,'facealpha')
   case 1,  set(d.dcalpha1,'checked','on');
   case .8, set(d.dcalpha8,'checked','on');
   case .6, set(d.dcalpha6,'checked','on');
   case .5, set(d.dcalpha5,'checked','on');
   case .4, set(d.dcalpha4,'checked','on');
   case .2, set(d.dcalpha2,'checked','on');
   case 0,  set(d.dcalpha0,'checked','on');
  end

  if d.animincrement == 0
    set(d.dcanimstep,'checked','off');
  else
    set(d.dcanimstep,'checked','on');
  end
  
  if ~isempty(get(d.axx,'xticklabel'))
    set(d.dclabels,'checked','on');
  end
  
  if strcmp(get(d.axx,'visible'),'on')
    set(d.dcvis,'checked','on');
  end
  
  if 0
    xt = get(get(d.axx,'title'),'string');
    switch xt
     case 'X Slice Controller'
      set(d.dcslice,'checked','on');
    end
    
    xt = get(get(d.axiso,'title'),'string');
    switch xt
     case 'Iso Surface Controller'
      set(d.dciso,'checked','on');
    end
  end
  
function defaultmenu(fig, action)
% Handle toggling bits on the slice defaults menu
  
  d=getappdata(gcf,'sliceomatic');
  
  set([d.dfacet d.dflat d.dinterp d.dtex d.dtnone d.dtflat d.dtinterp ...
       d.dttex d.dcflat d.dcinterp d.dcblack d.dcwhite d.dcnone ...
       d.dlflat d.dlsmooth ...
       d.smcl1 d.smcl2 d.smcl3 d.smcl4 d.smcl5 d.smcl6 ], 'checked','off');
  switch d.defcolor
   case 'faceted'
    set(d.dfacet,'checked','on');
   case 'flat'
    set(d.dflat,'checked','on');
   case 'interp'
    set(d.dinterp,'checked','on');
   case 'texture'
    set(d.dtex,'checked','on');
   case 'none'
    set(d.dcnone,'checked','on');
  end
  switch d.defalpha
   case 'none'
    set(d.dtnone,'checked','on');
   case 'flat'
    set(d.dtflat,'checked','on');
   case 'interp'
    set(d.dtinterp,'checked','on');
   case 'texture'
    set(d.dttex,'checked','on');
  end
  switch d.deflight
   case 'flat'
    set(d.dlflat,'checked','on');
   case 'smooth'
    set(d.dlsmooth,'checked','on');
  end
  switch d.defcontourcolor
   case 'flat'
    set(d.dcflat,'checked','on');
   case 'interp'
    set(d.dcinterp,'checked','on');
   case 'black'
    set(d.dcblack,'checked','on');
   case 'white'
    set(d.dcwhite,'checked','on');
  end
  %set(d.dcsmooth,'checked',d.defcontoursmooth);
  switch d.defcontourlinewidth
   case 1, set(d.dcl1,'checked','on');
   case 2, set(d.dcl2,'checked','on');
   case 3, set(d.dcl3,'checked','on');
   case 4, set(d.dcl4,'checked','on');
   case 5, set(d.dcl5,'checked','on');
   case 6, set(d.dcl6,'checked','on');
  end

function slicecontextmenu(fig,action)
% Context menu state for slices

  d=getappdata(gcf,'sliceomatic');

  [a s]=getarrowslice;
  set([d.smfacet d.smflat d.sminterp d.smtex d.smtnone d.smtp5 ...
       d.smtflat d.smtinterp d.smttex d.smnone d.smcsmooth
      ],'checked','off');
  set(d.vistog,'checked',get(s,'visible'));
  
  if propcheck(s,'edgec',[0 0 0])
    set(d.smfacet,'checked','on');
  elseif propcheck(s,'facec','flat')
    set(d.smflat,'checked','on');
  end
  if propcheck(s,'facec','interp')
    set(d.sminterp,'checked','on');
  end
  if propcheck(s,'facec','texturemap')
    set(d.smtex,'checked','on');
  end
  if propcheck(s,'facec','none')
    set(d.smnone,'checked','on');
  end
  if propcheck(s,'facea',1)
    set(d.smtnone,'checked','on');
  end
  if propcheck(s,'facea',.5)
    set(d.smtp5,'checked','on');
  end
  if propcheck(s,'facea','flat')
    set(d.smtflat,'checked','on');
  end
  if propcheck(s,'facea','interp')
    set(d.smtinterp,'checked','on');
  end
  if propcheck(s,'facea','texturemap')
    set(d.smttex,'checked','on');
  end
  cm = [d.smcflat d.smcinterp d.smcblack d.smcwhite d.smccolor ...
       d.smcl1 d.smcl2 d.smcl3 d.smcl4 d.smcl5 d.smcl6 ];
  set(cm,'checked','off');
  if isempty(getappdata(s,'contour'))
    set(d.smcontour,'enable','on');
    set(d.smcsetauto,'enable','off');
    set(d.smcsetav,'enable','off');
    set(d.smclevels,'enable','off');
    set(d.smrcontour,'enable','off');
    set(d.smcsmooth,'enable','off');
    set(cm,'enable','off');
  else
    set(d.smcontour,'enable','off')
    set(d.smcsetauto,'enable','on');
    set(d.smcsetav,'enable','on');
    set(d.smclevels,'enable','on');
    set(d.smrcontour,'enable','on')
    set(d.smcsmooth,'enable','on');
    set(cm,'enable','on')
    c = getappdata(s,'contour');
    if propcheck(c,'linesmoothing','on')
      set(d.smcsmooth,'checked','on');
    end
    ec = get(c,'edgecolor');
    if isa(ec,'char')
      switch ec
       case 'flat'
        set(d.smcflat,'checked','on');
       case 'interp'
        set(d.smcinterp,'checked','on');
      end
    else
      if ec == [ 1 1 1 ]
        set(d.smcwhite,'checked','on');
      elseif ec == [ 0 0 0 ]
        set(d.smcblack,'checked','on');
      else
        set(d.smccolor,'checked','on');
      end
    end
    clw = get(c,'linewidth');
    switch clw
     case 1, set(d.smcl1,'checked','on');
     case 2, set(d.smcl2,'checked','on');
     case 3, set(d.smcl3,'checked','on');
     case 4, set(d.smcl4,'checked','on');
     case 5, set(d.smcl5,'checked','on');
     case 6, set(d.smcl6,'checked','on');
    end
  end
  
function isocontextmenu(fig,action)
% Context menu state for isosurfaces

  d=getappdata(gcf,'sliceomatic');
  
  [a s]=getarrowslice;
  if propcheck(s,'facelighting','flat')
    set(d.isoflatlight,'checked','on');
    set(d.isosmoothlight,'checked','off');
  else
    set(d.isoflatlight,'checked','off');
    set(d.isosmoothlight,'checked','on');
  end
  set(d.vistogiso,'checked',get(s,'visible'));
  if ~isempty(getappdata(s,'isosurfacecap'))
    set(d.isocap,'checked','on');
  else
    set(d.isocap,'checked','off');
  end


%%% start patch 3of3 RAB 2/18/05 %%%
%----------------------------------------------------------------------
function SavePrefs(obj,event)

%appdata structure knows everything about the implementation
d = getappdata(gcf,'sliceomatic');

%extract only preferences that need to be sticky
prefs.anntoolbar_Checked = get(d.toolbar,'Visible');
prefs.defcolor = d.defcolor;
prefs.defalpha = d.defalpha;
prefs.deflight = d.deflight;
prefs.defcontourcolor = d.defcontourcolor;
prefs.defcontourlinewidth = d.defcontourlinewidth;
prefs.defcontoursmooth = d.defcontoursmooth;

prefs.camtoolbar_checked = cameratoolbar('getvisible');
prefs.ticklabels = get(d.axx,'xticklabelmode');
prefs.animincrement = d.animincrement;
prefs.controlalpha = get(d.pxx,'facealpha');

%store mini structure (locally where Slice-O-Matic installed)
fileName = UserStickyPrefsFileName;
save(fileName,'prefs')
disp([ 'Saved: ' fileName])


%----------------------------------------------------------------------
function dOut = OverrideStickyUserPreferences(d)

%characteristic prefs file (stored locally where Slice-O-Matic installed)
fileName = UserStickyPrefsFileName;

%override particular field values (if file exists)
if exist(fileName,'file')
  load(fileName)
  set(d.toolbar,'visible',prefs.anntoolbar_Checked)
  if prefs.camtoolbar_checked
    cameratoolbar('show');
  else
    cameratoolbar('hide');
  end
  d.defcolor = prefs.defcolor;
  d.defalpha = prefs.defalpha;
  d.deflight = prefs.deflight;
  d.defcontourcolor = prefs.defcontourcolor;
  d.defcontourlinewidth = prefs.defcontourlinewidth;
  d.defcontoursmooth = prefs.defcontoursmooth;
  
  if strcmp('auto',prefs.ticklabels)
    set([d.axx d.axiso],'xticklabelmode','auto');
    set([d.axy d.axz],'yticklabelmode','auto');
  else
    set([d.axx d.axiso],'xticklabel',[]);
    set([d.axy d.axz],'yticklabel',[]);
  end
  d.animincrement = prefs.animincrement;
  set([d.pxx d.pxy d.pxz] , 'facealpha',prefs.controlalpha);
  iso = findobj(d.axiso,'type','image');
  set(iso,'alphadata',prefs.controlalpha);

  disp('Sticky preferences loaded.')
end

%return modified structure
dOut = d;


%----------------------------------------------------------------------
function fileName = UserStickyPrefsFileName
localPath = fileparts(which(mfilename));
fileName = fullfile(localPath,'Sliceomatic.Prefs.mat');

%%% end patch 3of3 RAB 2/18/05 %%%
