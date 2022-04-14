function sliceomatic1p0(p1,p2,p3,p4,p5)
% SLICEOMATIC - Slice and isosurface volume exploration GUI
%
% SLICEOMATIC(DATA) - Use 3D double matrix DATA as a volume data
% SLICEOMATIC(X,Y,Z,DATA) - Use 3d double matrix DATA as a volume
%       Set the volume dimensions as X, Y, and Z.
%       (Not currently used.)
%
% Example:
%
%       [x,y,z] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2);
%       v = x .* exp(-x.^2 - y.^2 - z.^2);
%       sliceomatic(v)
%
% Using the GUI:
% -------------
% The white bars on the top, left, and right allow insertion of
% new slices on the X, Y, and Z planes.  Click in an empty area to
% add a new slice or surface.  Click on a control arrow to add a
% new slice or surface.
%
% The colored bar at the bottom is used to place and position an
% isosurface.  The color in the bar indicates a position (as seen
% in the slice) where the isosurface will go.
%
% When the rotate camera button is on, the popup menu will control
% the camera.  Turn off camera rotation in order to get individual
% control over properties of the slices and isosurfaces.
%
% The defaults menu provides default features of newly created
% slices and surfaces.  The AllSlices menu controls properties of
% all the slices and surfaces in the scene.  Use popup menus on the
% objects themselves, or the control arrows to change individual
% properties.
%
% If the data is very large, a reduced model of the data is created.
% This reduced data set is used for interactively positioning
% isosurfaces.
%
% The Colormap popdown controls the currently active colormap.
% This map is used to color the slices.  The Alphamap popdown
% controls the alphamap used on the slices.
%
% Doing Cool Stuff:
% ----------------
%
% Exploration:
% You can get a quick feel of the current data set by adding a
% slice using the ColorTexture option.  Such a slice can be dragged
% through the data very quickly.
%
% Highlight an Area:
% If certain values in your data are interesting (very large, very
% small, or very median values) you can use transparency to make
% parts of your slices disappear.  Choose AlphaTexture options from
% the defaults, and sweep some slices across your data.  Use the
% AlphaMap to pick out certain data sets.  The example given here
% looks best with the `vdown' alphamap.
%
% Hidden Shapes:
% Use the isosurface control bar to create an isosurface.  Be
% patient with this control.  It often takes a while for the
% surface to be created.  Click and hold the mouse button down
% until the first surface appears.  Drag the surface through the
% values until you get something you like, then let go.  If your
% data set is very large, you will need to wait while the new and
% more accurate isosurface is created.
%
% Volumes:
% You can simulate a volume object by creating lots of stacked
% slices.  Simply use the proper Alphamap and transparent textures
% to highlight the correct data, and a large stack of slices will
% let you simulate a volume object.
%
% BUGS:
% ----
%
% 1) Sliceomatic does not use the `slice' command.  All slices are
%    created by explicitly extracting data from the volume.  As such,
%    only slices at integer values are allowed.
% 

%
%
% Sliceomatic is a tool I wrote for fun.  There are no warrenties
% expressed or implied.

% Written by Eric Ludlam <eludlam@mathworks.com>
% Copyright 2000, 2001 The MathWorks Inc

if isa(p1,'double')
  if nargin==4
    d.Xv=p1;
    d.Yv=p2;
    d.Zv=p3;
    p1=p4;
  else
    d.Yv=1:size(p1,1);
    d.Xv=1:size(p1,2);
    d.Zv=1:size(p1,3);
  end
  
  % Init sliceomatic
  figure('name','Sliceomatic','toolbar','none');
  d.data=p1;
  lim=[min(min(min(p1))) max(max(max(p1)))];
  d.axmain = axes('units','normal','pos',[.2  .2 .6 .6],'box','on',...
		  'ylim',[1 size(p1,1)],...
		  'xlim',[1 size(p1,2)],...
		  'zlim',[1 size(p1,3)],...
		  'clim',lim,...
		  'alim',lim);
  view(3);
  axis tight vis3d;
  hold on;
  grid on;
  xlabel X;
  ylabel Y;
  zlabel Z;
  daspect([1 1 1]);
  %title('SliceOMatic');
  cameratoolbar('togglescenelight');
  d.axx    = axes('units','normal','pos',[.2  .81 .6 .1],'box','on',...
		  'ytick',[],'xgrid','on','xaxislocation','top',...
		  'xlim',[1 size(p1,2)],...
		  'ylim',[1 5],...
		  'layer','top',...
		  'color','white');
  title('X Slice Controller');
  d.axy    = axes('units','normal','pos',[.05 .05 .1 .75],'box','on',...
		  'xtick',[],'ygrid','on',...
		  'xlim',[1 5],...
		  'ylim',[1 size(p1,1)],...
		  'layer','top',...
		  'color','white');
  title('Y Slice');
  d.axz    = axes('units','normal','pos',[.85 .05 .1 .75],'box','on',...
		  'xtick',[],'ygrid','on','yaxislocation','right',...
		  'xlim',[1 5],...
		  'ylim',[1 size(p1,3)],...
		  'layer','top',...
		  'color','white');
  title('Z Slice');
  d.axiso  = axes('units','normal','pos',[.2 .05 .6 .1],'box','on',...
		  'ytick',[],'xgrid','off','ygrid','off','xaxislocation','bottom',...
		  'xlim',lim,'ylim',[1 5],'clim',lim,'color','none',...
		  'layer','top');
  title('Iso Surface Controller');
  image('parent',d.axiso,'cdata',1:64,'cdatamapping','direct',...
	'xdata',lim,'ydata',[0 5],...
	'hittest','off');
  set([d.axx d.axy d.axz d.axiso],'handlevisibility','off');
  try
    cameratoolbar('show');
    cameratoolbar('setmode','orbit');
  end
  set(d.axx,'buttondownfcn','sliceomatic Xnew');
  set(d.axy,'buttondownfcn','sliceomatic Ynew');
  set(d.axz,'buttondownfcn','sliceomatic Znew');
  set(d.axiso,'buttondownfcn','sliceomatic ISO');
  
  d.uic=uicontextmenu('callback','sliceomatic contextmenu');
  d.vistog = uimenu(d.uic,'label','Visible','callback','sliceomatic togglevisible');
  uimenu(d.uic,'label','Delete','callback','sliceomatic deleteslice');
  d.smfacet   = uimenu(d.uic,'label','Color Faceted','callback','sliceomatic setfaceted');
  d.smflat    = uimenu(d.uic,'label','Color Flat','callback','sliceomatic setflat');
  d.sminterp  = uimenu(d.uic,'label','Color Interp','callback','sliceomatic setinterp');
  d.smtex     = uimenu(d.uic,'label','Color Texture','callback','sliceomatic settexture');
  d.smtnone   = uimenu(d.uic,'label','Transparency None','callback','sliceomatic setalphanone');
  d.smtp5     = uimenu(d.uic,'label','Transparency .5','callback','sliceomatic setalphapoint5');
  d.smtflat   = uimenu(d.uic,'label','Transparency Flat','callback','sliceomatic setalphaflat');
  d.smtinterp = uimenu(d.uic,'label','Transparency Interp','callback','sliceomatic setalphainterp');
  d.smttex    = uimenu(d.uic,'label','Transparency Texture','callback','sliceomatic setalphatexture');

  d.uiciso=uicontextmenu('callback','sliceomatic contextiso');
  d.vistogiso = uimenu(d.uiciso,'label','Visible','callback','sliceomatic isotogglevisible');
  uimenu(d.uiciso,'label','Delete','callback','sliceomatic isodelete');
  d.isoflatlight=uimenu(d.uiciso,'label','Lighting Flat','callback','sliceomatic isoflatlight');
  d.isosmoothlight=uimenu(d.uiciso,'label','Lighting Smooth','callback','sliceomatic isosmoothlight');
  uimenu(d.uiciso,'label','Change Color','callback','sliceomatic isocolor');
  d.isoalpha=uimenu(d.uiciso,'label','Change Transparency');
  uimenu(d.isoalpha,'label','.2','callback','sliceomatic isoalpha .2');
  uimenu(d.isoalpha,'label','.5','callback','sliceomatic isoalpha .5');
  uimenu(d.isoalpha,'label','.8','callback','sliceomatic isoalpha .8');
  uimenu(d.isoalpha,'label','1','callback','sliceomatic isoalpha 1');
  d.isocap=uimenu(d.uiciso,'label','Add IsoCaps','callback','sliceomatic isocaps');

  d.defmenu = uimenu(gcf,'label','Defaults', 'callback','sliceomatic defaultmenu');
  d.dfacet  = uimenu(d.defmenu,'label','Color Faceted','callback','sliceomatic defaultfaceted');
  d.dflat   = uimenu(d.defmenu,'label','Color Flat',   'callback','sliceomatic defaultflat');
  d.dinterp = uimenu(d.defmenu,'label','Color Interp', 'callback','sliceomatic defaultinterp');
  d.dtex    = uimenu(d.defmenu,'label','Color Texture','callback','sliceomatic defaulttexture');
  d.dtnone  = uimenu(d.defmenu,'label','Transparency None','callback','sliceomatic defaulttransnone');
  d.dtflat  = uimenu(d.defmenu,'label','Transparency Flat','callback','sliceomatic defaulttransflat');
  d.dtinterp= uimenu(d.defmenu,'label','Transparency Interp','callback','sliceomatic defaulttransinterp');
  d.dttex   = uimenu(d.defmenu,'label','Transparency Texture','callback','sliceomatic defaulttranstexture');
  d.dlflat  = uimenu(d.defmenu,'label','Lighting Flat','callback','sliceomatic defaultlightflat');
  d.dlsmooth  = uimenu(d.defmenu,'label','Lighting Smooth','callback','sliceomatic defaultlightsmooth');

  d.defcolor='flat';
  d.defalpha='none';
  d.deflight='smooth';

  d.allmenu = uimenu(gcf,'label','AllSlices');
  uimenu(d.allmenu,'label','Color Faceted','callback','sliceomatic allfacet');
  uimenu(d.allmenu,'label','Color Flat','callback','sliceomatic allflat');
  uimenu(d.allmenu,'label','Color Interp','callback','sliceomatic allinterp');
  uimenu(d.allmenu,'label','Color Texture','callback','sliceomatic alltex');
  uimenu(d.allmenu,'label','Transparency None','callback','sliceomatic alltnone');
  uimenu(d.allmenu,'label','Transparency .5','callback','sliceomatic alltp5');
  uimenu(d.allmenu,'label','Transparency Flat','callback','sliceomatic alltflat');
  uimenu(d.allmenu,'label','Transparency Interp','callback','sliceomatic alltinterp');
  uimenu(d.allmenu,'label','Transparency Texture','callback','sliceomatic allttex');
  
  uicontrol('style','text','string','ColorMap',...
	    'units','normal','pos',[0 .9 .19 .1]);
  uicontrol('style','popup','string',...
	    {'jet','hsv','cool','hot','pink','bone','copper','flag','prism','rand'},...
	    'callback','sliceomatic colormap',...
	    'units','normal','pos',[0 .85 .19 .1]);

  uicontrol('style','text','string','AlphaMap',...
	    'units','normal','pos',[.81 .9 .19 .1]);
  uicontrol('style','popup','string',{'rampup','rampdown','vup','vdown','rand'},...
	    'callback','sliceomatic alphamap',...
	    'units','normal','pos',[.81 .85 .19 .1]);

  disp('Smoothing for IsoNormals...');
  d.smooth=smooth3(d.data);  % ,'box',5);
  d.reducenumbers=[floor(size(d.data,2)/20)...
		   floor(size(d.data,1)/20)...
		   floor(size(d.data,3)/20) ];
  d.reducenumbers(d.reducenumbers==0)=1;
  % Vol vis suite takes numbers in X/Y form.
  ly = 1:d.reducenumbers(1):size(d.data,2);
  lx = 1:d.reducenumbers(2):size(d.data,1);
  lz = 1:d.reducenumbers(3):size(d.data,3);

  d.reducelims={ ly lx lz };
  disp('Generating reduction volume...');
  d.reduce= reducevolume(d.data,d.reducenumbers);
  d.reducesmooth=smooth3(d.reduce,'box',5);

  setappdata(gcf,'sliceomatic',d);
else
  % Interpret commands
  d=getappdata(gcf,'sliceomatic');
  try
    switch p1
     case 'Xnew'
      if strcmp(get(gcf,'selectiontype'),'normal')
	pt=get(gcbo,'currentpoint');
	axis(gcbo);
	X=pt(1,1);
	newa=arrow(gcbo,'down',[X 0]);
	set(gcf,'currentaxes',d.axmain);
	new=localslice(d.data, X, [], []);
	setappdata(new,'controlarrow',newa);
	setappdata(newa(2),'arrowslice',new);
	set(new,'alphadata',get(new,'cdata'),'alphadatamapping','scaled');
	set(newa,'buttondownfcn','sliceomatic Xmove');
	set([new newa],'uicontextmenu',d.uic);
	% Make sure whatever buttonupfcn on the figure is run now to "turn
	% off" whatever was going on before we got our callback on the
	% arrow.
	buf = get(gcf,'windowbuttonupfcn');
	if ~strcmp(buf,'')
	  eval(buf);
	end
	d.draggedarrow=newa(2);
	dragprep(newa(2));
      end
     case 'Ynew'
      if strcmp(get(gcf,'selectiontype'),'normal')
	pt=get(gcbo,'currentpoint');
	Y=pt(1,2);
	newa=arrow(gcbo,'right',[0 Y]);
	set(gcf,'currentaxes',d.axmain);
	new=localslice(d.data, [], Y, []);
	setappdata(new,'controlarrow',newa);
	setappdata(newa(2),'arrowslice',new);
	set(new,'alphadata',get(new,'cdata'),'alphadatamapping','scaled');
	set(newa,'buttondownfcn','sliceomatic Ymove');
	set([new newa],'uicontextmenu',d.uic);
	% Make sure whatever buttonupfcn on the figure is run now to "turn
	% off" whatever was going on before we got our callback on the
	% arrow.
	buf = get(gcf,'windowbuttonupfcn');
	if ~strcmp(buf,'')
	  eval(buf);
	end
	d.draggedarrow=newa(2);
	dragprep(newa(2));
      end % if strcmp(get(gcf,
     case 'Znew'
      if strcmp(get(gcf,'selectiontype'),'normal')
	pt=get(gcbo,'currentpoint');
	Y=pt(1,2);
	newa=arrow(gcbo,'left', [0 Y]);
	set(gcf,'currentaxes',d.axmain);
	new=localslice(d.data, [], [], Y);
	set(new,'alphadata',get(new,'cdata'),'alphadatamapping','scaled');
	setappdata(new,'controlarrow',newa);
	setappdata(newa(2),'arrowslice',new);
	set(newa,'buttondownfcn','sliceomatic Zmove');
	set([new newa],'uicontextmenu',d.uic);
	% Make sure whatever buttonupfcn on the figure is run now to "turn
	% off" whatever was going on before we got our callback on the
	% arrow.
	buf = get(gcf,'windowbuttonupfcn');
	if ~strcmp(buf,'')
	  eval(buf);
	end
	d.draggedarrow=newa(2);
	dragprep(newa(2));
      end % if strcmp(get(gcf,
     case 'ISO'
      if strcmp(get(gcf,'selectiontype'),'normal')
	pt=get(gcbo,'currentpoint');
	V=pt(1,1);
	newa=arrow(gcbo,'up',[V 0]);
	set(gcf,'currentaxes',d.axmain);
	new=localisosurface(d.reducelims,d.reduce,d.reducesmooth,V);
	set([newa new],'uicontextmenu',d.uiciso);
	setappdata(new,'controlarrow',newa);
	setappdata(new,'reduced',1);
	setappdata(newa(2),'arrowiso',new);
	set(newa,'buttondownfcn','sliceomatic ISOmove');
	% Make sure whatever buttonupfcn on the figure is run now to "turn
	% off" whatever was going on before we got our callback on the
	% arrow.
	buf = get(gcf,'windowbuttonupfcn');
	if ~strcmp(buf,'')
	  eval(buf);
	end	
	d.draggedarrow=newa(2);
	dragprep(newa(2));
      end % if strcmp(get(gcf,
     case 'Xmove'
      if strcmp(get(gcf,'selectiontype'),'normal')
	[a s]=getArrowSlice;
	d.draggedarrow=a;
	dragprep(a);
      end
     case 'Ymove'
      if strcmp(get(gcf,'selectiontype'),'normal')
	[a s]=getArrowSlice;
	d.draggedarrow=a;
	dragprep(a);
      end
     case 'Zmove'
      if strcmp(get(gcf,'selectiontype'),'normal')
	[a s]=getArrowSlice;
	d.draggedarrow=a;
	dragprep(a);
      end
     case 'ISOmove'
      if strcmp(get(gcf,'selectiontype'),'normal')
	[a s]=getArrowSlice;
	d.draggedarrow=a;
	dragprep(a);
      end      
     case 'up'
      if strcmp(get(gcf,'selectiontype'),'normal')
	dragfinis(d.draggedarrow);
      end
     case 'motion'
      a=d.draggedarrow;			% The arrow being dragged
      s=getappdata(a,'arrowslice');	% The slice to 'move'
      if isempty(s)
	s=getappdata(a,'arrowiso');	% or the isosurface
      end
      aa=get(a,'parent');		% arrow's parent axes
      pos=getappdata(a,'arrowcenter');	% the line the arrow points at.
      apos=get(aa,'currentpoint');
      if aa==d.axx | aa==d.axiso
	% We are moving an X slice
	xdiff=apos(1,1)-pos(1,1);
	v=get(a,'vertices');
	v(:,1)=v(:,1)+xdiff;
	set([a getappdata(a,'arrowedge')],'vertices',v);
	np=[ apos(1,1) 0 ];
	% This might be a slice, or an isosurface!
	if aa==d.axiso
	  new=localisosurface(d.reducelims,d.reduce,d.reducesmooth,...
			      apos(1,1),s);
	  setappdata(new,'reduced',1);
	else
	  if round(apos(1,1))~=round(pos(1,1))
	    new=localslice(d.data, apos(1,1), [], [],s);
	  end
	end
      else
	% We are moving a Y or Z slice
	ydiff=apos(1,2)-pos(1,2);
	v=get(a,'vertices');
	v(:,2)=v(:,2)+ydiff;
	set([a getappdata(a,'arrowedge')],'vertices',v);
	np=[ 0 apos(1,2) ];
	if round(apos(1,2))~=round(pos(1,2))
	  if aa==d.axy
	    new=localslice(d.data, [], apos(1,2), [], s);
	  else
	    new=localslice(d.data, [], [], apos(1,2), s);
	  end
	end
      end
      setappdata(a,'arrowcenter',np);
      drawnow;
      %
      % IsoSurface context menu items
      %
     case 'contextiso'
      [a s]=getArrowSlice;
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
     case 'isotogglevisible'
      [a s]=getArrowSlice;
      if propcheck(s,'visible','on')
	set(s,'visible','off');
      else
	set(s,'visible','on');
      end
     case 'isodelete'
      [a s]=getArrowSlice;
      delete(s);
      if prod(size(a))==1
	delete(getappdata(a,'arrowedge'));
      end
      cap=getappdata(s,'sliceomaticisocap');
      if ~isempty(cap)
	delete(cap);
      end
      delete(a);
     case 'isoflatlight'
      [a s]=getArrowSlice;
      set(s,'facelighting','flat');
     case 'isosmoothlight'
      [a s]=getArrowSlice;
      set(s,'facelighting','phong');
     case 'isocolor'
      [a s]=getArrowSlice;
      c=uisetcolor(get(s,'facecolor'));
      set(s,'facecolor',c);
     case 'isoalpha'
      [a s]=getArrowSlice;
      if nargin ~= 2
	error('Not enough arguments to sliceomatic.');
      end
      set(s,'facealpha',eval(p2));
     case 'isocaps'
      [a s]=getArrowSlice;
      cap=getappdata(s,'isosurfacecap');
      if isempty(cap)
	new=localisocaps(s);
	set(new,'uicontextmenu',d.uiciso);
      else
	delete(cap);
	setappdata(s,'isosurfacecap',[]);
      end
      %
      % Now for slice context menu items
      %
     case 'contextmenu'
      [a s]=getArrowSlice;
      set([d.smfacet d.smflat d.sminterp d.smtex d.smtnone d.smtp5 ...
	   d.smtflat d.smtinterp d.smttex],'checked','off');
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
     case 'togglevisible'
      [a s]=getArrowSlice;
      switch get(s,'visible')
       case 'on'
	set(s,'visible','off');
	pushset(a,'facealpha',.2);
       case 'off'
	set(s,'visible','on');
	popset(a,'facealpha');
      end
     case 'setfaceted'
      [a s]=getArrowSlice;
      set(s,'edgec','k','facec','flat');
      if ischar(get(s,'facea')) & strcmp(get(s,'facea'),'texturemap')
	set(s,'facea','flat');
      end
      textureizeslice(s,'off');
     case 'setflat'
      [a s]=getArrowSlice;
      set(s,'edgec','n','facec','flat');
      if ischar(get(s,'facea')) & strcmp(get(s,'facea'),'texturemap')
	set(s,'facea','flat');
      end
      textureizeslice(s,'off');
     case 'setinterp'
      [a s]=getArrowSlice;
      set(s,'edgec','n','facec','interp');
      if ischar(get(s,'facea')) & strcmp(get(s,'facea'),'texturemap')
	set(s,'facea','interp');
      end
      textureizeslice(s,'off');
     case 'settexture'
      [a s]=getArrowSlice;
      set(s,'facecolor','texture','edgec','none');
      if ischar(get(s,'facea'))
	set(s,'facealpha','texturemap');
      end
      textureizeslice(s,'on');
     case 'setalphanone'
      [a s]=getArrowSlice;
      set(s,'facealpha',1);
     case 'setalphapoint5'
      [a s]=getArrowSlice;
      set(s,'facealpha',.5);
     case 'setalphaflat'
      [a s]=getArrowSlice;
      set(s,'facealpha','flat');
      if ischar(get(s,'facec')) & strcmp(get(s,'facec'),'texturemap')
	set(s,'facecolor','flat');
	textureizeslice(s,'off');
      end
     case 'setalphainterp'
      [a s]=getArrowSlice;
      set(s,'facealpha','interp');
      if ischar(get(s,'facec')) & strcmp(get(s,'facec'),'texturemap')
	set(s,'facecolor','interp');
	textureizeslice(s,'off');
      end
     case 'setalphatexture'
      [a s]=getArrowSlice;
      set(s,'facealpha','texturemap');
      if ischar(get(s,'facec'))
	set(s,'facecolor','texturemap');
	textureizeslice(s,'on');
      end
     case 'deleteslice'
      [a s]=getArrowSlice;
      delete(s);
      if prod(size(a))==1
	delete(getappdata(a,'arrowedge'));
      end
      delete(a);
      %
      % Menu All Slices
      %
     case 'allfacet'
      s=allSlices;
      set(s,'facec','flat','edgec','k');
      textureizeslice(s,'off');
     case 'allflat'
      s=allSlices;
      set(s,'facec','flat');
      textureizeslice(s,'off');
     case 'allinterp'
      s=allSlices;
      set(s,'facec','interp');
      textureizeslice(s,'off');
     case 'alltex'
      s=allSlices;
      set(s,'facec','texturemap');
      textureizeslice(s,'on');
     case 'alltnone'
      s=allSlices;
      set(s,'facea',1);
      textureizeslice(s,'off');
     case 'alltp5'
      s=allSlices;
      set(s,'facea',.5);
      textureizeslice(s,'off');
     case 'alltflat'
      s=allSlices;
      set(s,'facea','flat');
      textureizeslice(s,'off');
     case 'alltinterp'
      s=allSlices;
      set(s,'facea','interp');
      textureizeslice(s,'off');
     case 'allttex'
      s=allSlices;
      set(s,'facea','texturemap');
      textureizeslice(s,'on');
      %
      % Menu Defaults callbacks
      %
     case 'defaultmenu'
      set([d.dfacet d.dflat d.dinterp d.dtex d.dtnone d.dtflat d.dtinterp ...
	   d.dttex ], 'checked','off');
      switch d.defcolor
       case 'faceted'
	set(d.dfacet,'checked','on');
       case 'flat'
	set(d.dflat,'checked','on');
       case 'interp'
	set(d.dinterp,'checked','on');
       case 'texture'
	set(d.dtex,'checked','on');
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
     case	'defaultfaceted'
      d.defcolor='faceted';
     case	'defaultflat'
      d.defcolor='flat';
     case	'defaultinterp'
      d.defcolor='interp';
     case	'defaulttexture'
      d.defcolor='texture';
      if strcmp(d.defalpha,'flat') | strcmp(d.defalpha,'interp')
	d.defalpha='texture';
      end
     case	'defaulttransnone'
      d.defalpha='none';
     case	'defaulttransflat'
      d.defalpha='flat';
     case	'defaulttransinterp'
      d.defalpha='interp';
     case	'defaulttranstexture'
      d.defalpha='texture';
      d.defcolor='texture';
     case      'defaultlightflat'
      d.deflight='flat';
     case      'defaultlightsmooth'
      d.deflight='smooth';
      %
      % UICONTROL callbacks
      %
     case 'colormap'
      str=get(gcbo,'string');
      colormap(str{get(gcbo,'value')});
     case 'alphamap'
      str=get(gcbo,'string');
      alphamap(str{get(gcbo,'value')});
     otherwise
      error('Bad slice-o-matic command.');
    end
  catch
    disp(get(0,'errormessage'));
  end
  setappdata(gcf,'sliceomatic',d);
end

function dragprep(arrowtodrag)

arrows=findall(gcf,'tag','sliceomaticarrow');

pushset(arrows,'facecolor','r');
pushset(arrows,'facealpha',.2);

pushset(arrowtodrag,'facecolor','g');
pushset(arrowtodrag,'facealpha',.7);

slices=allSlices;

for i=1:length(slices)
  fa=get(slices(i),'facea');
  if isa(fa,'double') & fa>.3
    pushset(slices(i),'facealpha',.3);
    pushset(slices(i),'edgecolor','n');
  else
    pushset(slices(i),'facealpha',fa);
    pushset(slices(i),'edgecolor',get(slices(i),'edgec'));
  end
end

isosurfs=allIsos;

for i=1:length(isosurfs)
  fa=get(isosurfs(i),'facea');
  if isa(fa,'double') & fa>.3
    pushset(isosurfs(i),'facealpha',.3);
    pushset(isosurfs(i),'edgecolor','n');
  else
    pushset(isosurfs(i),'facealpha',fa);
    pushset(isosurfs(i),'edgecolor',get(isosurfs(i),'edgec'));
  end
  cap=getappdata(isosurfs(i),'isosurfacecap');
  if ~isempty(cap)
    pushset(cap,'visible','off');
  end
end

ss=getappdata(arrowtodrag,'arrowslice');

if isempty(ss)
  ss=getappdata(arrowtodrag,'arrowiso');
end

popset(ss,'facealpha');
popset(ss,'edgecolor');

pushset(gcf,'windowbuttonupfcn','sliceomatic up');
pushset(gcf,'windowbuttonmotionfcn','sliceomatic motion');

function dragfinis(arrowtodrag)

arrows=findall(gcf,'tag','sliceomaticarrow');

popset(arrowtodrag,'facecolor');
popset(arrowtodrag,'facealpha');

popset(arrows,'facecolor');
popset(arrows,'facealpha');

ss=getappdata(arrowtodrag,'arrowslice');
if isempty(ss)
  ss=getappdata(arrowtodrag,'arrowiso');
end

% These pushes are junk which will be undone when all slices or
% isosurfs are reset below.
pushset(ss,'facealpha',1);
pushset(ss,'edgecolor','k');

slices=allSlices;

if ~isempty(slices)
  popset(slices,'facealpha');
  popset(slices,'edgecolor');
end

isosurfs=allIsos;

if ~isempty(isosurfs)
  popset(isosurfs,'facealpha');
  popset(isosurfs,'edgecolor');
end

d=getappdata(gcf,'sliceomatic');

for i=1:length(isosurfs)
  cap=getappdata(isosurfs(i),'isosurfacecap');
  if ~isempty(cap)
    popset(cap,'visible');
    localisocaps(isosurfs(i),cap);
  end
  if getappdata(isosurfs(i), 'reduced')
    setappdata(isosurfs(i),'reduced',0);
    localisosurface({},d.data,d.smooth,...
		    getappdata(isosurfs(i),'isosurfacevalue'),...
		    isosurfs(i));
  end
end

popset(gcf,'windowbuttonupfcn');
popset(gcf,'windowbuttonmotionfcn');

% Make sure whatever buttonupfcn on the figure is run now to "turn
% off" whatever was going on before we got our callback on the
% arrow.

buf = get(gcf,'windowbuttonupfcn');
if ~strcmp(buf,'')
  eval(buf);
end

function [a, s]=getArrowSlice
% Return the Arrow and Slice based on the GCO
if isempty(getappdata(gco,'controlarrow')) & ...
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

function p=arrow(parent,dir,pos)

%   21012    21012      12345     12345
% 5  *-*   5   *     2   *     2   *  
% 4  | |   4  / \    1 *-*\    1  /*-*
% 3 ** **  3 ** **   0 |   *   0 *   |
% 2  \ /   2  | |   -1 *-*/   -1  \*-*
% 1   *    1  *-*   -2   *    -2   *  

switch dir
 case 'down'
  pts=[ 0 1; -2 3; -1 3; -1 5; 1 5; 1 3; 2 3 ];
 case 'up'
  pts=[ 0 5; 2 3; 1 3; 1 1; -1 1; -1 3; -2 3; ];
 case 'right'
  pts=[ 5 0; 3 -2; 3 -1; 1 -1; 1 1; 3 1; 3 2 ];
 case 'left'
  pts=[ 1 0; 3 2; 3 1; 5 1; 5 -1; 3 -1; 3 -2 ];
end

f=[1 2 7; 3 4 5; 3 5 6 ];

if pos(1)
  lim=get(parent,'xlim');
  fivep=abs(lim(1)-lim(2))/15/5;
  pts(:,1)=pts(:,1)*fivep+pos(1);
elseif pos(2)
  lim=get(parent,'ylim');
  fivep=abs(lim(1)-lim(2))/15/5;
  pts(:,2)=pts(:,2)*fivep+pos(2);
end

p(1)=patch('vertices',pts,'faces',1:size(pts,1),'facec','n','edgec','k',...
	   'linewidth',2,'hittest','off',...
	   'parent',parent);
p(2)=patch('vertices',pts,'faces',f,'facec','g','facea',.5,'edgec','n',...
	   'parent',parent,'tag','sliceomaticarrow');
setappdata(p(2),'arrowcenter',pos);
setappdata(p(2),'arrowedge',p(1));

function p=localisocaps(isosurface,isocap)

% Get relevant info from the isosurface.
if nargin<2 | ~strcmp(get(isocap,'visible'),'off')
  data=getappdata(isosurface,'isosurfacedata');
  caps=isocaps(data,getappdata(isosurface,'isosurfacevalue'));
end

if nargin==2
  if ~strcmp(get(isocap,'visible'),'off')
    set(isocap,caps);
  end
  p=isocap;
else
  p=patch(caps,'edgecolor','none','facecolor','flat',...
	  'facelighting','none',...
	  'tag','sliceomaticisocap');

  setappdata(p,'isosurface',isosurface);
  setappdata(isosurface,'isosurfacecap',p);
  
  d=getappdata(gcf,'sliceomatic');
  
  switch d.defcolor
   case 'faceted'
    set(p,'facec','flat','edgec','k');
   case 'flat'
    set(p,'facec','flat','edgec','n');
   case 'interp'
    set(p,'facec','interp','edgec','n');
   case 'texture'
    set(p,'facec','flat','edgec','n');
  end
  switch d.defalpha
   case 'none'
    set(p,'facea',1);
   case 'flat'
    set(p,'facea','flat');
   case 'interp'
    set(p,'facea','interp');
   case 'texture'
    set(p,'facea','flat');
  end    
end


function p=localisosurface(volume, data, datanormals, value, oldiso)

pushset(gcf, 'pointer','watch');

fv = isosurface(volume{:},data, value);
clim=get(gca,'clim');
cmap=get(gcf,'colormap');
clen=clim(2)-clim(1);
idx=floor((value-clim(1))*length(cmap)/clen);

if nargin==5
  set(oldiso,fv,'facecolor',cmap(idx,:));
  p=oldiso;
  cap=getappdata(p,'isosurfacecap');
  if ~isempty(cap)
    localisocaps(p,cap);
  end
else
  p=patch(fv,'edgecolor','none','facecolor',cmap(idx,:),...
	  'tag', 'sliceomaticisosurface');
  d=getappdata(gcf,'sliceomatic');
  switch d.deflight
   case 'flat'
    set(p,'facelighting','flat');
   case 'smooth'
    set(p,'facelighting','phong');
  end
  setappdata(p,'isosurfacecap',[]);
end

setappdata(p,'isosurfacevalue',value);
setappdata(p,'isosurfacedata',data);

reducepatch(p,10000);
isonormals(volume{:},datanormals,p);

popset(gcf,'pointer');

function s=localslice(data, X, Y, Z, oldslice)

s=[];
d=getappdata(gcf,'sliceomatic');

ds=size(data);

if ~isempty(X)
  xi=round(X);
  if xi > 0 & xi < ds(2)
    cdata=reshape(data(:,xi,:),ds(1),ds(3));
    [xdata ydata zdata]=meshgrid(X,1:ds(1),1:ds(3));
  else
    return
  end
elseif ~isempty(Y)
  yi=round(Y);
  if yi > 0 & yi < ds(1)
    cdata=reshape(data(yi,:,:),ds(2),ds(3));
    [xdata ydata zdata]=meshgrid(1:ds(2),Y,1:ds(3));
  else
    return    
  end
elseif ~isempty(Z)
  zi=round(Z);
  if zi > 0 & zi < ds(3)
    cdata=reshape(data(:,:,zi),ds(1),ds(2));
    [xdata ydata zdata]=meshgrid(1:ds(2),1:ds(1),Z);
  else
    return
  end
else
  error('Nothing was passed into LOCALSLICE.');
end

cdata=squeeze(cdata);
xdata=squeeze(xdata);
ydata=squeeze(ydata);
zdata=squeeze(zdata);

if nargin == 5
  % Recycle the old slice
  set(oldslice,'cdata',cdata,'alphadata',cdata, 'xdata',xdata, ...
	       'ydata',ydata, 'zdata',zdata);
  s=oldslice;
  %delete(news);
  if propcheck(s,'facec','texturemap')
    textureizeslice(s,'on');
  end
  
else
  % setup the alphadata
  news=surface('cdata',cdata,'alphadata',cdata, 'xdata',xdata, ...
	       'ydata',ydata, 'zdata',zdata);
  set(news,'alphadata',cdata,'alphadatamapping','scaled','tag','sliceomaticslice',...
	   'facelighting','none',...
	   'uicontextmenu',d.uic);
  s=news;
  switch d.defcolor
   case 'faceted'
    set(s,'facec','flat','edgec','k');
   case 'flat'
    set(s,'facec','flat','edgec','n');
   case 'interp'
    set(s,'facec','interp','edgec','n');
   case 'texture'
    set(s,'facec','texture','edgec','n');
  end
  switch d.defalpha
   case 'none'
    set(s,'facea',1);
   case 'flat'
    set(s,'facea','flat');
   case 'interp'
    set(s,'facea','interp');
   case 'texture'
    set(s,'facea','texture');
  end    

  if strcmp(d.defcolor,'texture')
    textureizeslice(s,'on');
  end
end


function textureizeslice(slice,onoff)

for k=1:prod(size(slice))

  d=getappdata(slice(k),'textureoptimizeations');

  switch onoff
   case 'on'
    d.xdata=get(slice(k),'xdata');
    d.ydata=get(slice(k),'ydata');
    d.zdata=get(slice(k),'zdata');
    setappdata(slice(k),'textureoptimizeations',d);
    if max(size(d.xdata)==1)
      nx=[d.xdata(1) d.xdata(end)];
    else
      nx=[d.xdata(1,1)   d.xdata(1,end);
	  d.xdata(end,1) d.xdata(end,end)];
    end
    if max(size(d.ydata)==1)
      ny=[d.ydata(1) d.ydata(end)];
    else
      ny=[d.ydata(1,1)   d.ydata(1,end);
	  d.ydata(end,1) d.ydata(end,end)];
    end
    if max(size(d.zdata)==1)
      nz=[d.zdata(1) d.zdata(end)];
    else
      nz=[d.zdata(1,1)   d.zdata(1,end);
	  d.zdata(end,1) d.zdata(end,end)];
    end
    set(slice(k),'xdata',nx, 'ydata', ny, 'zdata', nz,...
		 'facec','texturemap');
    if ischar(get(slice(k),'facea'))
      set(slice(k),'facea','texturemap');
    end
    if ischar(get(slice(k),'facec'))
      set(slice(k),'facec','texturemap');
    end
   case 'off'
    if ~isempty(d)
      set(slice(k),'xdata',d.xdata,'ydata',d.ydata,'zdata',d.zdata);
      setappdata(slice(k),'textureoptimizeations',[]);
    end
    if ischar(get(slice(k),'facea')) & strcmp(get(slice(k),'facea'),'texturemap')
      set(slice(k),'facea','flat');
    end
    if ischar(get(slice(k),'facec')) & strcmp(get(slice(k),'facec'),'texturemap')
      set(slice(k),'facec','flat');
    end
  end
end

function tf=propcheck(obj, prop, value)

v=get(obj,prop);

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

function ss=allSlices
ss=findobj(gcf,'type','surface','tag','sliceomaticslice');

function ss=allIsos
ss=findobj(gcf,'type','patch','tag','sliceomaticisosurface');

function ss=allCaps
ss=findobj(gcf,'type','patch','tag','sliceomaticisocap');

function working(onoff)

ax=getappdata(gcf,'workingaxis');

if isempty(ax)
  ax=axes('units','norm','pos',[.3 .4 .4 .2],...
	  'box','on','ytick',[],'xtick',[],...
	  'xlim',[-1 1],'ylim',[-1 1],...
	  'color','none','handlevis','off');
  text('parent',ax,'string','Working...','fontsize',64,...
       'pos',[0 0], ...
       'horizontalalignment','center',...
       'verticalalignment','middle',...
       'erasemode','xor');
  setappdata(gcf,'workingaxis',ax);
end

disp(['Working...' onoff]);
set([ax get(ax,'children')],'vis',onoff);

%
% Property Stack Stuff
%
function pushset(handle,prop,value)
% PUSHSET - push new properties onto a value stack.
%
% PUSHSET(HANDLE, PROP, VALUE) will take the old value of PROP for
% HANDLE, and save it on a stack associated with HANDLE.  It then
% assigns VALUE as the new value.
%

nargchk(3,3,'wrong number of arguments.');

proplist=fieldnames(get(handle(1)));
prop=proplist{strcmpi(prop,proplist)};
appstr = [prop '_hgstack'];

for k=1:prod(size(handle))

  oldv = get(handle(k),prop);
  olds = getappdata(handle(k),appstr);
  set(handle(k),prop,value);
  setappdata(handle(k),appstr,{ oldv olds });
  
end

function popset(handle,prop)
% POPSET - pop values for a property from a value stack.
%
% POPSET(HANDLE, PROP) will restore a prevously HGPUSHED property value.
%


nargchk(2,2,'wrong number of arguments.');

proplist=fieldnames(get(handle(1)));
prop=proplist{strcmpi(prop,proplist)};

appstr = [prop '_hgstack'];

for k=1:prod(size(handle))
  
  olds = getappdata(handle(k),appstr);

  if length(olds) <= 1
    error(['Nothing left to pop for property ' prop '.']);
  end

  set(handle(k),prop,olds{1});
  setappdata(handle(k),appstr,olds{2:end});

end
