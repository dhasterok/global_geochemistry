function setpointer(fig, ptr)
% Set the pointer on the current figure to PTR
% has several specialized SOM (SliceOMatic) pointers
  
% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  switch ptr
   case 'SOM left'
    pd = [ nan nan nan nan 1   nan nan nan nan nan nan nan nan nan nan nan
           nan nan nan 1   1   nan nan nan nan nan nan nan nan nan nan nan
           nan nan nan 1   1   nan nan nan nan nan nan nan nan nan nan nan
           nan nan 1   2   1   nan nan nan nan nan nan nan nan nan nan nan
           nan nan 1   2   1   1   1   1   1   1   1   1   1   1   1   1
           nan 1   2   2   2   2   2   2   2   2   2   2   2   2   2   1
           nan 1   2   2   2   2   2   2   2   2   2   2   2   2   2   1
           1   2   2   2   2   2   2   2   2   2   2   2   2   2   2   1
           1   2   2   2   2   2   2   2   2   2   2   2   2   2   2   1
           nan 1   2   2   2   2   2   2   2   2   2   2   2   2   2   1
           nan 1   2   2   2   2   2   2   2   2   2   2   2   2   2   1
           nan nan 1   2   1   1   1   1   1   1   1   1   1   1   1   1
           nan nan 1   2   1   nan nan nan nan nan nan nan nan nan nan nan
           nan nan nan 1   1   nan nan nan nan nan nan nan nan nan nan nan
           nan nan nan 1   1   nan nan nan nan nan nan nan nan nan nan nan
           nan nan nan nan 1   nan nan nan nan nan nan nan nan nan nan nan ];
    set(fig,'pointershapecdata', pd,...
            'pointershapehotspot', [ 8 1 ] , ...
            'pointer','custom');
   case 'SOM right'
    pd = [ nan nan nan nan nan nan nan nan nan nan nan 1  nan nan nan nan
           nan nan nan nan nan nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan nan nan nan nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan nan nan nan nan nan nan nan nan nan 1  2   1   nan nan 
           1   1   1   1   1   1   1   1   1   1   1   1  2   1   nan nan 
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   2   1   
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   2   1   
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           1   1   1   1   1   1   1   1   1   1   1   1  2   1   nan nan 
           nan nan nan nan nan nan nan nan nan nan nan 1  2   1   nan nan 
           nan nan nan nan nan nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan nan nan nan nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan nan nan nan nan nan nan nan nan nan 1  nan nan nan nan ];
    set(fig,'pointershapecdata', pd,...
            'pointershapehotspot', [ 8 16 ] , ...
            'pointer','custom');
   case 'SOM bottom'
    pd = [ nan nan nan nan 1   1   1   1 1 1   1   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           1   1   1   1   1   2   2   2 2 2   2   1   1   1   1   1  
           nan 1   1   2   2   2   2   2 2 2   2   2   2   1   1   nan
           nan nan nan 1   1   2   2   2 2 2   2   1   1   nan nan nan
           nan nan nan nan nan 1   1   2 2 1   1   nan nan nan nan nan
           nan nan nan nan nan nan nan 1 1 nan nan nan nan nan nan nan ];

    set(fig,'pointershapecdata', pd,...
            'pointershapehotspot', [ 16 8 ] , ...
            'pointer','custom');
   
   case 'SOM top'
    pd = [ nan nan nan nan nan nan nan 1 1 nan nan nan nan nan nan nan
           nan nan nan nan nan 1   1   2 2 1   1   nan nan nan nan nan
           nan nan nan 1   1   2   2   2 2 2   2   1   1   nan nan nan
           nan 1   1   2   2   2   2   2 2 2   2   2   2   1   1   nan
           1   1   1   1   1   2   2   2 2 2   2   1   1   1   1   1  
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   1   1   1 1 1   1   1   nan nan nan nan ];

    set(fig,'pointershapecdata', pd,...
            'pointershapehotspot', [ 1 8 ] , ...
            'pointer','custom');
   
   case 'SOM leftright'
    pd = [ nan nan nan nan 1   nan nan nan nan nan nan 1  nan nan nan nan
           nan nan nan 1   1   nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan nan 1   1   nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan 1   2   1   nan nan nan nan nan nan 1  2   1   nan nan 
           nan nan 1   2   1   1   1   1   1   1   1   1  2   1   nan nan 
           nan 1   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           nan 1   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   2   1   
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   2   1   
           nan 1   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           nan 1   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           nan nan 1   2   1   1   1   1   1   1   1   1  2   1   nan nan 
           nan nan 1   2   1   nan nan nan nan nan nan 1  2   1   nan nan 
           nan nan nan 1   1   nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan nan 1   1   nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan nan nan 1   nan nan nan nan nan nan 1  nan nan nan nan ];
    set(fig,'pointershapecdata', pd,...
            'pointershapehotspot', [ 8 8 ] , ...
            'pointer','custom');
   
   case 'SOM topbottom'
    pd = [ nan nan nan nan nan nan nan 1 1 nan nan nan nan nan nan nan
           nan nan nan nan nan 1   1   2 2 1   1   nan nan nan nan nan
           nan nan nan 1   1   2   2   2 2 2   2   1   1   nan nan nan
           nan 1   1   2   2   2   2   2 2 2   2   2   2   1   1   nan
           1   1   1   1   1   2   2   2 2 2   2   1   1   1   1   1  
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           1   1   1   1   1   2   2   2 2 2   2   1   1   1   1   1  
           nan 1   1   2   2   2   2   2 2 2   2   2   2   1   1   nan
           nan nan nan 1   1   2   2   2 2 2   2   1   1   nan nan nan
           nan nan nan nan nan 1   1   2 2 1   1   nan nan nan nan nan
           nan nan nan nan nan nan nan 1 1 nan nan nan nan nan nan nan ];


    set(fig,'pointershapecdata', pd,...
            'pointershapehotspot', [ 8 8 ] , ...
            'pointer','custom');
    
   otherwise
    % Set it to the string passed in
     set(fig,'pointer', ptr);
  end
