function[h]=nocontours(h)
%NOCONTOURS  Removes contours from a CONTOURF plot.
%
%   NOCONTOURS(H) where H is an array of handles to lines, output by
%   CONTOURF, sets the linestyle for all the lines to 'none'.  This
%   removes the contours from the CONTOURF plot.
%
%   NOCONTOURS with no input arguments sets the linestyle of all patch
%   objects (for Matlab versions < 2014b) or all contour objects (for 
%   versions 2014b and later) in the current axes to 'none'.
%
%   H=NOCONTOURS also returns the handles to the objects.  
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2015 J.M. Lilly --- type 'help jlab_license' for details        


if nargin==0
   h=contourhandles(gca);
end

set(h,'linestyle','none');
    
if nargout ==0
  clear h
end

function[h]=contourhandles(axh)
%CONTOURHANDLES  Finds all contour handles from a given set of axes.
%
%   H=CONTOURHANDLES returns handles to all contoures associated with the
%   current axes.
%	 
%   H=CONTOURHANDLES(AX) returns handles to all contoures associated with
%   the set of axes whose handle is AX.
%
%   See also LINEHANDLES, PATCHHANDLES, AXESHANDLES.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014 J.M. Lilly --- type 'help jlab_license' for details        
  
if nargin==0
   axh=gca;
end


h=get(axh,'children');
bool=false(size(h));
for j=1:length(h)
    if strcmpi(get(h(j),'type'),'contour')
       bool(j)=true;
    end
end
h=h(bool);     
