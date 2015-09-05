function[h]=axeshandles(fignum)
%AXESHANDLES  Returns handles to all axes children.
%
%   AXESHANDLES returns a vector of handles to all axes childen
%   regardless of which figure is their parent.
%
%   AXESHANDLES(FIGNUM) returns a vector of handles to axes children of
%   Figure FIGNUM.
%
%   See also ALLCHILD, FINDALL, LINEHANDLES, CONTOURHANDLES, PATCHHANDLES.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 1998, 2004 J.M. Lilly --- type 'help jlab_license' for details

if nargin==1
	h=get(fignum,'children');
else
  	figs=get(0,'children');
	for i=1:length(figs)
		h=[h;get(figs(i),'children')];
	end
end
bool=false(size(h));
for i=1:length(h)
	if strcmpi(get(h(i),'type'),'axes')
		bool(i)=1;
	end
end
h=h(bool); 

