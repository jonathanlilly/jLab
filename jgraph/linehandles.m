function[h]=linehandles(axh)
%LINEHANDLES  Finds all line and patch handles from a given set of axes.
%
%   H=LINEHANDLES returns handles to all lines associated with the current
%   axes.
%	 
%   H=LINEHANDLES(AX) returns handles to all lines associated with the set
%   of axes whose handle is AX.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2015 J.M. Lilly --- type 'help jlab_license' for details        
  
if nargin==0
   axh=gca;
end

h=get(axh,'children');
bool=false(size(h));
for j=1:length(h)
    if strcmpi(get(h(j),'type'),'line')
       bool(j)=true;
    end
end
h=h(bool);



