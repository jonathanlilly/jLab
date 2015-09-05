function[x]=nonnan(x)
%NONNAN  Return all non-NAN elements of an array.
%
%   NONNAN(X) is equivalent to X(FIND(~ISNAN(X))).
%   The empty set is returned if all elements of X are NANS.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2007 J.M. Lilly --- type 'help jlab_license' for details        

index=find(~isnan(x)); 
if ~isempty(index)
  x=x(index);
else
  x=[];
end
