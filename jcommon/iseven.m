function[bool]=iseven(x)
%ISEVEN  True for even integer values; false otherwise. 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2001--2015 J.M. Lilly --- type 'help jlab_license' for details
  
bool=floor(x/2)==x/2;
