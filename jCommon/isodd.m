function[bool]=isodd(x)
%ISODD  True for odd integer values; false otherwise. 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002--2015 J.M. Lilly --- type 'help jlab_license' for details
bool=iseven(x+1);
