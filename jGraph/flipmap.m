function[]=flipmap
%FLIPMAP Flips the current colormap upside-down.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details  
ZZmap=colormap;
colormap(flipud(ZZmap));
clear ZZmap
