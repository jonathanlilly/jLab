function map = haxby(m)
%HAXBY  The Haxby colormap.
%
%   HAXBY(M) returns an M-by-3 matrix containing the Haxby colormap, see
%
%        https://www.mathworks.com/matlabcentral/fileexchange/
%                       25690-haxby-color-map?s_tid=prof_contriblnk
%
%   HAXBY with no arguments returns a colormap having the same number of 
%   colors as the colormap of the current figure.
%
%   To make HAXBY your default colormap, add to your startup.m file the
%   line "set(0,'DefaultFigureColormap',haxby)".
%
%   This implementation was written by Kelsey Jordahl, modified and 
%   redistributed in accordance with the copyright policies given in 
%   HAXBY_COPYRIGHT.
%
%   Kelsey's original notes:
%
%   Colormap is based on the colors used by W. F. Haxby's Gravity Field of
%   World's oceans, 1985, developed for geoid and gravity maps.  The 
%   version used here is formed from a linear interpolation of the GMT 
%   color table used by MB-System by David W. Caress and Dale N. Chayes.
%
%        https://www.mbari.org/products/research-software/mb-system/
%
%   Usage: h=haxby(M);
%          colormap haxby
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2020 J.M. Lilly --- type 'help jlab_license' for details
 

if nargin < 1, m = size(get(gcf,'colormap'),1); end
% mbm_grdplot Haxby color pallette
ncolors=11;
c=[ 37    57   175;    40   127   251;    50   190   255;   106   235   255;
    138   236   174;   205   255   162;   240   236   121;   255   189    87;
    255   161    68;   255   186   133;   255   255   255];
pp=1:(m-1)/(ncolors-1):m;
r=interp1(pp,c(:,1),1:m);
g=interp1(pp,c(:,2),1:m);
b=interp1(pp,c(:,3),1:m);
map=[r' g' b']/255;
