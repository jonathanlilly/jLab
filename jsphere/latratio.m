function[ar]=latratio(lat,h)
%LATRATIO  Set plot aspect ratio for latitude / longitude plot.
%
%   LATRATIO(LAT) sets the aspect ratio of the current axis correctly a 
%   Cartesian plot centered about the latitude LAT, i.e. the x/y aspect 
%   ratio is set to COS(LAT).
%
%   Equal distances along the x- and y-axes then correspond to the same 
%   physical distance.
%
%   LAT is measured in degrees.
%
%   LATRATIO(LAT,H) does the same for the axis with handle H.
%
%   LATRATIO with no input arguments sets the aspect ratio of the current
%   axis using the midpoint of the y-axis limits.
%
%   AR=LATRATIO(LAT) alternately returns the correct aspect ratio in a 
%   two-element array AR, without modifying any plots. 
%
%   Note when used with QUIVER, one should multiply the u-velocity input to
%   QUIVER by AR(1) in order to obtain a correct scaling.  For example,
%   
%      figure,axis([-5 5 55 65]);hold on
%      latratio(60);ar=latratio(60);quiver(0,60,ar(1)*1,1);
%  
%   produces an arrow oriented to the northeast, as expected. 
%
%   Usage: latratio(lat,h);
%          latratio(lat);
%          latratio;
%          ar=latratio(lat);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if nargin==0
   ax=axis;
   lat=frac(1,2)*(ax(3)+ax(4));
end

if nargin<2
    h=gca;
end

ar=[1./cos(jdeg2rad(lat)) 1 1];
if nargout==1
    ar=ar(1:2);
else
    set(h,'dataaspectratio',ar)
end
if nargout==0
    clear ar
end

