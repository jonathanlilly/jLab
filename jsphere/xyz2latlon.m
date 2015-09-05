function[lat,lon]=xyz2latlon(x,y,z,R)
%XYZ2LATLON  Converts 3D Cartesian coordinates into latitude and longitude.
%
%   [LAT,LON]=XYZ2LATLON(X,Y,Z) converts Cartesian coordinates of a
%   vector in three dimensions into latitude and longitude. 
%
%   LAT and LON are in degrees, and X, Y, and Z can be in any units.  The 
%   output LON will be in the range [-180,180].
%
%   The length of the vector (X,Y,Z) does not matter since X, Y, and Z are 
%   rescaled in order to lie on the unit sphere.
%
%   The Cartesian coordinate system is a right-handed system whose origin
%   lies at the center of the sphere.  It is oriented with the Z-axis 
%   passing though the poles and the X-axis passing through LAT=0, LON=0.  
%
%   XYZ2LATLON is inverted by LATLON2XYZ.
%
%   See JSPHERE for related functions.
%
%   'xyz2latlon --t' runs a test.
%
%   Usage: [lat,lon]=xyz2latlon(x,y,z);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2014 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(x, '--t')
    xyz2latlon_test,return
end

R=sqrt(abs(x).^2+abs(y).^2+abs(z).^2);
x=x./R;
y=y./R;
z=z./R;

phi=asin(z);
th=imlog(x+sqrt(-1).*y);

[lat,lon]=jrad2deg(phi,th);

function[]=xyz2latlon_test
 
lon=(1e-10:2:360)-180;
lat=(-89:1:89);
[lon,lat]=meshgrid(lon,lat);

[x,y,z]=latlon2xyz(lat,lon);
[lat1,lon1]=xyz2latlon(x,y,z);

reporttest('XYZ2LATLON inverts LATLON2XY',aresame(lon1,lon,1e-10) && aresame(lat1,lat,1e-10))


