function[u,v,w]=sphere2uvw(lat,lon,dr,dth,dphi)
%SPHERE2UVW Converts a 3D spherical vector to a 3D Cartesian vector.
%
%   [U,V,W]=XYZ2SPHERE(LAT,LON,V1,V2,V3) converts a vector in spherical
%   coordinates with components V1, V2, and V3 located at point 
%   (LAT, LON) into a Cartesian 3-vector with components U, V, and W.
%
%   LAT and LON are in degrees.
%
%   The vector in the spherical coordinate system has radial component
%   V1, longitudinal component V2, and latitudinal component V3.   
%
%   The Cartesian vector [U,V,W] is in a reference frame with the
%   X-axis at zero degrees longitude and the Z-axis at the North Pole.  
%
%   All input arguments should be arrays of the same size.
%
%   SPHERE2UVW is inverted by UVW2SPHERE.
%
%   See JSPHERE for related functions.
%
%   Usage: [x,y,z]=sphere2uvw(lat,lon,v1,v2,v3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2014 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(lat, '--t')
    sphere2uvw_test,return
end
[phi,theta]=jdeg2rad(lat,lon);

u      =  cos(phi).*cos(theta).*dr  -sin(theta).*dth -  sin(phi).*cos(theta).*dphi;
v      =  cos(phi).*sin(theta).*dr + cos(theta).*dth -  sin(phi).*sin(theta).*dphi ;
w      =  sin(phi).*dr                                + cos(phi).*dphi;


function[]=sphere2uvw_test
 
lat=[0  0  45         0        -90]';
lon=[0  90 0         90         0 ]';
v1= [1  1  1          0         0]'; 
v2= [0  0  0         -1         0]';
v3= [0  0  0          0         1]';
u=  [1  0  sqrt(2)/2  1         1 ]';
v=  [0  1  0          0         0]';
w=  [0  0  sqrt(2)/2  0         0]';

[u2,v2,w2]=sphere2uvw(lat,lon,v1,v2,v3);
tol=1e-6;
reporttest('SPHERE2XYZ example points',aresame(u,u2,tol)&&aresame(v,v2,tol)&&aresame(w,w2,tol))
