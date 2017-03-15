function[curl,uy,vx]=spherecurl(varargin)
%SPHERECURL  Curl of a vector field on the surface of a sphere.
%
%   CURL=SPHERECURL(LAT,LON,U,V) computes the vertical component of the 
%   curl of the vector field (U,V) on the surface of the sphere.
%
%   U and V are zonal and meridional velocities in meters per second, and
%   CURL is in inverse seconds.
%
%   LAT and LON are vectors specifing an evenly-spaced grid, and U and V
%   are arrays of size LENGTH(LAT) x LENGTH(LON) x M, where M is greater
%   than or equal to one. LAT and LON are in degrees.
%
%   The radius of the Earth as specified by RADEARTH is used by default.
%   SPHERECURL(...,R) uses a sphere of radius R, in kilometers, instead.
%
%   Derivatives are computed using the first central difference.
%
%   [CURL,UY,VX]=SPHERECURL(LAT,LON,U,V) optionally returns the zonal and
%   meridional contributions UY and VX, such that CURL=VX-UY.
%   ___________________________________________________________________
%
%   First and last points
%
%   SPHERECURL can use different boundary conditions in the numerical
%   computation of derivatives for the first and last points.
% 
%   SPHERECURL(...,'periodic') uses a periodic derivative with respect to
%   longitude.  This is the default behavior.  
%
%   SPHERECURL(...,'endpoint') uses the forwards / first backwards 
%   difference at the first and last longitude, respectively. 
%
%   For both of these options, the endpoint condition is used for 
%   differentiation with respect to latitude.
%
%   SPHERECURL(...,'nans')  fills in the first and last values on all four
%   sides of the region with NANs. 
%
%   See VDIFF for more information.  
%   ___________________________________________________________________
%
%   See also SPHEREGRAD, SPHERELAP, SPHERDIV, JSPHERE.
%
%   'spherecurl --t' runs a test.
%
%   Usage: curl=spherecurl(lat,lon,u,v);
%          [curl,uy,vx]=spherecurl(lat,lon,u,v);
%          curl=spherecurl(lat,lon,u,v,R);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2017 J.M. Lilly --- type 'help jlab_license' for details
 

%http://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates

if strcmpi(varargin{1}, '--t')
    spherecurl_test,return
end
 
if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
else
    str='periodic';
end

if strcmpi(str(1:3),'per')
    strlat='endpoint';
else
    strlat=str;
end

if length(varargin)==4
    R=radearth;
else
    R=varargin{5};
end

lat=varargin{1};
lon=varargin{2};
uh=varargin{3};
vh=varargin{4};

R=R*1000;

[phi,theta]=jdeg2rad(lat,lon);
dphi=phi(2)-phi(1);
dth=theta(2)-theta(1);
[lon,lat]=meshgrid(lon,lat);

if ~aresame(size(lon),[size(uh,1),size(uh,2)])
    error('U and V must be oriented with longitude in columns and latitude in rows.')
end
if ~aresame(size(uh),size(vh))
    error('U and V must be the same size.')
end
if dphi<0
    error('LAT should be ordered as an increasing array.')
end
    
if size(uh,3)>1
    vrep(lon,lat,size(uh,3),3);
end

[phi,theta]=jdeg2rad(lat,lon);
%div=frac(1,R.*cos(phi).*dphi).*vdiff(vh.*cos(phi),1,strlat)+frac(1,R.*cos(phi).*dth).*vdiff(uh,2,str);
vx=frac(1,R.*cos(phi).*dth).*vdiff(vh,2,str);
uy=frac(1,R.*cos(phi).*dphi).*vdiff(uh.*cos(phi),1,strlat);
curl=vx-uy;


function[]=spherecurl_test
lon=(0:2:358)-180;
lat=(-90:1:90);
[long,latg]=meshgrid(lon,lat);

f=randn(size(long));
[gradx,grady]=spheregrad(lat,lon,f);
curl=spherecurl(lat,lon,gradx,grady);


tol=1e-6;
reporttest('SPHERECURL curl grad vanishes',maxmax(abs(curl))<tol)
