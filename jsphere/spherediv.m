function[div,ux,vy]=spherediv(varargin)
%SPHEREDIV  Divergence of a vector field on the surface of a sphere.
%
%   DIV=SPHEREDIV(LAT,LON,U,V) computes the divergence of the vector field
%   (U,V) on the surface of the sphere.
%
%   U and V are zonal and meridional velocities in meters per second, and
%   DIV is in inverse seconds.
%
%   LAT and LON are vectors specifing an evenly-spaced grid, and U and V 
%   are arrays of size LENGTH(LON) x LENGTH(LAT) x M, where M is greater
%   than or equal to one. LAT and LON are in degrees.
%
%   The radius of the Earth as specified by RADEARTH is used by default.
%   SPHEREDIV(...,R) uses a sphere of radius R, in kilometers, instead.
%
%   Derivatives are computed using the first central difference.
%
%   [DIV,UX,VY]=SPHEREDIV(LAT,LON,U,V) optionally returns the zonal and
%   meridional contributions UX and VY, such that DIV=UX+VY.
%   ___________________________________________________________________
%
%   First and last points
%
%   SPHEREDIV can use different boundary conditions in the numerical
%   computation of derivatives for the first and last points.
% 
%   SPHEREDIV(...,'periodic') uses a periodic derivative with respect
%   to longitude.  This is the default behavior.  
%
%   SPHEREDIV(...,'endpoint') uses the forwards / first backwards 
%   difference at the first and last longitude, respectively. 
%
%   For both of these options, the endpoint condition is used for 
%   differentiation with respect to latitude.
%
%   SPHEREDIV(...,'nans')  fills in the first and last values on all four
%   sides of the region with NANs. 
%
%   See VDIFF for more information.  
%   ___________________________________________________________________
%
%   See also SPHEREGRAD, SPHERELAP, SPHERECURL, JSPHERE.
%
%   'spherediv --t' runs a test.
%
%   Usage: div=spherediv(lat,lon,u,v);
%          [div,ux,vy]=spherediv(lat,lon,u,v);
%          div=spherediv(lat,lon,u,v,R);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2009 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    spherediv_test,return
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
vy=frac(1,R.*cos(phi).*dphi).*vdiff(vh.*cos(phi),1,strlat);
ux=frac(1,R.*cos(phi).*dth).*vdiff(uh,2,str);
div=ux+vy;

%DIV = 1 / (R cos(phi) dphi)  * d/dphi vphi*cos(phi) + 1 / (R cos(phi) dth)* d/dth vth
%Longer way, but still correct -- I put this as a test into hor2xyz
%[u,v,w]=hor2xyz(lat,lon,uh,vh);
%[v1,v2,v3]=xyz2sphere(lat,lon,u,v,w);
%div=frac(1,R.*cos(phi).*dphi).*vdiff(v3.*cos(phi),1,strlat)+frac(1,R.*cos(phi).*dth).*vdiff(v2,2,str);

function[]=spherediv_test
 
lon=(1e-10:2:360)-180;
lat=(-90:1:90);
[long,latg]=meshgrid(lon,lat);

u=1+0*latg;
v=0*latg;

div=spherediv(lat,lon,u,v,'nans');
tol=1e-10;
reporttest('SPHEREDIV divergenceless eastward velocity field',maxmax(abs(div))<tol)

[phi,theta]=jdeg2rad(latg,long);
u=0*latg;
v=1./cos(phi);

div=spherediv(lat,lon,u,v,'nans');
tol=1e-10;
reporttest('SPHEREDIV divergenceless northward velocity field',maxmax(abs(div))<tol)

u=0*latg;
v=1+0*latg;

[phi,theta]=jdeg2rad(lat,lon);
dphi=phi(2)-phi(1);
[phi,theta]=jdeg2rad(latg,long);
div2=frac(1,1000*radearth.*cos(phi).*dphi).*vdiff(v.*cos(phi),1,'nans');

div=spherediv(lat,lon,u,v,'nans');
tol=1e-12;
reporttest('SPHEREDIV constant northward velocity field',maxmax(abs(div2-div))<tol)
