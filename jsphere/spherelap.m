function[lap]=spherelap(varargin)
%SPHERELAP  Laplacian of a field on the surface of a sphere.
%
%   DEL2=SPHERELAT(LAT,LON,F) computes the Laplacian of the scalar field F
%   on the surface of the sphere.
%
%   DEL2 has units of F per inverse meter squared.
%
%   LAT and LON are vectors specifing an evenly-spaced grid, and F is an 
%   array of size LENGTH(LON) x LENGTH(LAT) x M, where M is greater than
%   or equal to one.  LAT and LON are in degrees.
%
%   The radius of the Earth as specified by RADEARTH is used by default.
%   SPHERELAP(...,R) uses a sphere of radius R, in kilometers, instead.
%
%   Derivatives are computed using the second central difference.
%   ___________________________________________________________________
%
%   First and last points
%
%   SPHERELAP can use different boundary conditions in the numerical
%   computation of derivatives for the first and last points.
% 
%   SPHERELAP(...,'periodic') uses a periodic derivative with respect
%   to longitude.  This is the default behavior.  
%
%   SPHERELAP(...,'endpoint') uses the forwards / second backwards 
%   difference at the first and last longitude, respectively. 
%
%   For both of these options, the endpoint condition is used for 
%   differentiation with respect to latitude.
%
%   SPHERELAP(...,'nans')  fills in the first and last values on all four
%   sides of the region with NANs. 
%
%   See VDIFF for more information.  
%   ___________________________________________________________________
%
%   See also SPHEREDIV, SPHEREGRAD, SPHERECURL, JSPHERE.
%
%   'spherelap --t' runs a test.
%
%   Usage: del2=spherelap(lat,lon,f);
%          del2=spherelap(lat,lon,f,R);
%          del2=spherelap(lat,lon,f,'endpoint');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2009 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    spherelap_test,return
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

if length(varargin)==3
    R=radearth;
else
    R=varargin{4};
end

lat=varargin{1};
lon=varargin{2};
f=varargin{3};

R=R*1000;  %Convert to meters

[phi,theta]=jdeg2rad(lat,lon);
dphi=phi(2)-phi(1);
dth=theta(2)-theta(1);
[lon,lat]=meshgrid(lon,lat);
[phi,theta]=jdeg2rad(lat,lon);

if ~aresame(size(lon),[size(f,1),size(f,2)])
    error('F must be oriented with longitude in columns and latitude in rows.')
end

if dphi<0
    error('LAT should be ordered as an increasing array.')
end
    
if size(f,3)>1
    vrep(lon,lat,size(f,3),3);
end

lapa=frac(1,R.^2.*cos(phi).^2).*frac(1,dphi.^2).*vdiff(vdiff(f,1,strlat),1,strlat);
lapb=frac(1,R.^2.*cos(phi)).*frac(1,dth.^2).*vdiff(cos(phi).*vdiff(f,2,str),2,str);

lap=lapa+lapb;
tol=1e-6;
index=find(abs(cos(phi))<tol);
if ~isempty(index)
    lap(index)=nan;
end

function[]=spherelap_test
 
lon=(0:2:358)-180;
lat=(-90:1:90);
[long,latg]=meshgrid(lon,lat);
f=randn(size(latg));

del2fa=spherelap(lat,lon,f,'nans');

[gradx,grady]=spheregrad(lat,lon,f,'nans');
del2fb=spherediv(lat,lon,gradx,grady,'nans');

reporttest('SPHERELAP matches SPHEREDIV of SPHEREGRAD',aresame(del2fa,del2fb,1e-6))

