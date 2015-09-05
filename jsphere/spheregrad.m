function[gradx,grady]=spheregrad(varargin)
%SPHEREGRAD  Gradient of a field on the surface of a sphere.
%
%   [FX,FY]=SPHEREGRAD(LAT,LON,F) computes the gradient of the scalar field
%   F on the surface of the sphere.
%
%   FX and FY are the components of the gradient of F in the zonal and
%   meridional directions, respectively, with units of F per meter.
%
%   LAT and LON are vectors specifing an evenly-spaced grid, and F is an 
%   array of size LENGTH(LON) x LENGTH(LAT) x M, where M is greater than or
%   equal to one.  LAT and LON are in degrees.
%
%   The radius of the Earth as specified by RADEARTH is used by default.
%   SPHEREDIV(...,R) uses a sphere of radius R, in kilometers, instead.
%
%   Derivatives are computed using the first central difference.
%
%   FZ=SPHEREGRAD(LAT,LON,F) with one output argument returns the 
%   complex-valued array FZ=FX+i*FY.
%   ___________________________________________________________________
%
%   First and last points
%
%   SPHEREGRAD can use different boundary conditions in the numerical
%   computation of derivatives for the first and last points.
% 
%   SPHEREGRAD(...,'periodic') uses a periodic derivative with respect to
%   longitude.  This is the default behavior.  
%
%   SPHEREGRAD(...,'endpoint') uses the forwards / first backwards 
%   difference at the first and last longitude, respectively. 
%
%   For both of these options, the endpoint condition is used for 
%   differentiation with respect to latitude.
%
%   SPHEREGRAD(...,'nans')  fills in the first and last values on all four
%   sides of the region with NANs. 
%
%   See VDIFF for more information.  
%   ___________________________________________________________________
%
%   See also SPHEREDIV, SPHERELAP, SPHERECURL, JSPHERE.
%
%   'spheregrad --t' runs a test.
%
%   Usage: [fx,fy]=spheregrad(lat,lon,f);
%          [fx,fy]=spheregrad(lat,lon,f,R);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2009 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    spheregrad_test,return
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

if ~aresame(size(lon),[size(f,1),size(f,2)])
    error('F must be oriented with longitude in columns and latitude in rows.')
end

if dphi<0
    error('LAT should be ordered as an increasing array.')
end
    
if size(f,3)>1
    vrep(lon,lat,size(f,3),3);
end
[phi,theta]=jdeg2rad(lat,lon);

gradx=frac(1,R.*cos(phi).*dth).*vdiff(f,2,str);
grady=frac(1,R.*dphi).*vdiff(f,1,strlat);

tol=1e-6;
index=find(abs(cos(phi))<tol);
if ~isempty(index)
    gradx(index)=nan;
end

if nargout==1
    gradx=gradx+sqrt(-1)*grady;
end

function[]=spheregrad_test
 
lon=(0:2:358)-180;
lat=(-90:1:90);
[long,latg]=meshgrid(lon,lat);

f=latg*1e6;

df=2*1e6/(2*2*pi/360)/(1000*radearth);
 
[gradx,grady]=spheregrad(lat,lon,f,'nans');

tol=1e-6;
reporttest('SPHEREGRAD purely meridional gradient',maxmax(abs(grady-df))<tol)

f=long*1e6;
[gradx,grady]=spheregrad(lat,lon,f,'nans');
df=2*1e6./(2*2*pi/360)./(1000*radearth.*cos(jdeg2rad(latg)));
tol=1e-6;
reporttest('SPHEREGRAD purely zonal gradient',maxmax(abs(gradx-df))<tol)


f=latg.^2+cos(jdeg2rad(long)).^3;
u=1+0*latg;
v=1+0*latg;
div=spherediv(lat,lon,u,v);
[gradx,grady]=spheregrad(lat,lon,f,'nans');

x1=spherediv(lat,lon,f.*u,f.*v,'nans');
x2=f.*div+u.*gradx+v.*grady;

tol=1e-5;
reporttest('SPHEREGRAD, SPHEREDIV flux identity',aresame(x1,x2,tol))


