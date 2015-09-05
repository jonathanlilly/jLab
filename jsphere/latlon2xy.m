function[varargout]=latlon2xy(varargin)
% LATLON2XY  Converts latitude and longitude into local Cartesian coordinates.
%
%   [X,Y]=LATLON2XY(LAT,LON,LATO,LONO) converts (LAT,LON) with units of
%   degrees into displacements (X,Y) in a plane tangent to the earth at the
%   point (LATO, LONO). X and Y have units of kilometers.
%
%   CX=LATLON2XY(LAT,LON,LATO,LONO) with one output argument returns the 
%   location as a complex-valued quantity X+SQRT(-1)*Y. NANs in LAT or LON
%   become NAN+SQRT(-1)*NAN.
%   
%   LAT and LON are arrays of the same size.  LATO and LONO are either also
%   arrays of this size, or else scalars.  X and Y have the same size as 
%   the input arrays.  
%
%   X and Y are computed by projecting the tangent plane onto the sphere
%   using full spherical geometry.  
%
%   The radius of the earth is given by the function RADEARTH.
%
%   Note that X and Y are set to NANs for points on the opposite side of
%   the earth from the tangent plane, that is, where the great circle 
%   distance would exceed RADEARTH * pi/2.  
%   ___________________________________________________________________
%
%   Great circle distance
%
%   [X,Y,D]=LATLON2XY(...) also returns the great circle distance D between
%   the two sets of points.  
%
%   For points on the same side on the earth from the tangent plane, i.e.
%   where the great circle distance is less than RADEARTH * pi/2, LATLON2XY
%   gives the same distance as SPHEREDIST.  
%
%   However, for points on the opposite side of the earth, LATLON2XY 
%   returns NANs whereas SPHEREDIST returns the correct distance.
%
%   The great circle distance computed here is useful because it is a fast 
%   computation if X and Y are already known.
%   ___________________________________________________________________
%
%   Cell array input / output
%
%   LATLON2XY returns cell array output given cell array input.  
%
%   That is, if LAT, LON, LATO, and LONO are all cell arrays of length K, 
%   containing K different numerical arrays, then the output will also be 
%   cell arrays of length K.  
%
%   This also works if LAT and LON are cell arrays but LATO and LONO are 
%   scalars.
%   ___________________________________________________________________
%
%   LATLON2XY is inverted by XY2LATLON.
%
%   See also XY2LATLON, LATLON2UV, SPHEREDIST.
%
%   'latlon2xy --t' runs some tests.
%
%   Usage:  [x,y]=latlon2xy(lat,lon,lato,lono);
%           [x,y,d]=latlon2xy(lat,lon,lato,lono);
%           cx=latlon2xy(lat,lon,lato,lono);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2014 J.M. Lilly --- type 'help jlab_license' for details        
  
if strcmpi(varargin{1}, '--t')
  latlon2xy_test,return
end

na=nargin;
if ischar(varargin{end})
    str=varargin{end};
    na=na-1;
else 
    str='sphere';
end

lat=varargin{1};
lon=varargin{2};
lato=varargin{3};
lono=varargin{4};

R=radearth;

if ~iscell(lat)
    [x,y,d]=latlon2xy_celloop(nargout,lat,lon,lato,lono,R,str);
else
   if ~iscell(lato)
       lato1=lato;
       clear lato
       for i=1:length(lat)
           lato{i}=lato1;
       end
   end
   if ~iscell(lono)
       lono1=lono;
       clear lono
       for i=1:length(lat)
           lono{i}=lono1;
       end
   end

    for i=1:length(lat)
       [x{i},y{i},d{i}]=latlon2xy_celloop(nargout,lat{i},lon{i},lato{i},lono{i},R,str); 
    end
end

if nargout==1
    if iscell(x)
        varargout{1}=celladd(x,cellmult(sqrt(-1),y));
    else
        varargout{1}=x+sqrt(-1)*y;
    end
elseif nargout==2
    varargout{1}=x;
    varargout{2}=y;
elseif nargout==3
    varargout{1}=x;
    varargout{2}=y;
    varargout{3}=d;
end


function[x,y,d]=latlon2xy_celloop(N,lat,lon,lato,lono,R,str)

[x,y,d]=vempty;
if ~isempty(lat)&&~isempty(lato)
    if numel(lato)==1
        lato=lato+0*lat;
    end
    if numel(lono)==1
        lono=lono+0*lat;
    end
    
    if strcmpi(str(1:3),'sma')
        [x,y]=latlon2xy_cartesian(lat,lon-lono,lato,R);
    elseif strcmpi(str(1:3),'sph')
        if N<3
            [x,y]=latlon2xy_sphere(lat,lon-lono,lato,R);
        else
            [x,y,d]=latlon2xy_sphere(lat,lon-lono,lato,R);
        end
    end
end
    
function[x,y]=latlon2xy_cartesian(lat,lon,lato,R)

[lat,lon,lato]=jdeg2rad(lat,lon,lato);
r1=R.*cos(lato);
dlon=angle(rot(lon));
x=dlon.*r1;
y=(lat-lato).*R;

function[x,y,d]=latlon2xy_sphere(lat,lon,lato,R)

x=nan*lat;
y=nan*lon;

%  Check to see if I'm in the same hemisphere, or on the other side
%  Here I'm implementing latlon2xyz to find dot product of two vectors
%  and, I avoid recomputing angles

%if ~isreal(lato),lato,end
coslat=cosd(lat);
coslon=cosd(lon);
coslato=cosd(lato);
sinlato=sind(lato);
sinlat=sind(lat);

%  The next line is the dot product, i.e. x.*x2+z.*z2 with
%  x=cos(lat).*cos(lon);z=sin(lat);
%  x2=cos(lato);z2=sin(lato);

index=find(coslat.*coslon.*coslato+sinlato.*sinlat>0);

%  Now implement these equations
%    x=R*cos(lat).*sin(lon);
%    y=-R*cos(lat).*sin(lato).*cos(lon)+R.*cos(lato).*sin(lat);
%  but only for points in the same hemisphere.

%  To derive these equations, calculate the xyz position in space of both 
%  a lat/lon point, and a point on the tangent plane.  It's helpful to
%  define d(x,y)=perpendicular distance to Earth from tangent plane 
%  where  d(x,y)=R - sqrt(R.^2-x^2-y^2)

%  Also note, correctly, that y>0 at the north pole for delta lon=180 
%                         but y>0 at the south pole for delta lon=0;

if ~isempty(index)
    x(index)=R*coslat(index).*sind(lon(index));
    y(index)=-R*coslat(index).*sinlato(index).*coslon(index)+R.*coslato(index).*sinlat(index);
end

%d=2*R.*asin(frac(1,sqrt(2)).*sqrt(1-sqrt(1-frac(x.^2+y.^2,R.^2)))); 
if nargout==3
    if ~isempty(index)
        d=nan*lon;
        d(index)=2*R.*asin(frac(1,sqrt(2)).*sqrt(1-sqrt(1-frac(x(index).^2+y(index).^2,R.^2)))); 
    end
end


function[]=latlon2xy_test
latlon2xy_sphere_test1
latlon2xy_sphere_test2
latlon2xy_sphere_test3
latlon2xy_sphere_test4

function[]=latlon2xy_sphere_test1
rng(1);
N=1000;
tol=1e-3;

lon=2*pi*rand(N,1)-pi;
lat=pi*rand(N,1)-pi/2;
[lat,lon]=jrad2deg(lat,lon);

lat=lat/1000;
lon=lon/1000;

tic;[x,y]=latlon2xy(lat,lon,0,0,'sphere');etime1=toc;
tic;[x2,y2]=latlon2xy(lat,lon,0,0,'small');etime2=toc;

b=aresame(x,x2,tol) && aresame(y,y2,tol);
reporttest('LATLON2XY Cartesian and spherical algorithms match for small LAT and LON about zero',b);

function[]=latlon2xy_sphere_test2
rng(2);
N=100;
tol1=1e-1;
tol2=1e-1;

lon=2*pi*rand(N,1)-pi;
lat=pi*rand(N,1)-pi/2;
[lat,lon]=jrad2deg(lat,lon);

lat=lat/1000;
lon=lon/1000;

lono=2*pi*rand(N,1)-pi;
lato=pi*rand(N,1)-pi/2;
[lato,lono]=jrad2deg(lato,lono);

lat=lat+lato;
lon=lon+lono;

clear x y lat2 lon2

for i=1:length(lato)
    [x(i,1),y(i,1)]=latlon2xy(lat(i),lon(i),lato(i),lono(i),'sphere');
    [lat2(i,1),lon2(i,1)]=xy2latlon(x(i),y(i),lato(i),lono(i),'sphere');
end

[x2,y2]=latlon2xy(lat,lon,lato,lono);

b=aresame(x,x2,tol1) && aresame(y,y2,tol1);
reporttest('LATLON2XY Cartesian and spherical algorithms match for small LAT and LON perturbations',b);


%figure,plot([lon lon2 lon-lon2])
%This is not the best test as some random numbers can lead to a irrelevant 
%360 degree offset... so that's why I set rng
b=aresame(lat,lat2,tol2) && aresame(lon,lon2,tol2);
reporttest('XY2LATLON Cartesian and spherical algorithms match for small LAT and LON perturbations',b);



function[]=xy2latlon_test
rng(1);
latc=44;
lonc=0;
N=100;

x=randn(N,1)*5;
y=randn(N,1)*5;
[lat,lon]=xy2latlon(x,y,latc,lonc,'small');
[x2,y2]=latlon2xy(lat,lon,latc,lonc,'small');

tol=1e-10;
bool=aresame(x2,x,tol).*aresame(y2,y,tol);
reporttest('XY2LATLON / LATLON2XY conversion', bool)

latc=-44;
lonc=180;
N=100;

x=randn(N,1)*5;
y=randn(N,1)*5;
[lat,lon]=xy2latlon(x,y,latc,lonc,'small');
[x2,y2]=latlon2xy(lat,lon,latc,lonc,'small');

tol=1e-10;
bool(2)=aresame(x2,x,tol).*aresame(y2,y,tol);
reporttest('XY2LATLON / LATLON2XY conversion at 180', bool)

N=100;
tol=1e-6;

R=radearth;
x=frac(R,sqrt(2)).*(2.*rand(N,1)-1);
y=frac(R,sqrt(2)).*(2.*rand(N,1)-1);

lato=2*pi*rand(N,1)-pi;
lono=pi*rand(N,1)-pi/2;
[lato,lono]=jrad2deg(lato,lono);

clear lat lon x2 y2

for i=1:length(lato)
    [lat(i,1),lon(i,1)]=xy2latlon(x(i),y(i),lato(i),lono(i),'sphere');
    [x2(i,1),y2(i,1)]=latlon2xy(lat(i),lon(i),lato(i),lono(i),'sphere');
    [x3(i,1),y3(i,1)]=latlon2xy(lat(i),lon(i),lato(i),lono(i));
end


b=aresame(x,x2,tol) && aresame(y,y2,tol);
reporttest('LATLON2XY spherical algorithm inverts XY2LATLON spherical algorithm',b);

function[]=latlon2xy_sphere_test3
rng(1);
lon=(-180:5:175);
lat=(-90:5:90);
[latg,long]=meshgrid(lat,lon);

N=100;
lat1=180*rand(1,1,N)-90;    
lon1=360*rand(1,1,N);   

latg=vrep(latg,N,3);
long=vrep(long,N,3);

lat1=vrep(vrep(lat1,size(latg,1),1),size(latg,2),2);
lon1=vrep(vrep(lon1,size(latg,1),1),size(latg,2),2);

[x,y,d]=latlon2xy(lat1,lon1,latg,long);
d2=spheredist(lat1,lon1,latg,long);

b1=1;
b2=1;
index=find(d2<radearth*(pi/2));
if ~isempty(index)
    b1=aresame(d2(index),d(index),1e-4);
end
index=find(d2>radearth*(pi/2));
if ~isempty(index)
    b2=allall(isnan(d(index)));
end
        
reporttest('LATLON2XY matches SPHEREDIST for points in same hemisphere',b1);
reporttest('LATLON2XY gives NANs for points in other hemisphere',b2);


function[]=latlon2xy_sphere_test4
rng(1);
N=100;
tol=1e-3;

lon=2*pi*rand(N,1)-pi;
lat=pi*rand(N,1)-pi/2;
[lat,lon]=jrad2deg(lat,lon);

lat=lat/1000;
lon=lon/1000;

[x,y]=latlon2xy(lat,lon,0,0,'sphere');

latc{1}=lat;lonc{1}=lon;lonoc{1}=0;latoc{1}=0;
latc{2}=lat;lonc{2}=lon;lonoc{2}=0;latoc{2}=0;

[xc,yc]=latlon2xy(latc,lonc,latoc,lonoc,'sphere');

b=aresame(xc{1},x,tol) && aresame(yc{1},y,tol)&&aresame(xc{2},x,tol) && aresame(yc{2},y,tol);
reporttest('LATLON2XY cell input / output',b);

