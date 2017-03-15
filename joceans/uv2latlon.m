function[lat,lon]=uv2latlon(varargin)
%UV2LATLON  Integrates horizontal velocity to give latitude and longitude.
%
%   [LAT,LON]=UV2LATLON(NUM,U,V,LATO,LONO) where NUM is the data in DATENUM
%   format, integates velocities U and V from intial location (LATO,LONO)
%   to form a trajectory on the earth described by LAT and LON.  
%
%   LATO and LONO are an initial latitude and longitude in degrees,
%   and U and V are eastward and northward velocity components in cm/s. 
%
%   NUM, U, and V may be column vectors, in which case LATO and LONO are 
%   both scalars.  Alternatively U and V may be matrices with time oriented
%   in rows.  NUM is then either an array of length SIZE(U,1), or a matrix
%   of the same size as U and V, while LATO and LONO are either scalars or 
%   arrays of length SIZE(U,2). 
%
%   [LAT,LON]=UV2LATLON(NUM,CV,LATO,LONO), where CV is the complex-valued
%   velocity CV=U+SQRT(-1)*V, also works.
%
%   UV2LATLON is inverted by LATLON2UV.  
%   ___________________________________________________________________
%
%   Cell array input / output
%
%   UV2LATLON returns cell array output given cell array input.  
%
%   That is, if NUM, U, and V, are all cell arrays of length K, containing
%   K different numerical arrays, then the output will also be cell arrays
%   of length K.  
%
%   In this case LATO and LONO are either scalars, or arrays of length K.
%   ___________________________________________________________________
%
%   Algorithm
%
%   UV2LATLON works by converting the velocity at the current point into 
%   a vector in 3D Cartesian coordinates, using SPHERE2UVW, finding the new
%   location in 3D space, and then converting this back into latitude 
%   and longitude using XYZ2LATLON, and then interating for the next point.
%
%   By default, UV2LATLON uses a forward integration from an initial point. 
%   This inverts LATLON2UV(...,'forward') to a high degree of precision for
%   typical drifter and float velocity values and sampling intervals.
%
%   UV2LATLON(NUM,CV,LATF,LONF,'backward') instead uses a *backward* 
%   integration, and LATF and LONF are now the *final* points rather than 
%   the initial points.  This inverts LATLON2UV(...,'backward'). 
%   ___________________________________________________________________
% 
%   See also XY2LATLON, LATLON2XY, LATLON2UV.
%
%   'uv2latlon --t' runs a test.
%
%   Usage: [lat,lon]=uv2latlon(num,u,v,lato,lono);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014 J.M. Lilly --- type 'help jlab_license' for details

if strcmpi(varargin{1}, '--t')
    uv2latlon_test,return
end

if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
else
    str='forward';
end

na=length(varargin);
num=varargin{1};
if na==4
    if iscell(varargin{2})
        u=cellreal(varargin{2});
        v=cellimag(varargin{2});
    else
        u=real(varargin{2});
        v=imag(varargin{2});
    end
    varargin=varargin(3:end);
else
    u=varargin{2};
    v=varargin{3};
    varargin=varargin(4:end);
end
lato=varargin{1};
lono=varargin{2};

if ~iscell(num)
    if ~aresame(size(num),size(u))
        num=vrep(num(:),size(u,2),2);
    end
    if length(lato)==1
        lato=lato+zeros(size(u(1,:)));
    end
    if length(lono)==1
        lono=lono+zeros(size(u(1,:)));
    end
    %vsize(num,lato,lono,u,v)
    [lat,lon]=uv2latlon_one(num,lato,lono,u,v,str);
else
    if length(lato)==1
        lato=lato+zeros(size(u));
    end
    if length(lono)==1
        lono=lono+zeros(size(u));
    end
    for i=1:length(num)
        [lat{i,1},lon{i,1}]=uv2latlon_one(num{i},lato(i),lono(i),u{i},v{i},str);
    end
end

function[lat,lon]=uv2latlon_one(num,lato,lono,u,v,str)

switch str(1:3)
    case 'for'
        [lat,lon]=uv2latlon_integrator(num,lato,lono,u,v);
    case 'bac'
        [lat,lon]=uv2latlon_integrator(num,lato,lono,-flipud(u),-flipud(v));
        lat=flipud(lat);
        lon=flipud(lon);
%     case 'cen'
%         [lat1,lon1]=uv2latlon_integrator(num,lato(1),lono(1),u,v);
%         [lat2,lon2]=uv2latlon_integrator(num,lato(2),lono(2),-flipud(u),-flipud(v));
%         lat2=flipud(lat2);
%         lon2=flipud(lon2);
%         [x1,y1,z1]=latlon2xyz(lat1,lon1);
%         [x2,y2,z2]=latlon2xyz(lat2,lon2);
%         [lat,lon]=xyz2latlon((x1+x2)/2,(y1+y2)/2,(z1+z2)/2);
end
        
function[lat,lon]=uv2latlon_integrator(num,lato,lono,u,v)

dt=num(2)-num(1);
[xo,yo,zo]=latlon2xyz(lato,lono);
[lat,lon]=vzeros(size(u));
lat(1,:)=lato;
lon(1,:)=lono;

u=u(1:end-1,:);
v=v(1:end-1,:);

% if strcmpi(str(1:3),'for')
%     u=u(1:end-1,:);
%     u=u(1:end-1,:);
% elseif strcmpi(str(1:3),'for')
%     u=u(1:end-1,:);
%     u=u(1:end-1,:);
% 
% u=1/2*(u+circshift(u,-1));
% u=u(1:end-1,:);
% 
% v=1/2*(v+circshift(v,-1));
% v=v(1:end-1,:);

u=u*(dt*24*3600)/100/1000;%Now it's km
v=v*(dt*24*3600)/100/1000;%Now it's km

%vsize(lat,lon,u,v)
for i=1:size(u,1)
    [dx,dy,dz]=sphere2uvw(lat(i,:),lon(i,:),0,u(i,:),v(i,:));
    xo=xo+dx;
    yo=yo+dy;
    zo=zo+dz;
    [lat(i+1,:),lon(i+1,:)]=xyz2latlon(xo,yo,zo);
end

 
function[]=uv2latlon_test
load npg2006

use npg2006
vindex(num,lat,lon,1:max(find(~isnan(lat))),1);
cv=latlon2uv(num,lat,lon,'forward');
[lat2,lon2]=uv2latlon(num,cv,lat(1),lon(1),'forward');

bool1=allall(abs(frac(lon-lon2,lon))<5e-6);
bool2=allall(abs(frac(lat-lat2,lat))<5e-6);
reporttest('UV2LATLON inverts LATLON2UV to within 5e-6 for NPG2006 data, forward algorithm',bool1&&bool2)

[lat2,lon2]=uv2latlon([num num],[cv cv],[lat(1) lat(1)],[lon(1) lon(1)],'forward');
bool1=allall(abs(frac([lon lon]-lon2,[lon lon]))<5e-6);
bool2=allall(abs(frac([lat lat]-lat2,[lat lat]))<5e-6);
reporttest('UV2LATLON inverts LATLON2UV to within 5e-6 for NPG2006 data, forward algorithm, matrix input',bool1&&bool2)


use npg2006
vindex(num,lat,lon,1:max(find(~isnan(lat))),1);
cv=latlon2uv(num,lat,lon,'backward');
[lat3,lon3]=uv2latlon(num,cv,lat(end),lon(end),'backward');

bool1=allall(abs(frac(lon-lon3,lon))<5e-6);
bool2=allall(abs(frac(lat-lat3,lat))<5e-6);
reporttest('UV2LATLON inverts LATLON2UV to within 5e-6 for NPG2006 data, backward algorithm',bool1&&bool2)

% use npg2006
% vindex(num,lat,lon,1:max(find(~isnan(lat))),1);
% cv=latlon2uv(num,lat,lon);
% [lat4,lon4]=uv2latlon(num,cv,[lat(1) lat(end)],[lon(1) lon(end)],'central');
% 
% bool1=allall(abs(frac(lon-lon4,lon))<5e-6);
% bool2=allall(abs(frac(lat-lat4,lat))<5e-6);
% reporttest('UV2LATLON inverts LATLON2UV to within 5e-6 for NPG2006 data, central algorithm',bool1&&bool2)
