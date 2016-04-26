function[cv,ca]=latlon2uv(varargin)
%LATLON2UV  Converts latitude and longitude to horizontal velocity.
%
%   CV=LATLON2UV(NUM,LAT,LON) where NUM is the data in DATENUM format
%   and LAT and LON are the latitude and longitude in degrees, outputs the
%   complex-valued velocity CV=U+SQRT(-1)*V in cm/s, computed using the 
%   first central difference.  
%
%   NUM is a column vector or a matrix of the same size as LAT and LON.  
%   LAT and LON are matices having SIZE(NUM,1) rows.  CV is the same size
%   as LAT and LON. NANs in LAT or LON become NAN+SQRT(-1)*NAN.
%  
%   LATLON2UV computes the velocity components from the distance travelled
%   across the surface of the sphere and the heading, taking account of 
%   the sphericity of the earth.
%
%   The first and last points must be treated differently, as the central 
%   difference is not defined there.  At the first and last point we use,
%   respectively, the first forward and first backward difference.  
%
%   The radius of the earth is given by the function RADEARTH.
%
%   [CV,CA]=LATLON2UV(...,'acceleration') also returns CA, a measure of the
%   complex-valued acceleration based on the difference between the first
%   forward difference velocity and the first backward difference velocity.
%   The units of acceleration are centimeters per inverse second squared.
%   ___________________________________________________________________
%
%   Cell array input / output
%
%   LATLON2UV returns cell array output given cell array input.  
%
%   That is, if NUM, LAT, and LON, are all cell arrays of length K, 
%   containing K different numerical arrays, then the output will also be 
%   cell arrays of length K.  
%   ___________________________________________________________________
%
%   Difference algorithm
%
%   By default, LATLON2UV uses a first central difference algorithm which 
%   is defined to be the average of a forward and a backward difference.
%
%   CV=LATLON2UV(...,'forward') instead uses a first forward difference. 
%   Velocity is computed based on the great circle distance and bearing 
%   from each point to the next point.  
%
%   CV=LATLON2UV(...,'backward') instead uses a first backward difference.  
%   Velocity is computed based on the great circle distance and bearing 
%   from each point to the previous point.  This is equivalent to the first 
%   forward difference computed backward in time.
%
%   The first forward and first backward difference can be thought of as 
%   giving the Cartesian velocities of a departing and arriving particle,
%   respectively.
%
%   Note that acceleration is not output when the 'forward' or 'backward'
%   algorithms are specified, as acceleration is computed as a part of the 
%   first central difference operation.
%   ___________________________________________________________________
%
%   Parallelization
%
%   LATLON2UV(...,'parallel') with cell array input parallelizes the 
%   computation using a PARFOR loop.  This requires that Matlab's Parallel
%   Computing Toolbox be installed, and is useful for very large datasets.
%   __________________________________________________________________
%
%   See also XY2LATLON, LATLON2XY, UV2LATLON.
%
%   Usage:  cv=latlon2uv(num,lat,lon);
%           [cv,ca]=latlon2uv(num,lat,lon);
%
%   'latlon2uv --t' runs a test.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details        
  

if strcmpi(varargin, '--t')
    latlon2uv_test,return
end

str='centered';
cores='serial';

for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'par')||strcmpi(varargin{end}(1:3),'ser')
            cores=varargin{end};
        else
            str=varargin{end};
        end
    end
end

num=varargin{1};
lat=varargin{2};
lon=varargin{3};
varargin=varargin(4:end);

if strcmpi(cores(1:3),'par')
    if exist('parpool')~=2
        disp('Sorry, parallel algorithm requires the Parallel Computing Toolbox.')
        disp('Defaulting to the standard algorithm.')
        cores='serial';
    end
end

if nargout==2
    str='acc';  %Output acceleration
end

if ~iscell(num)
    [cv,ca]=latlon2uv_celloop(num,lat,lon,str);
else
    if strcmpi(cores(1:3),'ser')
        for i=1:length(num)
            [cv{i,1},ca{i,1}]=latlon2uv_celloop(num{i},lat{i},lon{i},str);
        end
    elseif strcmpi(cores(1:3),'par')
        parfor i=1:length(num)
            [cv{i,1},ca{i,1}]=latlon2uv_celloop(num{i},lat{i},lon{i},str);
        end
    end
end

if ~strcmpi(str(1:3),'acc')
    clear ca
end


function[cv,ca]=latlon2uv_celloop(num,lat,lon,str)

ca=[];
     
if length(lat)<=2
    cv=nan+0*lat;
    ca=nan+0*lat;
else
    if size(num,1)==1
        if size(num')==size(lat)
            num=num';
        end
    end
    if size(num,2)==1
        num=vrep(num,size(lat,2),2);
    end
    if strcmpi(str(1:3),'cen')||strcmpi(str(1:3),'acc')
        [u1,v1]=latlon2uv_celloop_one(num,lat,lon,'forward');
        [u2,v2]=latlon2uv_celloop_one(num,lat,lon,'backward');
        
        %Correction for first and last points
        u2(1,:)=u1(1,:);
        v2(1,:)=v1(1,:);
        u1(end,:)=u2(end,:);
        v1(end,:)=v2(end,:);
   
        cv=frac(1,2)*(u1+u2)+1i*frac(1,2)*(v1+v2);
   
        if strfind(str,'acc')
            dt=vdiff(num*24*3600,1,1);
            ca=((u1-u2)+1i*(v1-v2))./dt;
            %Correction for first and last points
            ca([1 end],:)=ca([2 end-1],:);
        end
    else
        [u,v]=latlon2uv_celloop_one(num,lat,lon,str);
        cv=u+1i*v;
    end
end


function[u,v]=latlon2uv_celloop_one(num,lat,lon,str)


if strfind(str,'old')
    [dr,dt,gamma]=latlon2uv_old(num,lat,lon,str);
elseif strfind(str,'for')
    [dr,dt,gamma]=latlon2uv_forward(num,lat,lon,str);
elseif strfind(str,'bac')
    [dr,dt,gamma]=latlon2uv_forward(flip(num,1),flip(lat,1),flip(lon,1),str);
    dr=flip(dr,1);
    dt=flip(dt,1);
    gamma=flip(gamma,1);
end


%Convert to centimeters
c=100*1000;
u=c.*dr./dt.*cos(gamma);
v=c.*dr./dt.*sin(gamma);
index=find(isnan(u)|isnan(v)); 
if ~isempty(index)
    u(index)=nan;
    v(index)=nan;
end 

function[dr,dt,gamma]=latlon2uv_old(num,lat,lon,str)
%Former version of a first central difference algorithm

%[phi,th]=jdeg2rad(lat,lon);
[x,y,z]=latlon2xyz(lat,lon);


dt=vdiff(num*24*3600,1);
[dx,dy,dz]=vdiff(x,y,z,1);
xm=vshift(x,1,1)/2+vshift(x,-1,1)/2;
ym=vshift(y,1,1)/2+vshift(y,-1,1)/2;
zm=vshift(z,1,1)/2+vshift(z,-1,1)/2;

%Correction for first and last point
xm([1 end],:)=x([1 end],:);
ym([1 end],:)=y([1 end],:);
zm([1 end],:)=z([1 end],:);


[latm,lonm]=xyz2latlon(xm,ym,zm);
[phim,thm]=jdeg2rad(latm,lonm);

%Now perform a local rotation
dx2=dx.*cos(-thm)-dy.*sin(-thm);
dy2=dx.*sin(-thm)+dy.*cos(-thm);
dz2=dz;

%dx3=dx2.*cos(-phim)-dz2.*sin(-phim);    %Correct but not needed
dy3=dy2; 
dz3=dx2.*sin(-phim)+dz2.*cos(-phim);

[lat1,lon1]=vshift(lat,lon,1,1);
[lat2,lon2]=vshift(lat,lon,-1,1);
dr=spheredist(lat1,lon1,lat2,lon2)/2;  

%Treat first and last points differently
dr(1,:)=spheredist(lat(1,:),lon(1,:),lat(2,:),lon(2,:));  
dr(end,:)=spheredist(lat(end-1,:),lon(end-1,:),lat(end,:),lon(end,:));  

gamma=imlog(dy3+sqrt(-1)*dz3);


function[dr,dt,gamma]=latlon2uv_forward(num,lat,lon,str)

%lon=degunwrap(lon);
%First forward difference
dt=diff(num*24*3600,1,1);

[lat2,lon2]=vshift(lat,lon,1,1);
dr=spheredist(lat,lon,lat2,lon2); %figure,plot(dr)

xx=sind(lon2-lon).*cosd(lat2);
yy=cosd(lat).*sind(lat2)-sind(lat).*cosd(lat2).*cosd(lon2-lon);

gamma=atan2(yy,xx);
%figure,uvplot(diff(rot(gamma)))
%gamma2=imag(log(xx+1i*yy));
%aresame(gamma,gamma2,1e-10)

%From http://www.movable-type.co.uk/scripts/latlong.html
%ATAN2(COS(lat1)*SIN(lat2)-SIN(lat1)*COS(lat2)*COS(lon2-lon1),SIN(lon2-lon1)*COS(lat2)) 
%But note, that source actually reverses the inputs to the inverse tangent function 

dt(end+1,:)=dt(end,:);
dr(end,:)=dr(end-1,:);
gamma(end,:)=gamma(end-1,:);

function[]=latlon2uv_test
latlon2uv_dlat_test
latlon2uv_dlon_test
latlon2uv_small_test
latlon2uv_displacement_test
latlon2uv_npg2006_test
latlon2uv_npg2006_nonuniform_test
latlon2uv_testformer

function[]=latlon2uv_dlat_test


N=100;
tol=1e-3;

lono=2*pi*rand(1,N)-pi;
lato=pi*rand(1,N)-pi/2;
[lato,lono]=jrad2deg(lato,lono);

dlat=(rand(1,N)-1/2)/10;

lat=[lato-dlat;lato;lato+dlat];
lon=[lono;lono;lono];

num=[-1+0*lato;0*lato;1+0*lato];
index=find(max(lat)<90&min(lat)>-90);
vindex(num,lat,lon,dlat,index,2);

cv1=0*dlat+1i*100*1000*2*pi*radearth/360.*dlat/(3600*24);
cv2=latlon2uv(num,lat,lon);
b=aresame(cv2(2,:),cv1,tol);
reporttest('LATLON2UV small delta latitude',b);



function[]=latlon2uv_dlon_test


N=100;
tol=1e-3;

lono=2*pi*rand(1,N)-pi;
lato=pi*rand(1,N)-pi/2;
[lato,lono]=jrad2deg(lato,lono);

dlon=(rand(1,N)-1/2)/10;

lat=[lato;lato;lato];
lon=[lono-dlon;lono;lono+dlon];

num=[-1+0*lato;0*lato;1+0*lato];

cv1=100*1000*2*pi*radearth/360.*dlon/(3600*24).*cos(jdeg2rad(lato))+1i*0*lato;
cv2=latlon2uv(num,lat,lon);
b=aresame(cv2(2,:),cv1,tol);
reporttest('LATLON2UV small delta longitude',b);


function[]=latlon2uv_small_test


N=100;
tol=5e-3;

lono=2*pi*rand(1,N)-pi;
lato=pi*rand(1,N)-pi/2;
[lato,lono]=jrad2deg(lato,lono);

dlon=(rand(1,N)-1/2)/10;
dlat=(rand(1,N)-1/2)/10;

lat=[lato-dlat;lato;lato+dlat];
lon=[lono-dlon;lono;lono+dlon];

num=[-1+0*lato;0*lato;1+0*lato];
index=find(max(lat)<90&min(lat)>-90);
vindex(num,lat,lon,dlat,dlon,lato,index,2);

u1=100*1000*2*pi*radearth/360.*dlon/(3600*24).*cos(jdeg2rad(lato));
v1=100*1000*2*pi*radearth/360.*dlat/(3600*24);
cv1=u1+1i*v1;
cv=latlon2uv(num,lat,lon);

b=aresame(cv(2,:),cv1,tol);
reporttest('LATLON2UV small displacements',b);

function[]=latlon2uv_displacement_test


N=100;
tol=2e-3;

lon=2*pi*rand(N,1)-pi;
lat=pi*rand(N,1)-pi/2;

num=(1:N)';
cv=latlon2uv(num,lat,lon);

dr=abs(cv).*(3600.*24)./(100*1000);

[lat1,lon1]=vshift(lat,lon,1,1);
[lat2,lon2]=vshift(lat,lon,-1,1);

dr1=spheredist(lat1,lon1,lat2,lon2)/2;

b=max(frac(abs(dr(2:end-1)-dr1(2:end-1)),dr(2:end-1)))<tol;
reporttest('LATLON2UV displacement matches SPHEREDIST',b);


function[]=latlon2uv_npg2006_test
load npg2006
use npg2006
cv1=latlon2uv(num,lat,lon);
cv2=vdiff(cx,1).*100.*1000./(3600.*24.*dt);
tol=1/2;
b=aresame(cv1(2:end-1,:),cv2(2:end-1,:),tol);
reporttest('LATLON2UV matches VDIFF of Cartesian position for NPG2006 data',b);


function[]=latlon2uv_npg2006_nonuniform_test
load npg2006
use npg2006
cv1=latlon2uv([num 2*num],[lat lat],[lon lon]);
cv2=vdiff(cx,1).*100.*1000./(3600.*24.*dt);

tol=1/2;
b=aresame(cv1(2:end-1,1),cv2(2:end-1,1),tol)&&aresame(cv1(2:end-1,2),cv2(2:end-1,1)./2,tol);
reporttest('LATLON2UV matches NPG2006 data for DT different by columns, central diff',b);


cv1=latlon2uv([num 2*num],[lat lat],[lon lon],'forward');
cv2=diff(cx,1,1).*100.*1000./(3600.*24.*dt);

tol=1/2;
b=aresame(cv1(1:end-1,1),cv2(:,1),tol)&&aresame(cv1(1:end-1,2),cv2(:,1)./2,tol);
reporttest('LATLON2UV matches NPG2006 data for DT different by columns, forward diff',b);

cv1=latlon2uv(num,lat,lon,'forward');
cv2=latlon2uv(num,lat,lon,'backward');

tol=1/2;
b=aresame(cv1(1:end-1,:),cv2(2:end,:),tol);
reporttest('LATLON2UV forward and shifted backward difference agree approximately',b);

cv3=latlon2uv(num,lat,lon,'central');

b=aresame(frac(cv1(2:end-1,:)+cv2(2:end-1,:),2),cv3(2:end-1,:));
reporttest('LATLON2UV central difference is exactly average of forward and backward',b);

cv4=latlon2uv(num,lat,lon,'old');
b=aresame(cv3(2:end-1,:),cv4(2:end-1,:),0.02);
reporttest('LATLON2UV central difference matches old algorithm for NPG2006 data',b);


% cx=10*rot([0:0.001:4*pi]');
% [lat,lon]=xy2latlon(cx,65,0);
% num=[1:length(cx)]'/24/100;
% cv=latlon2uv(num,lat,lon);
% ca=latlon2uv(num,lat,lon,'acc');
% ca2=vdiff(cv,1)./(1/24/100*24*3600);
% 
% uvplot(num,ca), hold on, uvplot(num,ca2)
%
%But there is a problem with eastward acceleration for very small amplitude
%circles.  Keep an eye out for noise coming from latlon2uv.  

function[]=latlon2uv_testformer

lon1=[1:360]';
lat1=60 +0*lon1;

lon2=180+0*lon1;
lat2=0+0*lon1;

lon3=[360:-1:1]';
lat3=-60 +0*lon1;


lon=[lon1';lon2';lon3'];lon=lon(:);
lat=[lat1';lat2';lat3'];lat=lat(:);

num=[1:length(lat)]';

clear cv
cv(:,1)=latlon2uv(num,lat,lon,'forward');
cv(:,2)=latlon2uv(num,lat,lon,'backward');
cv(:,3)=latlon2uv(num,lat,lon,'central');
cv(:,4)=latlon2uv(num,lat,lon,'old');

cv=cv(2:3:end,:);
lat=lat(1:3:end,:);
lon=lon(1:3:end,:);


ii=[91:269]';
b=aresame(1+0*ii,cv(ii,4)./cv(ii,3),0.0001);

reporttest('LATLON2UV former algorithm and new algorthm agree to 1e-4 when Delta Lon<90',b);








