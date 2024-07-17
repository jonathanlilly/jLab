function[wxh,wyh]=spheretrans(varargin)
%SPHERETRANS  Wavelet transform for oscillations on the surface of a sphere.
%
%   [WX,WY]=SPHERETRANS(LAT,LON,PSI) implements a version of the wavelet
%   transform WAVETRANS appropriate for analyzing oscillations in latitude
%   LAT and longitude LON on the surface of a sphere using wavelet PSI.
%
%   By default, SPHERETRANS uses a sphere with the Earth's radius as given
%   by RADEARTH.  SPHERETRANS(R,...) insteads uses radius R, in kilometers.
%
%   SPHERETRANS works by first converting latitude and longitude to three-
%   dimensional position in space, computing the wavelet transform in 3D,
%   and then projecting this back onto a tangent plane centered on the 
%   time-varying center of the oscillation in each wavelet band. 
%   
%   See WAVETRANS for the format of PSI.  It can be either a matrix or a
%   cell array of parameters. 
%
%   SPHERETRANS accepts any of the trailing strings accepted by WAVETRANS, 
%   such as SPHERETRANS(...,'parallel') to specify using a PARFOR loop.
%
%   'spheretrans --t' runs a test.
%
%   Usage: [wx,wy]=spheretrans(lat,lon,psi);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2016--2023 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    %spheretrans_test;hor2uvw_test;uvw2hor_test,return
    hor2uvw_test;uvw2hor_test,return
end
 
if ~iscell(varargin{1})&&length(varargin{1})==1
    R=varargin{1};
    varargin=varargin(2:end);
end
lat=varargin{1};
lon=varargin{2};
psi=varargin{3};

%String arguments to wavetrans
str1='periodic';
str2='nodetrend';
str3='series';
for i=1:3
    varend=varargin{end};
    if ischar(varend)
        if strcmpi(varend(1:3),'ser')||strcmpi(varend(1:3),'par')
            str3=varend;
        elseif strcmpi(varend(1:3),'det')||strcmpi(varend(1:3),'nod')
            str2=varend;
        else
            str1=varend;
        end
        varargin=varargin(1:end-1);
    end
end

if ~iscell(lon)
    lon=frac(360,2*pi)*unwrap(frac(2*pi,360)*lon);
    [x,y,z]=latlon2xyz(lat,lon);%Position as a 3-vector with origin at center of earth
else
    for i=1:length(lon)
        lon{i}=frac(360,2*pi)*unwrap(frac(2*pi,360)*lon{i});
        [x{i},y{i},z{i}]=latlon2xyz(lat{i},lon{i});%Position as a 3-vector with origin at center of earth
    end
end

[wx,wy,wz]=wavetrans(x,y,z,psi,str1,str2,str3);


if ~iscell(lon)
    [wxh,wyh]=spheretrans_xyz2xy(lat,lon,x,y,z,wx,wy,wz);
else
    for i=1:length(lon)
        [wxh{i,1},wyh{i,1}]=spheretrans_xyz2xy(lat{i},lon{i},x{i},y{i},z{i},wx{i},wy{i},wz{i});
    end
end

%I think that the real answer is that I should decompose velocity, 
%not displacement. 

function[wxh,wyh]=spheretrans_xyz2xy(lat,lon,x,y,z,wx,wy,wz)

%Project wavelet transform onto local horizontal
latmat=vrep(lat,size(wx,2),2);
lonmat=vrep(lon,size(wx,2),2);

%Subtract myself
xmat=vrep(x,size(wx,2),2);
ymat=vrep(y,size(wx,2),2);
zmat=vrep(z,size(wx,2),2);
[latnew,lonnew]=xyz2latlon(xmat-real(wx),ymat-real(wy),zmat-real(wz));
[wxh,wyh]=uvw2hor(latnew,lonnew,wx,wy,wz);
%[wxh,wyh]=uvw2hor(vrep(lat,size(wx,2),2),vrep(lon,size(wx,2),2),wx,wy,wz);


function[uh,vh]=uvw2hor(lat,lon,u,v,w)
%UVW2HOR  Projects a 3D Cartesian vector into a horizontal vector on a sphere.
%
%   [UH,VH]=UVW2HOR(LAT,LON,U,V,W) takes the 3D Cartesian vector with
%   components U, V, and W located at point (LAT,LON) and projects it to
%   find the local horizontal components of the vector UH and VH.  
%
%   LAT and LON are in degrees.
%
%   U, V, and W are in a reference frame with the X-axis at zero 
%   degrees longitude and the Z-axis at the North Pole.  
%
%   All input arguments should be arrays of the same size.   Note U, V, and
%   W may be complex-valued, in which case UH and VH will be also.
%
%   UVW2HOR inverts HOR2UVW, but the reverse is not true, since some
%   information is lost during the projection.
%
%   See JSPHERE for related functions.
%  
%   'uvw2hor --t' runs a test.
%
%   Usage: [uh,vh]=uvw2hor(lat,lon,u,v,w);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2014 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(lat, '--t')
    uvw2hor_test,return
end
 
[phi,theta]=jdeg2rad(lat,lon);

uh=v.*cos(theta)-u.*sin(theta);%
%vh=w./cos(phi);%
vh=w.*cos(phi)-u.*cos(theta).*sin(phi)-v.*sin(theta).*sin(phi);

%aresame(vh,w./cos(phi),1e-11)

function[]=uvw2hor_test

lon=(1e-10:2:360)-180;
%lon=0;
lat=(-90:2:90);
[lon,lat]=meshgrid(lon,lat);

rng(1)
uh=randn(size(lat));
vh=randn(size(lat));
%generate a random velocity in the tangent plane and convert to uvw
[u,v,w]=hor2uvw(lat,lon,uh,vh);
%convert uvw back to velocity in the horizontal plane 
[uh2,vh2]=uvw2hor(lat,lon,u,v,w);

%clf,plot(lat,vh,'b'),hold on,plot(lat,vh2,'r'),plot(lat,vh,'b')

tol=1e-6;
reporttest('UVW2HOR inverts HOR2UVW',aresame(uh,uh2,tol)&&aresame(vh,vh2,tol))


function[u,v,w]=hor2uvw(lat,lon,uh,vh)
%HOR2UVW  Converts a horizontal vector on a sphere into a 3D Cartesian vector.
%
%   [U,V,W]=HOR2UVW(LAT,LON,UH,VH) converts the vector with local 
%   horizontal components UH and VH at point (LAT,LON) on a sphere 
%   into a 3D Cartesian vector having components U, V, and W.
%
%   LAT and LON are in degrees.  
%
%   U, V, and W are in a reference frame with the X-axis at zero 
%   degrees longitude and the Z-axis at the North Pole.  
%
%   All input arguments should be arrays of the same size.
%
%   HOR2UVW is inverted by UVW2HOR.
%
%   See JSPHERE for related functions.
%
%   'hor2uvw --t' runs a test.
%
%   Usage: [u,v,w]=hor2uvw(lat,lon,uh,vh);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2014 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(lat, '--t')
    hor2uvw_test,return
end

[phi,theta]=jdeg2rad(lat,lon);

u=-uh.*sin(theta)-vh.*cos(theta).*sin(phi);%ok;ok
v= uh.*cos(theta)-vh.*sin(theta).*sin(phi);%ok;ok
w= vh.*cos(phi);%ok

function[]=hor2uvw_test
 
lat=[0  0  45         0        -90    45]';
lon=[0  90 0         90         0     0]';
uh= [1  0  1          1         1     0]';
vh= [0  1  0          0         0     1]';
u=  [0  0  0         -1         0     -sqrt(2)/2]'; 
v=  [1  0  1          0         1    0]';
w=  [0  1  0          0         0    sqrt(2)/2]';


[u2,v2,w2]=hor2uvw(lat,lon,uh,vh);

tol=1e-6;
reporttest('HOR2UVW example points',aresame(u,u2,tol)&&aresame(v,v2,tol)&&aresame(w,w2,tol))

lon=(0:2:358)-180;
lat=(-90:1:90);
[long,latg]=meshgrid(lon,lat);

uh=randn(size(latg));
vh=randn(size(latg));
[u,v,w]=hor2uvw(latg,long,uh,vh);
%[v1,v2,v3]=uvw2sphere(latg,long,u,v,w);

tol=1e-10;
%reporttest('HOR2UVW plus UVW2SPHERE for vanishing w',aresame(uh,v2,tol)&&aresame(vh,v3,tol)&&aresame(v1,0*v1,tol))

function[]=spheretrans_test

%Should test this with Jeffrey's inertial oscillation


lato=[2.5:5:87.5]';
[time,lat,lon,u,v]=vzeros(2881,length(lato));

for i=1:length(lato)
    [time(:,i),lon(:,i),lat(:,i),u(:,i),v(:,i),frequency]=ExactInertialSolution(lato(i),20,24*3600*30,15*60);
end


[time,lon,lat,u,v,frequency]=ExactInertialSolution(45,20,24*3600*30,15*60);

gamma=3;beta=3;
fs=morsespace(gamma,beta,{0.05,pi},1/1000,4);
%[wlat,wlon]=wavetrans(lat,lon,{gamma,beta,fs,'bandpass'},'mirror');
[wx,wy]=spheretrans(lat,0*lon,{gamma,beta,fs,'bandpass'},'mirror');






load ebasnfloats
use ebasnfloats
num=num{33};lat=lat{33};lon=lon{33};p=p{33};t=t{33};
vindex(num,lat,lon,1:547,1);
cx=latlon2xy(lat,lon,vmean(lat(:),1),vmean(lon(:),1));

%Decide on frequencies
fs=2*pi./(logspace(log10(1),log10(100),50)');

ga=3;
be=2;

[wx,wy]=spheretrans(lat,lon,{ga,be,fs,'bandpass'},'reverse');
[ir,jr,wr]=ridgewalk(1,wx,wy,fs,{1.5,0,'amplitude'});


%Compute wavelet transforms using generalized Morse wavelets
[wx2,wy2]=wavetrans(real(cx),imag(cx),{1,ga,be,fs,'bandpass'},'reverse');
[ir2,jr2,wr2]=ridgewalk(1,wx2,wy2,fs,{1.5,0,'amplitude'});


figure,
plot(ir2{1},imag(wr2{1})),hold on
plot(ir{1},imag(wr{1}))
plot(ir{1},real(wr{1}))
plot(ir2{1},real(wr2{1}))




%Form ridges of component time series



wr=wr(1:end-1,:);






[lathat,lonhat,latres,lonres]=xy2latlon(wr(:,1),wr(:,2),lat,lon);


%Should test this with Jeffrey's inertial oscillation
[time,lon,lat,u,v,frequency]=ExactInertialSolution(45,20,24*3600*30,15*60);

%[time,lon,lat,u,v,frequency]=ExactInertialSolution(85,20,24*3600*30,15*60);

gamma=3;beta=3;
fs=morsespace(gamma,beta,{0.05,pi},1/1000,4);
%[wlat,wlon]=wavetrans(lat,lon,{gamma,beta,fs,'bandpass'},'mirror');
[wx,wy]=spheretrans(lat,lon,{gamma,beta,fs,'bandpass'},'mirror');


%Still, 3rd possibility is to use velocity instead... 

region=[  -99   -80    18    31];
load highres
use highres
clear highres

[lat,lon,num,id]=floatregion(region,lat,lon,num,id);
cell2col(lat,lon,num,id);
col2mat(lat,lon,num,id);
vindex(lat,lon,num,id,39,2); 
ii=1:min(find(isnan(lat)))-1;
vindex(lat,lon,num,id,ii,1); 



gamma=3;beta=3;
fs=morsespace(gamma,beta,{0.05,pi},1/1000,4);
%[wlat,wlon]=wavetrans(lat,lon,{gamma,beta,fs,'bandpass'},'mirror');
[wx2,wy2]=spheretrans(lat+60,lon,{gamma,beta,fs,'bandpass'},'mirror','2D');
[wx3,wy3]=spheretrans(lat+60,lon,{gamma,beta,fs,'bandpass'},'mirror','3D');

[wlat,wlon]=wavetrans(lat,lon,{gamma,beta,fs,'bandpass'},'mirror');
 