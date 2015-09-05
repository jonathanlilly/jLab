function[Ro]=ellrossby(lat,lambda,omega)
%ELLROSSBY  Ellipse Rossby number, for oceanographic applications.
%
%   RO=ELLROSSBY(LAT,LAMBDA,OMEGA) returns the ellipse Rossby number RO 
%   for an ellipse with linearity LAMBDA and joint instantaneous frequency
%   OMEGA, in radians per day, at latitude LAT.
%  
%   The magnitude of the ellipse Rossby number is |RO|=|OMEGA/F| where F
%   is the signed Coriolis frequency at latitude LAT.  
%
%   The sign of RO is SIGN(RO)=-SIGN(LAMBDA/F), so that RO is positive for
%   cyclones, and negative for anticyclones, in both hemispheres. 
%
%   For circular eddies, this definition of Rossby number is equivalent to
%   RO=-V/(R*F), where R is the radial distance to the eddy center and V is
%   the signed instantaneous azimuthal velocity. 
%
%   The input arguments may also be cell arrays of numeric arrays, all 
%   having the same size.  RO will then be a similarly sized cell array.
%   Alternatively, LAMBDA and OMEGA may be cell arrays and LAT a constant.
%
%   'ellrossby --t' runs a test.
%
%   Usage: ro=ellrossby(lat,lambda,omega);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2014 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(lat, '--t')
    ellrossby_test,return
end

if ~isempty(omega)
    if ~iscell(omega)
        Ro=ellrossby_one(lat,lambda,omega);
        Ro(isinf(lat.*omega))=inf;
    else
        for i=1:length(omega)
            if iscell(lat)
                lati=lat{i};
            else
                lati=lat;
            end
            Ro{i,1}=ellrossby_one(lati,lambda{i},omega{i});
            Ro{i,1}(isinf(lati.*omega{i}))=inf;
        end
    end
else
    Ro=omega;
end

function[Ro]=ellrossby_one(lat,lambda,omega)

fcor=(corfreq(lat))*24;  %Coriolis frequency in radians per day at center 
Ro=sign(lambda).*frac(omega,fcor); %Rossby number under solid-body assumption
    
    
function[]=ellrossby_test
 

%/*************************************************
load ebasnfloats
use ebasnfloats
num=num{33};lat=lat{33};lon=lon{33};p=p{33};t=t{33};
vindex(num,lat,lon,1:549,1);

cx=latlon2xy(lat,lon,vmean(lat(:),1),vmean(lon(:),1));
cv=latlon2uv(num,lat,lon);

ga=3;be=3;

dt=(num(2)-num(1));
mlat=vmean(lat(:),1);

fmax=abs(corfreq(mlat))*frac(dt*24,2);  %One cycle per 2 inertial periods = 1 cycle per 2.6 day
fmin=abs(corfreq(mlat))*frac(dt*24,40); %One cycle per 40 inertial periods = 1 cycle per 53 days
fs=morsespace(ga,be,{0.2,fmax},fmin,8);

%Compute wavelet transforms using generalized Morse wavelets
[wx,wy]=wavetrans(real(cx),imag(cx),{ga,be,fs,'bandpass'},'mirror');

%Form ridges of component time series
[ir,jr,wxr,wyr,fxr,fyr]=ridgewalk(dt,wx,wy,fs,{3,0});   

%Map into time series locations
[wrx,frx]=ridgemap(length(cx),wxr,fxr,ir);
[wry,fry]=ridgemap(length(cx),wyr,fyr,ir);

[kappa,lambda,theta,phi]=ellparams(wrx,wry);
om=vmean([frx fry],2,squared([wrx wry]));
ro=ellrossby(lat,lambda,om);

rm=ellrad(kappa,lambda);
vm=ellvel(24*3600,kappa,lambda,theta,phi,1e5);
om2=vm./rm*frac(24*3600,100*1000);
ro2=om2./abs(corfreq(lat))/24;

reporttest('ELLROSSBY matches V/(Rf) form to within 2% for EBASN anticyclone', sqrt(vmean(squared(ro-ro2),1))./abs(vmean(ro,1))<0.02)


