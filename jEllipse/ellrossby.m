function[Ro]=ellrossby(lat,xi,omega)
%ELLROSSBY  Ellipse Rossby number, for oceanographic applications.
%
%   RO=ELLROSSBY(LAT,XI,OMEGA) returns the vorticity Rossby number RO for 
%   an ellipse with signed circularity XI and joint instantaneous
%   frequency OMEGA, given in radians per day, at latitude LAT.
%  
%   The vorticity Rossby number is defined as RO = (2/XI)*OMEGA/F where F
%   is the signed Coriolis frequency at latitude LAT.  
%
%   For details, see Lilly and Perez-Brunius (2021b).
%
%   For steady, non-precessing ellipses, RO will be positive for cyclonic
%   motion and negative for anticyclonic motion in both hemispheres. 
%
%   The input arguments may also be cell arrays of numeric arrays, all 
%   having the same size.  RO will then be a similarly sized cell array.
%   Alternatively, XI and OMEGA may be cell arrays and LAT a constant.
%
%   'ellrossby --t' runs a test.
%
%   Usage: ro=ellrossby(lat,xi,omega);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2021 J.M. Lilly --- type 'help jlab_license' for details
 

%   Older definition
%   For circular eddies, this definition of Rossby number is equivalent to
%   RO=-2V/(RF), where R is the radial distance to the eddy center and V is
%   the signed instantaneous azimuthal velocity. 


if strcmpi(lat, '--t')
    ellrossby_test,return
end

if ~isempty(omega)
    if ~iscell(omega)
        Ro=ellrossby_one(lat,xi,omega);
        Ro(isinf(lat.*omega))=inf;
    else
        for i=1:length(omega)
            if iscell(lat)
                lati=lat{i};
            else
                lati=lat;
            end
            Ro{i,1}=ellrossby_one(lati,xi{i},omega{i});
            Ro{i,1}(isinf(lati.*omega{i}))=inf;
        end
    end
else
    Ro=omega;
end

function[Ro]=ellrossby_one(lat,xi,omega)

fcor=(corfreq(lat))*24;  %Coriolis frequency in radians per day at center 
Ro=frac(2*omega,xi.*fcor);

%older definition
%Ro=2*sign(xi).*frac(omega,fcor); %Rossby number under solid-body assumption
    
    
function[]=ellrossby_test
 

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
[wxr,wyr,ir,jr,fr]=ridgewalk(dt,wx,wy,fs,sqrt(be*ga),1);   

%Map into time series locations
[wxr,wyr,fr]=ridgemap(length(cx),wxr,wyr,fr,ir);

[kappa,lambda,theta,phi]=ellparams(wxr,wyr);
xi=sign(lambda).*sqrt(1-squared(lambda));
ro=ellrossby(lat,xi,fr);

rm=ellrad(kappa,xi);
vm=ellvel(24*3600,kappa,xi,theta,phi,1e5);
om2=abs(2*vm./rm*frac(24*3600,100*1000));
ro2=om2./corfreq(lat)/24./xi;
%ro3=2*fr./corfreq(lat)/24./xi;

%This worked with old definition, doesn't work anymore
%reporttest('ELLROSSBY matches 2V/(Rf) form to within 2% for EBASN anticyclone', sqrt(vmean(squared(ro-ro2),1))./abs(vmean(ro,1))<0.02)


