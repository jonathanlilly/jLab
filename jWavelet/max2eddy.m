function[varargout]=max2eddy(varargin)
%MAX2EDDY  Converts transform maxima into oceanic coherent eddy properties.
%
%   This function is part of 'element analysis' described in Lilly (2017), 
%   "Element analysis: a wavelet-based method for analyzing time-localized
%   events in noisy time series", available at www.jmlilly.net.
%  
%   [A,R,V,RO]=MAX2EDDY(DX,LAT,C,RHO) converts the wavelet transform 
%   maxmima properties C and RHO, at latitudes LAT and with alongtrack 
%   sample intervals of DX, into estimated eddy properties A, R, V, and Ro.
%
%   C and RHO are as output by MAXPROPS.
%
%   The real part of A is the *apparent* coherent eddy displacement in 
%   centimeters, while R is the apparent eddy radius in kilometers.  The 
%   real part of V is a velocity estimate formed from A and R, and the 
%   real part of Ro is a similarly formed Rossby number estimate.
%
%   Note that the imaginary parts of A, V, and Ro reflect the magnitude of
%   the locally odd, or sin-like, portion of the anomaly.  These will 
%   vanish for the noise-free transform of a Gaussian eddy.
%
%   The eddy properties are derived assuming a Gaussian form for the eddy.
%
%   For details, see Lilly (2017).
%   
%   See also TRANSMAX, TRANSMAXDIST, ISOMAX, MAXPROPS.
%
%   'max2eddy --t' runs a test.
%
%   Usage: [A,R,V,Ro]=max2eddy(dx,lat,C,rho);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2019--2023 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    max2eddy_test,return
end

dx=varargin{1};
lat=varargin{2};
C=varargin{3};
rho=varargin{4};

%A=frac(C,2*sqrt(pi));
A=frac(C,4*sqrt(pi));%additional factor of two is due to recent change in normalization
R=sqrt(2)*rho.*dx;
Ro=-1*frac(9.81*A/100,sqrt(exp(1)).*squared(corfreq(lat)/3600.*(R*1000)));
V =-1*frac(9.81*A/100,sqrt(exp(1)).*(corfreq(lat)/3600.*(R*1000)));

varargout{1}=A;
varargout{2}=R;
varargout{3}=V;
varargout{4}=Ro;

function[]=max2eddy_test

dx=1;
x=dx*[-1000:1000]';
y=exp(-frac(1,2)*squared(x/50));

ga=2;be=1;
fs=morsespace(ga,be,{0.2,pi},{3,length(y)},8);
w=wavetrans(y,{ga,be,fs},'mirror');
[index,ww,ff]=transmax(fs,w);
[C,rho,frho]=maxprops(ww,ff,ga,be,0);

[A,R,~,Ro]=max2eddy(dx,45,C,rho);
reporttest('MAX2EDDY matches test Gaussian with BETA=1',aresame(real(A),1,0.05)&&aresame(R,50,0.5))

ga=2;be=2;
fs=morsespace(ga,be,{0.2,pi},{3,length(y)},8);
w=wavetrans(y,{ga,be,fs},'mirror');
[index,ww,ff]=transmax(fs,w);
[C,rho,frho]=maxprops(ww,ff,ga,be,0);

[A,R,~,Ro]=max2eddy(dx,45,C,rho);

reporttest('MAX2EDDY matches test Gaussian with BETA=2',aresame(real(A),1,1e-2)&&aresame(R,50,2.5e-2))

ga=2;be=3;
fs=morsespace(ga,be,{0.2,pi},{3,length(y)},8);
w=wavetrans(y,{ga,be,fs},'mirror');
[index,ww,ff]=transmax(fs,w);
[C,rho,frho]=maxprops(ww,ff,ga,be,0);

[A,R,~,Ro]=max2eddy(dx,45,C,rho);
reporttest('MAX2EDDY matches test Gaussian with BETA=3',aresame(real(A),1,1e-2)&&aresame(R,50,3e-2))


