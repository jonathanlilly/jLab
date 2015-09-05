function[psi]=morsexpand(varargin)
%MORSEXPAND  Generalized Morse wavelets via time-domain Taylor series.
%
%   PSI=MORSEXPAND(T,GAMMA,BETA,FS) computes the lowest-order generalized
%   Morse wavelet specified by parameters GAMMA and BETA, and having a peak
%   frequency at *radian* frequency FS, using a time-domain Taylor series.
%
%   Normally one would use MORSEWAVE to compute the wavelets in the
%   frequency domain, but it is sometimes useful to have an explicit
%   time-domain representation.
%
%   There are no contraints on the size of T nor its ordering.  
%
%   All input parameters other than T should be scalars.
%   __________________________________________________________________
%
%   Algorithms
%
%   MORSEXPAND can use one of two different algorithms.
%
%   MORSEXPAND(,...,'moment'), the default behavior,  Taylor-expands the
%   time-domain wavelet directly.
%
%   MORSEXPAND(,...,'cumulant') expands the logarithm of the time-domain 
%   wavelet and then takes the exponential. 
%
%   MORSEXPAND(N,...) specifies the number of terms to include in the
%   Taylor expansion.  N defaults to 10 in the 'cumulant' case and 100 in
%   the 'moment' case.
%
%   These algorithms implement equations (30) and (31) of Lilly and
%   Olhede (2009).
%
%   For details see
%
%       Lilly and Olhede (2009).  Higher-order properties of analytic 
%           wavelets.  IEEE Trans. Sig. Proc., 57 (1), 146--160.
%   __________________________________________________________________
%
%   See also MORSEWAVE, WAVETRANS.
%   
%   'morsexpand --t' runs a test.
%
%   Usage: psi=morsexpand(t,gamma,beta,fs);
%          psi=morsexpand(n,t,gamma,beta,fs);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2015 J.M. Lilly --- type 'help jlab_license' for details
 

%   To remove effects of the Taylor series expansion leading to incorrect
%   values far from the wavelet center, any coefficients exceeding the
%   central maximum value by more than five percent are set to NaNs.

if strcmpi(varargin(1), '--t')
    morsexpand_test,return
end

str='mom';
if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
end

if length(varargin{1})==1
    nmax=varargin{1};
    varargin=varargin(2:end);
else 
    nmax=100;
end

t=varargin{1};
ga=varargin{2};
be=varargin{3};
%if anyany(be==0)
%    error('Sorry, BETA must be greater than zero.')
%end
if length(varargin)>3
    fs=varargin{4};
end
psi=zeros(size(t,1),size(t,2),nmax+1);
if be~=0
    s=morsefreq(ga,be)./fs;
else
    s=1;
end

if strcmpi(str(1:3),'mom')
    for n=0:nmax;
        psi(:,:,n+1)=frac((sqrt(-1)*(t./s)).^n,factorial(n)).*morsemom(n,ga,be);
    end
    psi=frac(1,s).*vsum(psi,3);
elseif strcmpi(str(1:3),'cum')
    [mn,nn,kn]=morsemom([0:nmax],ga,be);
    for n=0:nmax;
        psi(:,:,n+1)=frac((sqrt(-1)*(t./s)).^n,factorial(n)).*kn(n+1);
    end
    psi=frac(1,s).*exp(vsum(psi,3));
end

psi0=frac(1,s).*morsemom(0,ga,be);
index=find(abs(psi)>abs(psi0)*1.05);
if ~isempty(index)
    psi(index)=nan;
end

function[]=morsexpand_test
morsexpand_test_frequency

function[]=morsexpand_test_frequency
t=(-500:500)';
fs=2*pi/20; 
gamma=3;
beta=6;

psi=morsexpand(t,gamma,beta,fs,'mom');
psi2=morsewave(length(t),1,gamma,beta,fs,'bandpass');

dabs=vdiff(abs(psi),1);
index1=find(dabs.*t>0|isnan(dabs));
if ~isempty(index1)
    psi(index1)=nan;
end
psi([1 end])=nan;
err=vsum(abs(psi-psi2).^2,1);

tol=1e-6;
reporttest('MORSEXPAND moment versus frequency-domain definition',err<tol)
%figure,uvplot(psi),ylim([-1 1]/4),hold on, uvplot(psi2,'g')

psi2=morsexpand(t,gamma,beta,fs,'cum');
err=vsum(abs(psi-psi2).^2,1);
reporttest('MORSEXPAND moment vs. cumulant',err<tol)


