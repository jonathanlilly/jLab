function[P,alpha,beta,kbar,rbar]=ellpol(varargin)
%ELLPOL  Polarization parameters of an elliptical signal.
%
%   [P,ALPHA,BETA]=ELLPOL(KAPPA,LAMBDA,THETA,PHI) where KAPPA, LAMBDA,
%   THETA, and PHI are the time-varying parameters of an elliptical signal,
%   returns the time-averaged polarization parameters P, ALPHA, and BETA.
%
%   P, ALPHA, and BETA are properties of the frequency-integrated spectral
%   matrix associated with the elliptical signal.  P is related to the 
%   eigenvalues D1 and D2 of that matrix by (D1-D2)(D1+D2).
%
%   P gives the total polarization, ALPHA is the excess of positive 
%   to negative rotational energy, and BETA is the polarization of the real
%   part of the spectral matrix associated with linear motions.
%
%   These three quantities are related by P^2=ALPHA^2+BETA^2.
%
%   The input fields may be arrays of any dimension.  The averaging is
%   performed along rows, and the result is then squeezed.  If the 
%   input fields are column arrays, the output fields are scalars.
%   
%   If the input fields are cell arrays of column vectors, the output
%   will be column vectors with one entry per cell of the input arrays.
%
%   [P,ALPHA,BETA,KBAR,RBAR]=ELLPOL(...) also returns the average RMS axis
%   length KBAR and average geometric mean radius RBAR, given in terms of
%   the eigenvalues as DBAR=SQRT(D1^2/2+D2^2/2) and RBAR=SQRT(D1*D2).
%
%   See also POLPARAMS.
%
%   'ellpol --t' runs some tests.
%
%   Usage: [P,alpha,beta]=ellpol(kappa,lambda,theta,phi);
%          [P,alpha,beta,Kbar,Rbar]=ellpol(kappa,lambda,theta,phi);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2018 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    ellpol_test,return
end
 
kappa=varargin{1};
lambda=varargin{2};
theta=varargin{3};
phi=varargin{4};

if ~iscell(kappa)
    [P,alpha,beta,kbar,rbar]=ellpol_one(kappa,lambda,theta,phi);
else
    [P,alpha,beta,kbar,rbar]=vzeros(length(kappa),1);
    for i=1:length(kappa)
        [P(i),alpha(i),beta(i),kbar(i),rbar(i)]=ellpol_one(kappa{i},lambda{i},theta{i},phi{i});
    end
end


function[P,alpha,beta,kbar,rbar]=ellpol_one(kappa,lambda,theta,phi)

[xr,yr]=ellsig(kappa,lambda,theta,phi);
sxx=vmean(squared(xr),1);
syy=vmean(squared(yr),1);
sxy=vmean(xr.*conj(yr),1);
[d1,d2,th,nu]=specdiag(sxx,syy,sxy);
P=squeeze(frac(d1-d2,d1+d2));
%Frequency-integrated rototary ratio
alpha=squeeze(vmean(frac(2*imag(sxy),sxx+syy),1));  
[d1,d2,th,nu]=specdiag(sxx,syy,real(sxy));
beta=squeeze(frac(d1-d2,d1+d2));
kbar=sqrt(d1./2+d2./2);
rbar=sqrt(sqrt(d1.*d2));

if size(P,2)>1
    P=shiftdim(P);
    alpha=shiftdim(alpha);
    beta=shiftdim(beta);
    kbar=shiftdim(kbar);
    rbar=shiftdim(rbar);
end

function[]=ellpol_test

load npg2006
use npg2006

%Decide on frequencies
fs=2*pi./(logspace(log10(10),log10(100),50)');

%Compute wavelet transforms using generalized Morse wavelets
P=sqrt(2*4);
[wx,wy]=wavetrans(real(cx),imag(cx),{2,4,fs,'bandpass'},'mirror');
[wxr,wyr]=ridgewalk(dt,wx,wy,fs,P,3);   

[kappa,lambda,theta,phi]=ellparams(wxr,wyr);
[P,alpha,beta]=ellpol(kappa,lambda,theta,phi);
reporttest('ELLPOL P^2=ALPHA^2+BETA^2',aresame(P.^2,alpha.^2+beta.^2,1e-8))

[Po,alphao,betao]=ellpol(kappa,lambda,theta,phi);
[P,alpha,beta]=ellpol([kappa kappa],[lambda lambda],[theta theta],[phi phi]);
bool=aresame([P alpha beta],[Po alphao betao;Po alphao betao]);
reporttest('ELLPOL matrix input works as expected',bool)

[P,alpha,beta]=ellpol({kappa,kappa},{lambda,lambda},{theta,theta},{phi,phi});
bool=aresame([P alpha beta],[Po alphao betao;Po alphao betao]);
reporttest('ELLPOL cell array input works as expected',bool)
