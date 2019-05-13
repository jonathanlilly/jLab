function[k2,l2,theta2,phi2]=elldiff(varargin)
%ELLDIFF  Differentiation of modulated elliptical signals.
%
%   Given the properties of a time-varying elliptical signal, ELLDIFF 
%   finds the ellipse properties of the time derivative of that signal.  
% 
%   [K2,L2,TH2,PHI2]=ELLDIFF(K,L,TH,PHI) where the input arguments
%   specify a time-varying ellipse of the complex-valued time series 
%   Z=X+iY, returns the ellipse parameters for the associated time-
%   varing ellipse of the first central difference of Z.  
%
%   K and L are the input ellipse amplitude and linearity, respectively,
%   while TH and PHI are its orientation and orbital phase.  The output
%   arguments are the same parameters for the differentiated ellipse.
%    
%   [...]=ELLDIFF(DT,...) optionally specifies a time-step of DT
%   (default=1) for the differentiation. 
%
%   [...]=ELLDIFF(...,FACT) optionally multiplies the output amplitudes
%   by the scale factor FACT (default=1), e.g. for a unit conversion from 
%   kilometers into centimeters.
%
%   See also ELLPARAMS.
% 
%   Usage:  [k2,l2,theta2,phi2]=elldiff(k,l,theta,phi);
%           [k2,l2,theta2,phi2]=elldiff(4*3600,k,l,theta,phi,1e5);
%
%   'elldiff --t' runs a test
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2005--2015 J.M. Lilly --- type 'help jlab_license' for details    

  
if strcmpi(varargin, '--t')
  elldiff_test,return
end

%/********************************************************
%Sort out input arguments

if length(varargin{1})==1
    dt=varargin{1};
    varargin=varargin(2:end);
else 
    dt=1;
end
if length(varargin{end})==1
    fact=varargin{end};
    varargin=varargin(1:end-1);
else 
    fact=1;
end
k=varargin{1};
l=varargin{2};

[a,b]=kl2ab(k,l);

theta=varargin{3};
phi=varargin{4};

[p,n,pphi,nphi]=ellconv_ab2pn(a,b,theta,phi);
[p2,n2,pphi2,nphi2]=elldiff_pn(p,n,pphi,nphi,dt,fact);
[a2,b2,theta2,phi2]=ellconv_pn2ab(p2,n2,pphi2,nphi2);

[k2,l2]=ab2kl(a2,b2);

function[P2,N2,phip2,phin2]=elldiff_pn(P,N,phip,phin,dt,fact)

str='endpoint';

eps=1e-10;
vswap(P,0,eps);
vswap(N,0,eps);
%vsize(P,phip,dt,fact)

phip=unwrap(phip);
phin=unwrap(phin);

P2=fact*P.*sqrt(squared(frac(1,dt).*vdiff(log(P),1,str))+squared(frac(1,dt).*vdiff(phip,1,str)));
N2=fact*N.*sqrt(squared(frac(1,dt).*vdiff(log(N),1,str))+squared(frac(1,dt).*vdiff(phin,1,str)));
phip2=unwrap(phip+imlog(frac(1,dt).*vdiff(log(P),1,str)+sqrt(-1).*(frac(1,dt).*vdiff(phip,1,str))));
phin2=unwrap(phin+imlog(frac(1,dt).*vdiff(log(N),1,str)+sqrt(-1).*(frac(1,dt).*vdiff(phin,1,str))));


function[P,N,phip,phin]=ellconv_ab2pn(A,B,theta,phi)  
P=frac(A+B,2);
N=frac(A-B,2);

phip=unwrap(phi+theta);
phin=unwrap(phi-theta);

function[A,B,theta,phi]=ellconv_pn2ab(P,N,phip,phin)  
A=P+N;
B=P-N;

theta=unwrap(phip/2-phin/2);
phi=  unwrap(phip/2+phin/2);

function[]=elldiff_test
a=1;
b=(0:.05:.095);
phi=(0:.01:2*pi-0.01)';

[k,l]=ab2kl(a+0*b,b);
clear xi1 xi2
for i=1:length(l)
   [k2,l2,th2,phi2]=elldiff(k(i)+0*phi,l(i)+0*phi,0*phi,phi);
   xi1(:,i)=vdiff(ellsig(k(i),l(i),0,phi),1);
   xi2(:,i)=ellsig(k2,l2,th2,phi2);
end

%Put nans in same places
index=find(isnan(xi1)|isnan(xi2));
xi1(index)=nan;
xi2(index)=nan;

reporttest('ELLDIFF',aresame(xi1(2:end-1,:),xi2(2:end-1,:),1e-6));


