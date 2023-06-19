function[mu2,nu0,rho,nu2,mu4,nu4]=pm_kernmom(al,be,d)
%PM_KERNMOM  Moments of smoothing kernels for local polynomial fitting.
%   
%   [MU2,NU0]=PM_KERNMOM(ALPHA,BETA) returns the moments of the generalized 
%   beta family of smoothing kernels in 2D. These are of the form
%
%        K(R) = KAPPA * (1-R^ALPHA)^BETA,   R<1
%    
%   and are defined to vanish for R>=1.   MU2 is the second moment of K(R),
%   while NU0 is the first moment of the squared kernel K^2(R).
%
%   [MU2,NU0,RHO]=PM_KERNMOM(ALPHA,BETA) also returns RHO, a combination of
%   higher-order moments that appears in the variance of quadratic fits.
%
%   [MU2,NU0,RHO,NU2,MU4,NU4]=PM_KERNMOM(ALPHA,BETA) also returns further
%   higher-order moments.
%
%   'pm_kernmom --t' runs some tests.
%
%   Usage: [mu2,nu0]=pm_kernmom(alpha,beta);
%          [mu2,nu0,rho]=pm_kernmom(alpha,beta);
%          [mu2,nu0,rho,nu2,mu4,nu4]=pm_kernmom(alpha,beta);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2022 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(al, '--t')
    pm_kernmom_test,return
end

if nargin==2
    d=2;
end

[mu0,nu0]=pm_kernmom_one(al,be,d,0);
[mu2,nu2]=pm_kernmom_one(al,be,d,2);

if nargout>3
    [mu4,nu4]=pm_kernmom_one(al,be,d,4);
    numer=2.*mu4.^2-6*mu2.*mu4.*(nu2./nu0)+3*mu2.^2.*(nu4./nu0);
    rho=2*frac(numer,squared(3*mu2.^2-2*mu4));
end

function[mu,nu]=pm_kernmom_one(al,be,d,n)
if n==0
    A=2*pi*1;
elseif n==2
    A=2*pi*1/2;
elseif n==4
    A=2*pi*3/8;
end

kappa=frac(1,2*pi)*frac(al,beta(2./al,be+1));
if d==2
    mu=A.*frac(kappa,al)*beta(frac(n+2,al),be+1);
    nu=A.*frac(kappa.^2,al)*beta(frac(n+2,al),2*be+1);
end

function[A]=kern_A(n)
A=2*pi*frac(factorial(n),2.^(n).*squared(factorial(n/2)));

function[mu2,nu0,nu2,mu4,nu4]=pm_kernmom_numeric(al,be)
dr=0.00001;
r=[0:dr:1]';

kappa=frac(1,2*pi*sum(r.*((1-r.^al).^be)*dr));

kern=((1-r.^al).^be);
mat(:,1)=kappa*kern_A(2)*kern.*r.^3*dr;
mat(:,2)=kappa.^2*kern_A(0)*kern.^2.*r*dr;
mat(:,3)=kappa.^2*kern_A(2)*kern.^2.*r.^3*dr;
mat(:,4)=kappa*kern_A(4)*kern.*r.^5*dr;
mat(:,5)=kappa.^2*kern_A(4)*kern.^2.*r.^5*dr;

mat=sum(mat,1);
mu2=mat(:,1);
nu0=mat(:,2);
nu2=mat(:,3);
mu4=mat(:,4);
nu4=mat(:,5);


function[]=pm_kernmom_test

mom=zeros(1,6);
[mom(1) mom(2) mom(3) mom(4) mom(5) mom(6)]=pm_kernmom(2,1);
mom2=[1/6 4/3/pi 3.6 4/3/pi/8 1/16 4/3/pi*3/80];


tol=1e-10;
reporttest('PM_KERNMOM parabolic kernel',aresame(mom,mom2,tol))

[mu2,nu0,nu2,mu4,nu4]=pm_kernmom_numeric(2,1);
mom3=[mu2,nu0,nu2,mu4,nu4];
reporttest('PM_KERNMOM numeric integration with alpha=2, beta=1',...
    aresame(mom([1 2 4 5 6]),mom3,tol))

[mom(1) mom(2) mom(3) mom(4) mom(5) mom(6)]=pm_kernmom(2,2);
[mu2,nu0,nu2,mu4,nu4]=pm_kernmom_numeric(2,2);
mom3=[mu2,nu0,nu2,mu4,nu4];
reporttest('PM_KERNMOM numeric integration with alpha=2, beta=2',...
    aresame(mom([1 2 4 5 6]),mom3,tol))


%reporttest('PM_KERNMOM',aresame())
%    mu2=frac(1,2)*frac(beta(4./al,be+1),beta(2./al,be+1));
%    nu0=frac(al,2*pi)*frac(beta(2./al,2*be+1),squared(beta(2./al,be+1)));
%    mu2=frac(1,d).*frac(beta((d+2)./al,be+1),beta(d./al,be+1));
%    nu0=frac(al*gamma(d/2),2*pi.^(d/2))*frac(beta(d./al,2*be+1),squared(beta(d./al,be+1)));


