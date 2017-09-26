function[f]=simplepdf(x,mu,sig,flag)
%SIMPLEPDF  Gaussian, uniform, Cauchy, and exponential pdfs.
%
%   F=SIMPLEPDF(X,MU,SIG,'gaussian') computes a Gaussian pdf with mean
%   MU and standard deviation SIG.
%  
%   F=SIMPLEPDF(X,MU,SIG,'boxcar') computes a uniform pdf with mean MU
%   and standard deviation SIG.
%  
%   F=SIMPLEPDF(X,XO,ALPHA,'cauchy') computes a Cauchy pdf with location
%   parameter XO and scale parameter ALPHA.
%
%   F=SIMPLEPDF(X,BETA,'exponential') computes an exponential pdf with
%   scale parameter, hence mean and standard deviation, equal to BETA.
%
%   'simplepdf --f' generates a sample figure
%
%   Usage: f=simplepdf(x,mu,sig,'gaussian');
%          f=simplepdf(x,mu,sig,'boxcar');
%          f=simplepdf(x,xo,alpha,'cauchy');
%          f=simplepdf(x,beta,'exponential');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2014 J.M. Lilly --- type 'help jlab_license' for details    
  
  
if strcmpi(x,'--f')
  type  makefigs_simplepdf
  makefigs_simplepdf;
  return
end

dx=x(2)-x(1);

if nargin==3
    flag=sig;
end

if nargin<3&&strcmpi(flag,'exponential')||nargin<4&&~strcmpi(flag,'exponential')
    error('Not enough input arguments.')
end

if strcmpi(flag,'gaussian')
  f=exp(-(x-mu).^2./2./sig.^2)./sig./sqrt(2*pi);
elseif strcmpi(flag,'boxcar')
  f=0*x;
  ia=min(find(x-mu>-3.4641*sig/2))-1;
  ib=min(find(x-mu>3.4641*sig/2));
  f(ia:ib)=1;
  f=f./vsum(f*dx,1);
elseif strcmpi(flag,'cauchy')
  alpha=sig;
  f=frac(alpha./pi,(x-mu).^2 + alpha.^2);
elseif strcmpi(flag,'exponential')
  f=frac(1,abs(mu)).*exp(-abs(x./mu));
end

