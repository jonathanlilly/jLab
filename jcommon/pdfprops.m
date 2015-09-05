function[mu,sigma,skew,kurt]=pdfprops(x,fx,dim)
%PDFPROPS  Mean and variance associated with a probability distribution.
%
%   [MU,SIGMA]=PDFPROPS(X,FX), given a probability distribution function FX
%   over values X, returns the mean MU and the standard deviation SIGMA. 
%
%   The statistics are computed using a trapezoidal integration.  FX is 
%   multiplied by a constant so that it integrates to one.
%
%   [MU,SIGMA,SKEW,KURT]=PDFPROPS(X,FX) also retuns the skewness and the
%   kurtosis, which are the third and fourth central moments, respectively 
%   normalized by the third and fourth powers of the standard deviation.  
%
%   PDFPROPS(X,FX,DIM) alternately computes moments along dimension DIM.
%   By default, DIM=1, so the moments are computed along rows.
%   
%   'pdfprops --t' runs a test.
%
%   Usage:  [mu,sigma]=pdfprops(x,fx); 
%           [mu,sigma,skew,kurt]=pdfprops(x,fx);
%           [mu,sigma,skew,kurt]=pdfprops(x,fx,dim);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2013 J.M. Lilly --- type 'help jlab_license' for details    
  
if strcmpi(x,'--t')
    pdfprops_test;
    return
end

if nargin==2
    dim=1;
end

dx=vindex(x,2,dim)-vindex(x,1,dim);
%dx=x(2,:)-x(1,:);

fx=fx./vrep(trapint(fx,dx,dim),size(x,dim),dim);
mu=real(trapint(fx.*x,dx,dim));
%mu=trapint(fx.*x,dx);
murep=vrep(mu,size(x,dim),dim);
sigma=sqrt(trapint((x-murep).^2.*fx,dx,dim));
if nargout>=3
   skew=trapint((x-murep).^3.*fx,dx,dim);
   skew=skew./sigma.^3;
end
if nargout==4
   kurt=trapint((x-murep).^4.*fx,dx,dim);
   kurt=kurt./sigma.^4;
end


% for i=1:size(fx,2)
%     if trapint(fx(:,i),dx(i))~=1
%         %disp('Normalizing FX to unit area.')
%         fx(:,i)=fx(:,i)./trapint(fx(:,i),dx(i));
%     end
% end
% 
% for i=1:size(fx,2)
%   mu(i,1)=real(trapint(fx(:,i).*x(:,i),dx(i)));
%   sigma(i,1)=sqrt(trapint((x(:,i)-mu(i,1)).^2.*fx(:,i),dx(i)));
% end
% 


function[y]=trapint(f,dx,dim)
%Trapezoidal integration

fa=f;
fb=vshift(fa,1,dim);
fa=vindexinto(fa,0,1,dim);
fa=vindexinto(fa,0,size(fa,dim),dim);
fb=vindexinto(fb,0,1,dim);
fb=vindexinto(fb,0,size(fb,dim),dim);

% For the case of DIM=1, this becomes:
% fa(1,:)=0;
% fb(1,:)=0;
% fa(end,:)=0;
% fb(end,:)=0;
y=vsum(frac(fa+fb,2),dim).*dx;
vswap(y,0,1e-10);

function[]=pdfprops_test
x=(-30:.001:30)';

mu0=2;
sigma0=5;

f=simplepdf(x,mu0,sigma0,'gaussian');  %f=simplepdf(x,mu,sig,flag)  
[mug,sigmag,skewg,kurtg]=pdfprops(x,f);
[mug2,sigmag2,skewg2,kurtg2]=pdfprops(permute(x,[3 2 1]),permute(f,[3 2 1]),3);


f=simplepdf(x,mu0,sigma0,'boxcar');  %f=simplepdf(x,mu,sig,flag)  
[mu,sigma]=pdfprops(x,f);
tol=1e-3;

bool(1)=aresame(mu,mu0,tol).*aresame(sigma,sigma0,tol);
bool(2)=aresame(mug,mu0,tol).*aresame(sigmag,sigma0,tol);
bool(3)=aresame(skewg,0,tol).*aresame(kurtg,3,tol);

reporttest('PDFPROPS with uniform pdf', bool(1));
reporttest('PDFPROPS with Gaussian pdf', bool(2));
reporttest('PDFPROPS Gaussian skewness=0, kurtosis=3', bool(3));

bool=aresame([mug,sigmag,skewg,kurtg],[mug2,sigmag2,skewg2,kurtg2],1e-10);

reporttest('PDFPROPS calculation along alternate dimension', bool);



% %/********************************************************
% x=(-10:.001:10)';
% f=simplepdf(x,0,2,'gaussian');
% f(end/2:end)=2*f(end/2:end);
% f(1:end/2)=0;
% f=f./sum(f)./0.001;
% plot(x,cumsum(f*.001))
% %********************************************************
