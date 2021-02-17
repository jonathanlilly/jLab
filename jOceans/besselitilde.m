function[I]=besselitilde(varargin)
%BESSELITILDE  I-type Bessel function after factoring off exponential growth.
%
%   ITILDE=BESSELITILDE(NU,Z) returns the modified Bessel function of the
%   first kind of order NU at argument Z, after factoring off the  
%   asymptotic behavior of EXP(Z).
%
%   BESSELITILDE is useful for products of modified Bessel functions in
%   which the exponential behaviors cancel, but that cannot be evaluated
%   directly because of numerical overflow.
%
%   BESSELITILDE(NU,Z,N) use a summation truncated at N terms. The default
%   behavior uses N=30 and is highly accurate.
%
%   See 10.40.1 of https://dlmf.nist.gov/10.40.
%
%   This is low-level code used by WINDTRANS using an algorithm described 
%   in Lilly and Elipot (2021).
%
%   See also BESSELKTILDE.
%
%   'besselitilde --t' runs a test.
%
%   Usage: I=besselitilde(nu,z);
%          I=besselitilde(nu,z,nterms);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2019--2021 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    besselitilde_test,return
end

nu=varargin{1};
z=varargin{2};
nterms=30;
if nargin==3
    nterms=varargin{3};
end
 
sizez=size(z);
z=z(:);
%let's sum over the third dimensions
zk=vrep(-z,nterms,2);  %this makes the alternating sign
zk(:,1)=1;
zk=cumprod(zk,2);

k=(0:nterms-1)';
ak=(4*nu.^2-(2*k-1).^2);
ak(1)=1;
ak=cumprod(ak,1);
ak=frac(ak,factorial(k).*8.^k);
if any(~isfinite(ak))
    ak=ak(1:find(~isfinite(ak),1,'first')-1);
end
zk=zk(:,1:length(ak));

I=frac(1,sqrt(2*z*pi)).*((1./zk)*ak);
I=reshape(I,sizez);
%nterms

% %let's sum over the third dimensions
% zk=vrep(-z,nterms,3);  %this makes the alternating sign
% zk(:,:,1)=1;
% zk=cumprod(zk,3);
% 
% k=(0:nterms-1)';
% ak=(4*nu.^2-(2*k-1).^2);
% ak(1)=1;
% ak=cumprod(ak,1);
% ak=frac(ak,factorial(k).*8.^k);
% if any(~isfinite(ak))
%     ak=ak(1:find(~isfinite(ak),1,'first')-1);
% end
% zk=zk(:,:,1:length(ak));
% 
% ak=vrep(permute(ak,[3 2 1]),size(z),[1 2]);
% 
% I=frac(1,sqrt(2*z*pi)).*sum(frac(ak,zk),3);


function[]=besselitilde_test


for s=[1 -1]
    z=sqrt(s*1i)*[23:0.01:100]';
    %z=[20:0.01:100]';
    [bi0,bi]=vzeros(length(z),2);
    for i=1:2
        bi0(:,i)=exp(-z).*besseli(i-1,z);
        bi(:,i)=besselitilde(i-1,z);
    end
    
    if s==1
        reporttest('BESSELITILDE for z with phase of pi/4',allall(abs((bi0-bi)./bi)<1e-14))
    else
        reporttest('BESSELITILDE for z with phase of -pi/4',allall(abs((bi0-bi)./bi)<1e-14))
    end
    
	bi=vzeros(length(z),2);
    for i=1:2
        bi(:,i)=besselitilde(i-1,z,2);
    end
  
    if s==1
        reporttest('BESSELITILDE 2-term for z with phase of pi/4, order 0 and 1',allall(abs((bi0(:,1:2)-bi)./bi)<1e-3))
    else
        reporttest('BESSELITILDE 2-term for z with phase of -pi/4, order 0 and 1',allall(abs((bi0(:,1:2)-bi)./bi)<1e-3))
    end
      
    bi0=sqrt(frac(1,2*pi*z)).*(1+frac(1,8*z));
    bi1=sqrt(frac(1,2*pi*z)).*(1-frac(3,8*z));    
    
    if s==1
        reporttest('BESSELITILDE 2-term vs. analytic for z with phase of pi/4',allall(abs(([bi0 bi1]-bi)./bi)<1e-15))
    else
        reporttest('BESSELITILDE 2-term  vs. analytic for z with phase of -pi/4',allall(abs(([bi0 bi1]-bi)./bi)<1e-15))
    end
end

% figure,plot(abs(z),abs(bi-bi0)),ylog


% nu=5;
% z=[0:0.001:100]';
% tic;bi(:,1)=besselitilde(nu,z);toc
% tic;bi(:,2)=sqrt(2*z*pi).*exp(-z).*besseli(nu,z);toc
% figure,plot(z,log10(abs((bi(:,1)-bi(:,2))./bi(:,2))))
% %wow, really excellent
% 
% 
% tic;bi0=sqrt(2*z*pi).*exp(-z).*besseli(nu,z);toc
% for i=1:50
%     bi(:,i)=besselitilde(nu,z,i);
% end
% bi0=vrep(bi0,size(bi,2),2);
% plot(z,log10(abs((bi-bi0)./bi0)))
% 
% 




