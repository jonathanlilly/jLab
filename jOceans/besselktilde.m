function[K]=besselktilde(varargin)
%BESSELKTILDE  K-type Bessel function after factoring off exponential decay.
%
%   KTILDE=BESSELKTILDE(NU,Z) returns the modified Bessel function of the
%   second kind of order NU at argument Z, after factoring off the 
%   asymptotic behavior of EXP(-Z).
%
%   BESSELKTILDE is useful for products of modified Bessel functions in
%   which the exponential behaviors cancel, but that cannot be evaluated
%   directly because of numerical overflow.
%
%   BESSELITILDE(NU,Z,N) use a summation truncated at N terms. The default
%   behavior uses N=30 and is highly accurate. 
%
%   See 10.40.2 of https://dlmf.nist.gov/10.40.
%
%   This is low-level code used by WINDTRANS using an algorithm described 
%   in Lilly and Elipot (2021).
%
%   See also BESSELITILDE.
%
%   'besselktilde --t' runs a test.
%
%   Usage: K=besselktilde(nu,z);
%          K=besselktilde(nu,z,nterms);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2019--2021 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    besselktilde_test,return
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
zk=vrep(z,nterms,2);  
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

K=sqrt(frac(pi,2*z)).*((1./zk)*ak);
K=reshape(K,sizez);

%nterms
% %let's sum over the third dimensions
% zk=vrep(z,nterms,3);
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
% K=sqrt(frac(pi,2*z)).*sum(frac(ak,zk),3);

function[]=besselktilde_test
 
for s=[1 -1]
    z=sqrt(s*1i)*[15:0.01:100]';
    %z=[20:0.01:100]';
    [bk0,bk]=vzeros(length(z),2);
    for i=1:2
        bk0(:,i)=exp(z).*besselk(i-1,z);
        bk(:,i)=besselktilde(i-1,z);
    end
    
    if s==1
        reporttest('BESSELKTILDE for z with phase of pi/4',allall(abs((bk0-bk)./bk)<1e-14))
    else
        reporttest('BESSELKTILDE for z with phase of -pi/4',allall(abs((bk0-bk)./bk)<1e-14))
    end
    
	bk=vzeros(length(z),2);
    for i=1:2
        bk(:,i)=besselktilde(i-1,z,2);
    end
  
    if s==1
        reporttest('BESSELKTILDE 2-term for z wkth phase of pi/4, order 0 and 1',allall(abs((bk0(:,1:2)-bk)./bk)<1e-3))
    else
        reporttest('BESSELKTILDE 2-term for z wkth phase of -pi/4, order 0 and 1',allall(abs((bk0(:,1:2)-bk)./bk)<1e-3))
    end
      
    bk0=sqrt(frac(pi,2*z)).*(1-frac(1,8*z));
    bk1=sqrt(frac(pi,2*z)).*(1+frac(3,8*z));    
    
    if s==1
        reporttest('BESSELKTILDE 2-term vs. analytic for z with phase of pi/4',allall(abs(([bk0 bk1]-bk)./bk)<1e-15))
    else
        reporttest('BESSELKTILDE 2-term  vs. analytic for z with phase of -pi/4',allall(abs(([bk0 bk1]-bk)./bk)<1e-15))
    end
end

% z=sqrt(s*1i)*[17:0.01:100]';
%tic;bk=besselktilde(0,z,2);toc
%tic;bk=besselktilde(0,z,20);toc
%tic;bk=besselktilde(0,z,20);toc
%tic;bk=besselk(0,z);toc

        
% 
% z=sqrt(1i)*[1:0.01:100]';
% [bk0,bk]=vzeros(length(z),2);
% for i=1:2
%     bk0(:,i)=exp(z).*besselk(i-1,z);
%     bk(:,i)=besselktilde(i-1,z);
% end
% figure,plot(abs(z),abs(bk-bk0)),ylog

% for s=[1 -1]
%     z=sqrt(s*1i)*[18:0.01:100]';
%     [bk0,bk]=vzeros(length(z),2);
%     for i=1:2
%         bk0(:,i)=exp(z).*besselk(i-1,z);
%         bk(:,i)=besselktilde(i-1,z);
%     end
%     maxmax(abs(bk-bk0))<1e-15
% end

%reporttest('BESSELKTILDE',aresame())

% nu=5;
% z=[0:0.001:100]';
% tic;bk(:,1)=besselktilde(nu,z);toc
% tic;bk(:,2)=sqrt(frac(2*z,pi)).*exp(z).*besselk(nu,z);toc
% figure,plot(z,log10(abs((bk(:,1)-bk(:,2))./bk(:,2))))
% %wow, really excellent
% 
% 
% tic;bk0=sqrt(frac(2*z,pi)).*exp(z).*besselk(nu,z);toc
% for i=1:50
%     bk(:,i)=besselktilde(nu,z,i);
% end
% bk0=vrep(bk0,size(bk,2),2);
% plot(z,log10(abs((bk-bk0)./bk0)))
% 
% 



