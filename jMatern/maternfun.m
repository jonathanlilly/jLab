function[M]=maternfun(alpha,r)
%MATERNFUN  Returns the Matern function.  
%
%   MATERNFUN is a low-level function called by MATERNCOV.
%
%   M=MATERNFUN(ALPHA,X) returns the value of the Matern function of order
%   ALPHA at location X, with ALPHA>1/2.  The Matern function is defined as
%
%      M(X)=C*|X|^(ALPHA-1/2) K_(ALPHA-1/2)(|X|)
%   
%   where K_ALPHA is the modified Bessel function of the second kind of 
%   order ALPHA, and C=2/GAMMA(ALPHA-1/2)./2^(ALPHA-1/2) is a normalizing 
%   constant. This definition implies M(0)=1.
%
%   For the special case of ALPHA=1/2, the GAMMA function in C is omitted
%   in order to avoid developing a singularity. This change is acceptable 
%   because this function is used only in ratios for ALPHA<=1/2.  
%
%   'maternfun --f' generates a sample figure.
%
%   Usage:  M=maternfun(alpha,x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2017 J.M. Lilly --- type 'help jlab_license' for details
 
%   ALPHA and R can either be scalars or 1D arrays.  XI will be an array of
%   size LENGTH(R) x LENGTH(ALPHA).
%
%   [XI,XI1,XI2]=MATERNFUN(ALPHA,R) also returns its first two derivatives.
%

if strcmpi(alpha, '--f')
    maternfun_fig,return
end

if strcmpi(alpha, '--t')
    maternfun_test,return
end

r=abs(r);
if alpha==1/2
    fact=frac(2,pow2(abs(alpha-1/2)));
else
    fact=frac(2,gamma(abs(alpha-1/2)).*pow2(abs(alpha-1/2)));
end
M=fact.*(r.^(alpha-1/2)).*besselk(abs(alpha-1/2),r);
if alpha>1/2
    M(r==0)=1;   %fix for numerical problem at r=0
    %See 9.6.9 of Abramowitz and Stegun for the small-argument behavior of K
end

function[]=maternfun_fig
r1=[0:.01:10]';
alpha1=[1/2:0.1:4];
[alpha,r]=meshgrid(alpha1,r1);
xi=maternfun(alpha,r);
figure,plot(r1,xi)
hold on,h=plot(r1,xi(:,1)); linestyle -h h 3D
hold on,h=plot(r1,xi(:,end)); linestyle -h h 3k
title('The Matern function from \alpha = 1/2 (gray) to \alpha = 4 (black)')

function[]=maternfun_test
alpha=1/2;
r=[0:.01:10]';
xi=maternfun(-alpha,r);xi2=maternfun(alpha+1,r);
reporttest('MATERNFUN reflection relation for ALPHA = 1/2',aresame(xi,xi2./r.^(2*alpha+1),1e-10))

alpha=3/4;
r=[0:.01:10]';
xi=maternfun(-alpha,r);xi2=maternfun(alpha+1,r);
reporttest('MATERNFUN reflection relation for ALPHA = 3/4',aresame(1+0*xi,xi2./r.^(2*alpha+1)./xi,1e-10))

alpha=1;
r=[0:.01:10]';
xi=maternfun(-alpha,r);xi2=maternfun(alpha+1,r);
reporttest('MATERNFUN reflection relation for ALPHA = 1',aresame(1+0*xi,xi2./r.^(2*alpha+1)./xi,1e-10))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Old stuff
    
% if size(alpha,1)>1
%     error('MATERNFUN is expecting a scalar or a row vector for ALPHA.')
% end


%r(abs(r)<1e-10)=1e-10;
%[alpha,r]=arrayify(alpha(:)',r(:));

%fact=frac(2,gamma(alpha-1/2).*pow2(alpha-1/2));
%fact=1;

%if length(alpha)==1
    %if alpha~=0
%        M=fact.*(r.^(alpha-1/2)).*besselk(abs(alpha-1/2),r);
    %else
    %    xi=fact.*exp(-r);
    %end
%else
%    xi=fact.*(r.^(alpha-1/2)).*besselk(abs(alpha-1/2),r);
%    xi(alpha==0)=(fact(alpha==0).*exp(-r(alpha==0)));
%end
    


% %xi=real(xi);
% if nargout>1
%     xi1=-fact.*(r.^alpha).*besselk(alpha-1,r);
% end
% if nargout>2
%     xi2=fact.*(r.^alpha).*besselk(alpha-2,r)-fact.*(r.^(alpha-1)).*besselk(alpha-1,r);
% end

%function[]=maternfun_test

%alpha=[1/2+1/10 3/4:1/4:10];
%[xi,xi1,xi2]=maternfun(alpha,maternrad(alpha));
%reporttest('MATERNFUN second derivative vanishes at MATERNRAD',xi2(1)./xi(1)<1e-2&&allall(xi2(2:end)./xi(2:end)<1e-3))

%These are for an older definition of maternfun, differing by 1/2 alpha
% r1=[0:.01:10]';
% alpha1=[1/2:0.1:4];
% [alpha,r]=meshgrid(alpha1,r1);
% [xi,xi1]=maternfun(alpha,r);
% dxi=vdiff(xi,1)./0.01;
% reporttest('MATERNFUN iterative expression for derivative matches derivative',aresame(xi1(2:end,:),dxi(2:end,:),0.03))
% 
% r1=[0:.01:10]';
% alpha1=[1/2:0.1:0.9 1.1:0.1:4];
% [alpha,r]=meshgrid(alpha1,r1);
% [xi,xi1]=maternfun(alpha,r);
% dxi=-frac(r,2*(alpha-1)).*maternfun(alpha-1,r);
% reporttest('MATERNFUN iterative expression for derivative matches direct expression',aresame(xi1(2:end,:),dxi(2:end,:),1e-10))

