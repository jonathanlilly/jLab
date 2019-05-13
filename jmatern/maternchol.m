function[T,C]=maternchol(varargin)
%MATERNCHOL  Cholesky decomposition of Matern and fBm covariances. [with A. Sykulski]
%
%   MATERNCHOL is a low-level function called by MATERNOISE.
%
%   T=MATERNCHOL(DT,N,SIGMA,ALPHA,LAMBDA) returns the Cholesky 
%   decomposition of the N x N autocovariance matrix of a Matern process.
%
%   T=MATERNCHOL(DT,N,SIGMA,ALPHA,LAMBDA,NU,MU) does the same, but for the
%   extended Matern process.   See MATERNCOV for details.
%
%   T=MATERNCHOL(DT,N,SIGMA,ALPHA,LAMBDA,NU,MU,'composite') does the same, 
%   but for the composite Matern process.  See MATERNCOV for details.
%
%   T=MATERNCHOL(DT,N,A,ALPHA,0) returns the Cholesky decomposition for
%   fractional Brownian motion.  In this case, the third input argument is
%   the spectral amplitude A, not the standard deviation SIGMA.
%   __________________________________________________________________
%
%   Non positive-definiteness adjustment
%
%   It can be the case that the covariance matrix is not positive definite 
%   to numerical precision, owing to very small eigenvalues, and this
%   will cause the Cholesky decomposition to fail.
%
%   In such cases, an identity matrix with a very small amplitude is added
%   to the covariance matrix in order to ensure numerical stability, and a
%   notification is issued.  The amplitude begins at 1e-16 and increases
%   by powers of ten until positive definiteness is attained.
%   __________________________________________________________________
%
%   'maternchol --t' runs a test.
%
%   Usage: T=maternchol(dt,N,sigma,alpha,lambda,nu,mu);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2017  A.M. Sykulski and J.M. Lilly
%                                 --- type 'help jlab_license' for details 
if strcmp(varargin{1}, '--t')
    maternchol_test,return
end

alg='fast';
model='forward';
for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'loo')||strcmpi(varargin{end}(1:3),'fas')
            alg=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'for')||strcmpi(varargin{end}(1:3),'rev')||strcmpi(varargin{end}(1:3),'com')
            model=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end

dt=double(varargin{1});
N=varargin{2};
sigma=varargin{3};
alpha=varargin{4};

lambda=varargin{5};
nu=0*alpha;
mu=0*alpha;

if nargin>5
    nu=varargin{6};
    mu=varargin{7};
end


%dt,N,sigma,alpha,lambda,nu,mu

%/*************************************************************************
%Error checking
if strcmpi(model(1:3),'rev')&&strcmpi(alg(1:3),'loo')
    error('Sorry, MATERNCHOL cannot accept both the LOOP and REVERSE options.')
end

if strcmpi(model(1:3),'com')&&strcmpi(alg(1:3),'loo')
    error('Sorry, MATERNCHOL cannot accept both the LOOP and COMPOSITE options.')
end

if strcmpi(model(1:3),'rev')&&(lambda==0)
    error('Sorry, MATERNCHOL cannot employ REVERSE option with LAMBDA set to zero.')
end
%\*************************************************************************



if lambda==0
    C=fbm_matrix(dt,N,sigma,alpha,alg);
else
    [tau,tpz]=materncov(dt,N,sigma,alpha,lambda,nu,mu,model);
    C=toeplitz(tpz);
end

%Explicitly remove small imaginary parts along diagonal
for i=1:size(C,1)
    C(i,i)=real(C(i,i));
end

%figure,uvplot(diag(C))
%T=conj(chol(C,'lower'));
%figure,jpcolor(abs(C))

%Note lower triangular is a causal filter.  Plot rows (not columns).
try
    T=conj(chol(C,'lower'));
catch
    isdone=false;
    eps=1e-16;
    while ~isdone
        try
            C=(C+eye(size(C))*eps)./(1+eps);
            T=conj(chol(C,'lower'));
            isdone=true;
            disp(['MATERNCHOL covariance matrix is not strictly positive definite, adding ' num2str(eps) ' noise.'])
        catch
            eps=eps*10;
        end
    end
end


function [matrix]=fbm_matrix(dt,N,sigma,alpha,alg)

% This is coded up as in Barton and Poor (1988) 
%   "Signal Detection in Fractional Gaussian Noise"
% IEEE TRANSACTIONS ON INFORMATION THEORY, VOL. 34, NO. 5, SEPTEMBER 1988

% Note the Hurst parameter H = alpha  - 1/2
alpha(alpha==0.5)=0.5+1e-5;
alpha(alpha==1.5)=1.5-1e-5;
alpha(alpha==1)=1+1e-5;

%H = alpha - 1/2-sqrt(-1)*nu;
lambda = alpha - 1/2;
M=length(alpha);
t=dt*[1:N]';

if strcmpi(alg(1:3),'loo')
    matrix=zeros(N,N);
    
    %This is from Adam's original code, used only for testing purposes;
    %does not work when nu is nonzero
    for j=1:N
        for k=1:N
            matrix(j,k)=(1./2)...
                .*(abs(t(j)).^(2*lambda)+abs(t(k)).^(2*lambda)-abs(t(k)-t(j)).^(2*lambda));
        end
    end
    matrix=matrix.*(-gamma(2-2*lambda).*cos(pi*lambda))./(pi.*lambda.*(2*lambda-1));
else
    %My translation
    fact=frac(1,2).*frac(-gamma(2-2*real(lambda)).*cos(pi.*lambda),pi.*lambda.*(2.*lambda-1));
    tk2H=vrep(realpow(t,2*lambda),N,2);
    tk=vrep(t,size(t,1),2);
    dtk2H=realpow(abs(tk-tk'),2*lambda);
    matrix=fact*(tk2H+tk2H'-dtk2H);
end

matrix=matrix*sigma.^2;

function[]=maternchol_test


N = 100; 
alpha=1.6;
lambda=0.15;

t1=tic;C1 = fbm_matrix(1,N,1,alpha,'chol');etime1=toc(t1);
t2=tic;C2 = fbm_matrix(1,N,1,alpha,'fast');etime2=toc(t2);

reporttest('MATERNCHOL loop and loopless algorithms match for fBM', aresame(C1./C2,ones(size(C1)),1e-10))
%disp(['Loopless algorithm was ' num2str(etime1./etime2) ' times faster than loop.'])

%alpha=0.505+rand(1,N/2)/1.1;%min(alpha),max(alpha)
alpha=1+rand(1,N/2)/1.1;%min(alpha),max(alpha)
lambda=rand(1,N/2);

alpha(alpha<=0.5)=0.5+1e-5;
alpha(alpha>=1.5)=1.5-1e-5;

[C1,C2]=vzeros(N,N,length(alpha));
t1=tic;
for i=1:length(alpha)
    C1(:,:,i) =  fbm_matrix(1,N,1,alpha(i),'chol');
end
etime1=toc(t1);
t2=tic;
for i=1:length(alpha)
    C2(:,:,i) = fbm_matrix(1,N,1,alpha(i),'fast');
end
etime2=toc(t2);

reporttest('MATERNCHOL loop and loopless algorithms match for fBm with array input', aresame(C1./C2,ones(size(C1)),1e-10))
%disp(['Loopless algorithm was ' num2str(etime1./etime2) ' times faster than loop.'])

%T=maternchol(1000,1.5,1/100);

N = 500;
alpha1=[1:.1:3.5];
h1=logspace(-1,0,20)/5;
[alpha,h]=meshgrid(alpha1,h1);

[G,G1,G2]=vzeros(length(h1),length(alpha1),N);
T=zeros(length(h1),length(alpha1),N,N);
for j=1:length(alpha1)
    for i=1:length(h1)
        Ttemp=maternchol(1,N,1,alpha1(j),h1(i),0,0);
        G(i,j,:)=flipud(Ttemp(end,:)');
        [t,G1(i,j,:)]=maternimp(N,alpha1(j),h1(i));  
    end
end

err=frac(sum(squared(G-G1),3),sum(squared(G),3));
reporttest('MATERCHOL matches form of impulse response function to within 1%', allall(err<0.01))
%See MATERNOISE for computation of error term.
