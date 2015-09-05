function[varargout]=maternoise(varargin)
%MATERNOISE  Realizations of the Matern random process and variations.  [with A. Sykulski]
%
%   Z=MATERNOISE(DT,N,SIGMA,ALPHA,LAMBDA) simulates a length N complex-
%   valued Matern random process Z having variance SIGMA^2, slope parameter
%   ALPHA, and damping parameter LAMBDA.
%
%   DT is the sample interval.  Note that LAMBDA is understood to have the
%   same units as the inverse sample interval 1/DT.  For values of LAMBDA
%   greater than about 1/DT, the process becomes essentially white noise. 
%   Thus typical values will have LAMBDA<(1/DT).
%
%   Z=MATERNOISE(DT,[N,M],SIGMA,ALPHA,LAMBDA) generate M realizations of  
%   this process, so that Z is of size N x M.
%
%   See MATERNSPEC for details on the input arguments. 
%
%   Note that the process Z is explicitly set to have zero temporal mean. 
%
%   Z=MATERNOISE(DT,N,SIGMA,ALPHA,LAMBDA,NU,MU) simulates extensions of 
%   the Matern process.  MATERNOISE(...,'composite') also works.  See 
%   MATERNSPEC for details on these options.
%
%   For further details, see Sykulski, Olhede, Lilly, and Danioux (2015),
%   "Lagrangian time series models for ocean surface drifter trajectories."
%   __________________________________________________________________
%
%   Multiple parameter values
%
%   The input arguments SIGMA, ALPHA, etc. may all be either scalars, or
%   arrays of the same length, say K.  
%
%   Z=MATERNOISE(DT,N,SIGMA,ALPHA,...), with at least one of SIGMA, ALPHA,
%   etc. being an array of length K, returns a matrix Z that is N x K. 
%
%   Z=MATERNOISE(DT,[N M],SIGMA,ALPHA,...) then gives Z of size N x M x K.
%   __________________________________________________________________
%
%   Algorithm
%
%   MATERNOISE uses a Cholesky matrix decomposition method which makes the
%   autocovariance matrix of the generated process Z have exactly the form 
%   of a sampled Matern autocovariance function, for nonzero LAMBDA, or of 
%   fractional Brownian motion for H=0. 
%
%   Note that since the Cholesky decomposition requires O(N^3) operations,
%   generating very long time series (>2000 points or so) may be slow.
%   __________________________________________________________________
%
%   See also MATERNSPEC, MATERNCOV, MATERNIMP.
%
%   'maternoise --t' runs some tests.
%   'maternoise --f' generates a sample figure.
%
%   Usage: z=maternoise(dt,N,sigma,alpha,lambda);
%          z=maternoise(dt,[N,M],sigma,alpha,lambda,nu,mu);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2015  A.M. Sykulski and J.M. Lilly
%                                 --- type 'help jlab_license' for details
 
%   THe fast algorithm is not currently working.
%   Fast algorithm 
%
%   MATERNOISE(...,'fast') uses a fast, highly accurate approximate 
%   generation algorithm.  This method works by making use of the known
%   analytic form for the Matern impulse response function.
%
%   In the fast algorithm, oversampling is used to ensure that the 
%   structure of the Green's function is accurately resolved.  Note that
%   the fast algorithm is only accurate for ALPHA>=1.
%
%   [Z,ERR]=MATERNOISE(...,'fast') also returns the error ERR involved in 
%   the fast algorithm's approximation of the autocovariance sequence.  
%  
%   See the sample figures for a test of the fast algorithm.
%   __________________________________________________________________
%

if strcmpi(varargin{1}, '--t')
    maternoise_test,return
end
if strcmpi(varargin{1}, '--f')
    type makefigs_maternoise 
    makefigs_maternoise;
    return
end


%Sort out if string is input
alg='chol';
str='noerror';
model='forward';
for i=1:3
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'fas')||strcmpi(varargin{end}(1:3),'cho')
            alg=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'for')||strcmpi(varargin{end}(1:3),'rev')||strcmpi(varargin{end}(1:3),'com')
            model=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'err')||strcmpi(varargin{end}(1:3),'noe')
            str=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end

dt=double(varargin{1});
N=varargin{2};
sigma=varargin{3};
alpha=varargin{4};
lambda=varargin{5};

nu=0;
mu=0;
if length(varargin)>5
    nu=varargin{6};
end
if length(varargin)>6
    mu=varargin{7};
end

%dt, N, sigma, alpha, nu

%if strcmpi(model(1:3),'com')||strcmpi(model(1:3),'ext')
%H=real(varargin{4})*dt;
%nu=-imag(varargin{4})*dt;
%   G=varargin{5};
%    if strcmpi(model(1:3),'com')
%        G=G.*dt;
%%    elseif strcmpi(model(1:3),'ext')
%%        G=G/dt;
%    end
%else
% G=zeros(size(sigmama));
%  H=real(varargin{4});
%   nu=-imag(varargin{4});
%end

%H=real(varargin{4})*dt;
%nu=-imag(varargin{4})*dt;
%sigma=sigma.*dt.^(alpha-1/2);

if length(N)>1
    M=N(2);
    N=N(1);
else
    M=1;
end

arrayify(sigma,alpha,lambda,mu,nu);

if anyany(lambda==0)&&strcmpi(alg(1:3),'fas')
    disp('Sorry, MATERNOISE cannot use fast algorithm with fBm case of H=0.')
    disp('Reverting to Cholesky algorithm.')
    alg='chol';
end
% if anyany(alpha<1)&&strcmpi(alg(1:3),'fas')
%     disp('Sorry, MATERNOISE cannot use fast algorithm for ALPHA<1.')
%     disp('Reverting to Cholesky algorithm.')
%     alg='chol';
% end
%if lambda(1)==0&&(anyany(alpha<1/2)||anyany(alpha>1.5))
%    error('Sorry, alpha must be between 1/2 and 3/2 for fBm case of H=0.')
%end

if strcmpi(model(1:3),'rev')&&strcmpi(alg(1:3),'fas')
    error('Sorry, MATERNOISE cannot accept both the FAST and REVERSE options.')
end

if strcmpi(model(1:3),'rev')&&anyany(lambda==0)
    error('Sorry, MATERNOISE employ REVERSE option with H set to zero.')
end



z=zeros(N,M,length(sigma));
err=nan*zeros(size(sigma));
if strcmpi(alg(1:3),'fas')
    for i=1:length(sigma)
        [z(:,:,i),err(i)]=sim_noise_fast(dt,N,M,sigma(i),alpha(i),lambda(i),nu(i),str);
    end
% elseif strcmpi(alg(1:3),'cir')
%     for i=1:length(sigma)
%         [z(:,:,i),err(i)]=sim_noise_circulant(N,M,A(i),alpha(i),nu(i),H(i),str);
%     end
else
    [g,fact]=vempty;
    for i=1:length(sigma)
        z(:,:,i)=sim_noise(dt,N,M,sigma(i),alpha(i),lambda(i),nu(i),mu(i),model);
    end
end

varargout{1}=squeeze(z);
varargout{2}=err;

function [z]=sim_noise(dt,N,M,sigma,alpha,lambda,nu,mu,model)
T=maternchol(dt,N,sigma,alpha,lambda,nu,mu,model);
Z=frac(1,sqrt(2))*(randn(N,M)+sqrt(-1)*randn(N,M));
z=T*Z;

function [z,err]=sim_noise_fast(dt,N,M,sigma,alpha,lambda,nu,str)
epsilon=0.01;

lambda=lambda*dt;
nu=nu*dt;

Ne=ceil(maternedge(alpha,lambda,epsilon));

%Error increases as N*lambda increases
%Make sure N > 20 / h ... otherwise, super-resolve
fact=max(ceil(frac(lambda*N,20)),1);
N2=(N+Ne).*fact;

[t,g]=maternimp(N2,alpha,(lambda-1i*nu)./fact);
%sqrt(frac(materncfun(alpha),lambda.^(2*alpha-1)))
%frac(1,fact).^(alpha-1)
%g=g*(1,fact).^(alpha-1);
Z=frac(sqrt(N2./fact),sqrt(2))*(randn(N2,M)+sqrt(-1)*randn(N2,M));
G=sigma*vrep(fft(g),M,2);
z=ifft(G.*Z);
z=z(1:fact:end,:);
z=z(Ne+1:end,:);

err=nan;
if strcmpi(str(1:3),'err')
    R1=conv(g,flipud(conj(g))).*(fact).^(2*(alpha-1));
    ii=(length(R1)-1)/2+1;
    R1=R1(ii:fact:ii+N*fact-1,:);
    [tau,R]=materncov(dt,N,1,alpha,lambda/dt,nu/dt);
    %    vsize(R,R1)
    %figure,plot(tau,abs([R R1]))
    err=frac(sum(squared(R1-R),1),sum(squared(R),1));
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A few routines I'm experimenting with here
function [z,err]=sim_noise_circulant(N,M,A,alpha,nu,H,str)


N2 = 2*N-2; % Size of N2xN2 circulant matrix

[tau,R]=materncov(N,A,alpha,H-1i*nu);
cev = [R;R(N-1:-1:2)];   % Top row of circulant matrix, as column vector
Sk0=real(fft(cev));      % Take FFT to get eigenvalues,
Sk=max(0,Sk0)./N2;       % Set negatives as 0, normalize

Z=randn(N2,M)+1i*randn(N2,M);   % Generate random Gaussians mean 0, variance 1
Gk=Z.*vrep(sqrt(Sk),M,2);       % Multiply by square root of eigenvalue vector 
z=fft(Gk);   % FFT back into the time domain
z=z(1:N,:);  % Take first N values of real and imaginary components

err=nan;

function [z,g,fact]=sim_noise_fast_asymmetric(N,M,A,alpha,nu,H)

epsilon=0.01;

Ne=ceil(maternedge(alpha,H,epsilon));

%Error increases as N*H increases
%Make sure N > 20 / h ... otherwise, super-resolve
fact=max(ceil(frac(H*N,20)),1);  
N2=(N+Ne).*fact;

[t,g]=maternimp(N2,alpha,(H-1i*nu)./fact);
g=g.*frac(1,fact).^(alpha-1);
g=conj(flipud(g)');

%rho= 1.4464;
%rho=1.32;
%rho=1.1071;
%rho=pi/4;
rho=1;
[z,zeps]=vzeros(N2,M);
eps1=frac(1,sqrt(fact))*randn(N2,M);
eps2=frac(1,sqrt(fact))*randn(N2,M);

om(1)=randn(1);
om=zeros(N2,1);
for i=2:N2,om(i)=0.9*om(i-1)+randn(1);end
om=om./std(om)*H/100;
%figure,plot(om)
t=[1:N2];
om=sin(t*2*pi/N2)*H/50;

for i=1:N2 
    if i==1
        zeps(1,:)=eps1(1,:)+1i*eps2(1,:);
    else
        zhat=z(i-1,:)./abs(z(i-1,:));
        %zzeps(i,:)=(eps1(i,:).*cos(rho)-0.2*frac(1,sqrt(fact))+1i*eps2(i,:).*sin(rho)).*zhat;
        zeps(i,:)=(eps1(i,:).*cos(rho)-1i*eps2(i,:).*sin(rho)).*zhat;
    end
    %z(i,:)=A*g([end-i+1:end])*zeps(1:i,:);
    gtilde=g([end-i+1:end]).*rot(om(i).*t(end-i+1:end));
    z(i,:)=A*gtilde*zeps(1:i,:);
end

z=z(1:fact:end,:);
z=z(Ne+1:end,:);

zeps=zeps(1:fact:end,:);
g=zeps(Ne+1:end,:);



function [z,g,fact]=sim_noise_fast_asymmetric_former(N,M,A,alpha,nu,H)

epsilon=0.01;

Ne=ceil(maternedge(alpha,H,epsilon));

%Error increases as N*H increases
%Make sure N > 20 / h ... otherwise, super-resolve
fact=max(ceil(frac(H*N,20)),1);  
N2=(N+Ne).*fact;

[t,g]=maternimp(N2,alpha,(H-1i*nu)./fact);
g=g.*frac(1,fact).^(alpha-1);
g=conj(flipud(g)');

rho=1.32;
%rho=1.1071;
%rho=pi/4;
[u,v,ueps,veps]=vzeros(N2,M);
eps1=frac(1,sqrt(fact))*randn(N2,M);
eps2=frac(1,sqrt(fact))*randn(N2,M);

 
for i=1:N2 
    if i==1
        ueps(1,:)=eps1(1,:);
        veps(1,:)=eps2(1,:);
    else
        uhat=u(i-1,:)./sqrt(u(i-1,:).^2+v(i-1,:).^2);
        vhat=v(i-1,:)./sqrt(u(i-1,:).^2+v(i-1,:).^2);
        ueps(i,:)=eps1(i,:).*cos(rho).*uhat-eps2(i,:).*sin(rho).*vhat;
        veps(i,:)=eps1(i,:).*cos(rho).*vhat+eps2(i,:).*sin(rho).*uhat;      
    end
    u(i,:)=A*g([end-i+1:end])*ueps(1:i,:);
    v(i,:)=A*g([end-i+1:end])*veps(1:i,:);
end

% g=flipud(vshift(g,-N+1,1));
% g=vrep(conj(g'),N,1);
% for i=1:N
%     g(i,:)=vshift(g(i,:),1-i,2);
% end


u=u(1:fact:end,:);
v=v(1:fact:end,:);

z=u+1i*v;
z=z(Ne+1:end,:);

ueps=ueps(1:fact:end,:);
veps=veps(1:fact:end,:);

g=ueps+1i*veps;
g=g(Ne+1:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[]=maternoise_test

rng(1);
N=100;
alpha=[1 1.5 2 3 4];
h=[0.001 .01 .02 .05 .2];
[alpha,h]=meshgrid(alpha,h);
h=h.*alpha;
sigma=1+5*rand(size(h));
vcolon(sigma,h,alpha);

[z,z2]=vzeros(N,1000,length(h));
for i=1:length(alpha)
    z(:,:,i)=maternoise(1,[N,1000],sigma(i),alpha(i),h(i),'chol');
    z2(:,:,i)=maternoise(1,[N,1000],sigma(i),alpha(i),h(i),'fast');
end
z=reshape(z,N*1000,length(alpha));
z2=reshape(z2,N*1000,length(alpha));
stdz1=vstd(z,1)';
stdz2=vstd(z2,1)';
eps=abs(stdz1./sigma-1);
eps2=abs(stdz2./sigma-1);

reporttest('MATERNOISE predicted and realized standard deviations match to within 2.5%, unit sample rate, Cholesky',maxmax(eps<0.025))
reporttest('MATERNOISE predicted and realized standard deviations match to within 2.5%, unit sample rate, fast',maxmax(eps2<0.025))

dt=3600;
[z,z2]=vzeros(N,1000,length(h));
for i=1:length(alpha)
    z(:,:,i)=maternoise(dt,[N,1000],sigma(i),alpha(i),h(i)./dt,'chol');
    z2(:,:,i)=maternoise(dt,[N,1000],sigma(i),alpha(i),h(i)./dt,'fast');
end
z=reshape(z,N*1000,length(alpha));
z2=reshape(z2,N*1000,length(alpha));
stdz1=vstd(z,1)';
stdz2=vstd(z2,1)';
eps=abs(stdz1./sigma-1);
eps2=abs(stdz2./sigma-1);

reporttest('MATERNOISE predicted and realized standard deviations match to within 2.5%, non-unit sample rate, Cholesky',maxmax(eps<0.025))
reporttest('MATERNOISE predicted and realized standard deviations match to within 2.5%, non-unit sample rate, fast',maxmax(eps2<0.025))


N=1000;
sigma=16;
alpha=1.4;
[psi,psilambda]=sleptap(N,4);


clear rat
n=0;
rng(1);
for i=[-3:1/2:0];
    z1=maternoise(1,N,sigma,alpha,(10.^i));
    %Note when comparing spectra, it's important to  use adaptive 
    [fhat,Spphat,Snnhat]=mspec(z1-mean(z1),psi,psilambda,'adaptive');
    %This is how to compute the blurred spectrum
    %[tau,R]=materncov(N,A,alpha,(10.^i));
    %[f,Spp,Snn]=blurspec(R); 
    [f,Spp,Snn]=maternspec(1,N,sigma,alpha,(10.^i));
    n=n+1;
    rat(n)=median(Spphat./Spp);
    %figure,plot(f,Spp),hold on,plot(fhat,Spphat),xlog,ylog
end
reporttest('MATERNOISE spectrum matches MATERNSPEC for unit sample rate',allall(abs(rat-1)<0.2))

dt=3600;
clear rat
n=0;
rng(1);
for i=[-3:1/2:0];
    z1=maternoise(dt,N,sigma,alpha,(10.^i)/dt);
    [fhat,Spphat,Snnhat]=mspec(dt,z1-mean(z1),psi,psilambda,'adaptive'); 
    [f,Spp,Snn]=maternspec(dt,N,sigma,alpha,(10.^i)/dt);
    n=n+1;
    rat(n)=median(Spphat./Spp);
    %figure,plot(f,Spp),hold on,plot(fhat,Spphat),xlog,ylog
end
reporttest('MATERNOISE spectrum matches MATERNSPEC for non-unit sample rate',allall(abs(rat-1)<0.2))


N=1000;
sigma=16;
alpha=1.4;
[psi,psilambda]=sleptap(N,4);

clear rat
n=0;
rng(1);
for i=[-3:1/2:0];
    z1=maternoise(1,N,sigma,-1/2,(10.^i),0,alpha);
    [fhat,Spphat,Snnhat]=mspec(z1-mean(z1),psi,psilambda,'adaptive');
    [f,Spp,Snn]=maternspec(1,N,sigma,-1/2,(10.^i),0,alpha);
    n=n+1;
    rat(n)=median(Spphat./Spp);
    %figure,plot(f,Spp),hold on,plot(fhat,Spphat),xlog,ylog
end
reporttest('MATERNOISE spectrum matches MATERNSPEC for exponential spectrum and unit sample rate',allall(abs(rat-1)<0.2))


clear rat
n=0;
rng(1);
dt=3600;
for i=[-3:1/2:0];
    z1=maternoise(dt,N,sigma,-1/2,(10.^i)/dt,0,alpha*dt);
    [fhat,Spphat,Snnhat]=mspec(dt,z1-mean(z1),psi,psilambda,'adaptive'); 
    [f,Spp,Snn]=maternspec(dt,N,sigma,-1/2,(10.^i)/dt,0,alpha*dt);
    n=n+1;
    rat(n)=median(Spphat./Spp);
    %figure,plot(f,Spp),hold on,plot(fhat,Spphat),xlog,ylog
end
reporttest('MATERNOISE spectrum matches MATERNSPEC for exponential spectrum and non-unit sample rate',allall(abs(rat-1)<0.2))


N=1000;alpha=1.5;H=1/100;nu=2*H;
[z,err]=maternoise(1,N,10,alpha,H,nu,'fast','err');
reporttest('MATERNOISE error between approximated and true autocovariance < 1e-7 for ALPHA=1.5',err<1e-7)

%alpha=[0.6 0.7 0.8 0.9 1 1.1 1.2 1.5 2 3 4 5 6 7 8];
alpha=[1 1.1 1.2 1.5 2 3 4 5 6 7 8];
%alpha=[0.9 1 1.1 1.2 1.5 2 3 4];
N=1000;H=1/100;nu=2*H;
[z,err]=maternoise(1,N,10,alpha,H,nu,'fast','err');
reporttest('MATERNOISE error between approximated and true autocovariance < 1e-6 ALPHAs >=1',err<1e-6)

%N=1000;H=1/10;nu=2*H;
%[z,err]=maternoise(1,N,10,alpha,H,nu,'fast','err');
%reporttest('MATERNOISE error between approximated and true autocovariance < 1e-6 ALPHAs >=1, larger H',err<1e-6)


function[]=maternoise_figure1
 

% This goes after the script that is now in MAKEFIGS_MATERNOISE
% Find a better way to show spectra.
% [psi,lambda]=sleptap(N,8);
% [f,spp,snn,spn]=mspec(z,psi,lambda,'adaptive');  
% 
% [fo,sppo,snno]=maternspec(N,10./sigma,alpha,h);
% 
% logfo=log(fo);
% logfo=logfo./(maxmax(logfo)-minmin(logfo));
% 
% zs=vrep(logfo,size(spp,2),2)+sqrt(-1)*unitmax(log(spp));
% zso=vrep(logfo,size(spp,2),2)+sqrt(-1)*unitmax(log(sppo));
% 
% xos=xo(1:size(zs,1),:);
% 
% figure,
% plot(zso+xos*3),axis equal,xlabel('Increasing Delta'),ylabel('Increasing h')
% title('Spectra of Matern Processes')


function[]=maternoise_figure2
 
N=1000;

alpha=[1 1.5 2 3 4];
h=[.01 .02 .05 .2 1];
[alpha,h]=meshgrid(alpha,h);
h=h.*alpha;

d=frac(h.^(2*alpha-1),materncfun(alpha));

rng(1);  %set seed
tic;z=maternoise(1,N,10./sqrt(d),alpha,h);etime1=toc;
tic;z2=maternoise(1,N,10./sqrt(d),alpha,h,'fast');etime2=toc;

[psi,lambda]=sleptap(N,8);
[om,spp,snn]=mspec(z,psi,lambda,'adaptive');    
[om,spp2,snn2]=mspec(z2,psi,lambda,'adaptive');    

[f,sppo,snno]=maternspec(1,N,10./sqrt(d),alpha,h);

disp(['MATERNOISE fast algorithm was ' num2str(etime1./etime2) ' times faster than Cholesky algorithm.'])

figure
subplot(1,2,1),plot(om,spp),hold on,plot(om,sppo),xlog,ylog,title('Matern process, Cholesky algorithm')
axis tight,ylim(10.^[-4 14])
subplot(1,2,2),plot(om,spp2),hold on,plot(om,sppo),xlog,ylog,title('Matern process, Fast algorithm')
axis tight,ylim(10.^[-4 14])

function[]=maternoise_figure3
 
N=1000;

alpha=[1 1.5 2 3 4];
h=[.01 .02 .05 .2 1];

[alpha,h]=meshgrid(alpha,h);
h=h.*alpha;
nu=0.05.*alpha;

d=frac(h.^(2*alpha-1),materncfun(alpha));

rng(1);  %set seed
tic;z=maternoise(1,N,10./sqrt(d),alpha,h,nu);etime1=toc;
tic;z2=maternoise(1,N,10./sqrt(d),alpha,h,nu,'fast');etime2=toc;

[psi,lambda]=sleptap(N,8);
[om,spp,snn]=mspec(z,psi,lambda,'adaptive');    
[om,spp2,snn2]=mspec(z2,psi,lambda,'adaptive');    

[f,sppo,snno]=maternspec(1,N,10./sqrt(d),alpha,h,nu);

disp(['MATERNOISE fast algorithm was ' num2str(etime1./etime2) ' times faster than Cholesky algorithm.'])

figure
subplot(1,2,1),plot(om,spp),hold on,plot(om,sppo),xlog,ylog,title('Oscillatory Matern, Cholesky algorithm')
axis tight,ylim(10.^[-4 14])
subplot(1,2,2),plot(om,spp2),hold on,plot(om,sppo),xlog,ylog,title('Oscillatory Matern, Fast algorithm')
axis tight,ylim(10.^[-4 14])


function[]=maternoise_figure4

N=1000;
A=16;
alpha=1.4;
h=.1;

rng(1);
[psi,lambda]=sleptap(N,4);
z1=maternoise(1,N,A,alpha,h,'fast');
z2=maternoise(1,N,A,alpha,h);
[fhat,Spphat,Snnhat]=mspec(z1,psi);
[fhat2,Spphat2,Snnhat2]=mspec(z2,psi);
[f,Spp,Snn]=maternspec(1,N,A,alpha,h);
figure
subplot(1,2,1),plot(f,Spp),hold on,plot(fhat,Spphat),plot(fhat2,Spphat2),xlog,ylog
title('Noise spectra and true spectrum, DT=1'),axis tight

dt=3600;
rng(1);
[psi,lambda]=sleptap(N,4);
z1=maternoise(dt,N,A,alpha,h./dt,'fast');
z2=maternoise(dt,N,A,alpha,h./dt);
[fhat,Spphat,Snnhat]=mspec(dt,z1,psi);
[fhat2,Spphat2,Snnhat2]=mspec(dt,z2,psi);
[f,Spp,Snn]=maternspec(dt,N,A,alpha,h./dt);
subplot(1,2,2),
plot(f,Spp),hold on,plot(fhat,Spphat),plot(fhat2,Spphat2),xlog,ylog
title('Noise spectra and true spectrum, DT=3600'),axis tight


% N=1000;
% sigma=1;
% alpha=1.5;
% h=0.1;
% z=maternoise(N,sigma,alpha,h,h/5-1i*h*5,'composite');

