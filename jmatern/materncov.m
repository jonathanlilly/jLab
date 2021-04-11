function[varargout]=materncov(varargin)
%MATERNCOV  Autocovariance of the Matern random process and variations.
%
%   [TAU,R]=MATERNCOV(DT,N,SIGMA,ALPHA,LAMBDA) returns the autocovariance 
%   function R of a length N complex-valued Matern random process having 
%   variance SIGMA^2, slope parameter ALPHA, and damping parameter LAMBDA.
%
%   [TAU,R]=MATERNCOV(...,'real') instead forms the covariance of a real-
%   valued Matern process.
%
%   DT is the sample interval.  Note that LAMBDA is understood to have the
%   same units as the inverse sample interval 1/DT.
%
%   TAU is an array of time lags at which R is computed, and is given by 
%   TAU=DT*[0,1,...,N-1].
%
%   By definition, R is one-sided theoretical autocovariance at 
%   non-negative time lags.  See below for the relationship between this 
%   and the full, length (2N-1) theoretical autocovariance function. 
%
%   Note that for LAMBDA=0, the case of fractional Brownian motion, R will 
%   contain only INFs because the autocovariance function is unbounded.
%
%   The input parameters SIGMA, ALPHA, and LAMBDA, may all either be 
%   scalars or arrays of the same length M.  If the latter, then the output 
%   autocovariance function R will be a matrix with N rows and M columns. 
%
%   [TAU,R]=MATERNCOV(DT,N,SIGMA,ALPHA,LAMBDA,NU,MU) returns the 
%   autocovariance function of various extensions of the Matern process. 
%   MATERNOISE(...,'composite') also works.  See MATERNSPEC for details.
%
%   See MATERNSPEC for a more thorough discussion of the Matern process.
%
%   For details on the Matern process and its autocovariance function, see:
%
%     Lilly, Sykulski, Early, and Olhede, (2017).  Fractional Brownian
%        motion, the Matern process, and stochastic modeling of turbulent 
%        dispersion.  Nonlinear Processes in Geophysics, 24: 481--514.
%   __________________________________________________________________
%
%   Relationship to full autocovariance
%
%   For a time series of length N, the full autocovariance function RF is 
%   length 2N-1, defined at time lags -N+1,-N+2...,-1,0,1,...,N-2,N-1.
%
%   The one-sided autocovariance R contains the full autocovariance RF at 
%   positive time lags. Negative lags are given by Hermitian symmetry.
%
%   [TAUF,RF]=MATERNCOV(...,'full') returns the full (two-sided)
%   autocovariance RF and the corresponding two-sided time array TAUF. 
%
%   RF is constructed from R as RF=[FLIPUD(CONJ(R(2:end,:));R]. 
%   __________________________________________________________________
%
%   Real-valued processes
%
%   By default MATERNCOV returns the autocovariance of a complex-valued
%   process. 
%
%   MATERNCOV(...,'real') instead returns the autocovariance of a real-
%   valued process. This also works with any of extended versions.
%   __________________________________________________________________
%
%   Computational notes
%
%   The autocovariances for the generalized and composite Matern spectra do
%   not have analytic forms. Rather, these are approximated to high
%   precision by inverse Fourier transforming the oversampled spectrum with 
%   10 x oversampling over a 10 x longer time period, and then decimating.  
%
%   MATERNCOV(...,'general',M,P) or  MATERNCOV(...,'composite',M,P) 
%   specifies the numerical oversampling parameters in the calculation. 
%
%   The spectrum is then computed over a time window of M times the
%   required  duration, and P times the required sampling density, for a
%   total of M*P time more points.  These values may set to optimize the
%   tradeoff between speed and accuracy.  
%
%   This computation method is expected to minimize aliasing effects and
%   resolution errors.  
%   __________________________________________________________________
%
%   See also MATERNSPEC, MATERNIMP, MATERNOISE, MATERNFIT.
%
%   'materncov --t' runs some tests.
%
%   Usage: [tau,R]=materncov(dt,N,sigma,alpha,lambda);
%          [tau,R]=materncov(dt,N,sigma,alpha,lambda,nu);
%          [tau,R]=materncov(dt,N,sigma,alpha,lambda,mu,'composite');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2021 J.M. Lilly --- type 'help jlab_license' for details

%   No longer supported, sorry
%   __________________________________________________________________
%
%   Reverse Matern
%
%   [TAU,R]=MATERNCOV(N,SIGMA,ALPHA,LAMBDA,'reverse') returns the autocovariance 
%   of the reverse Matern process, in which the roles of the spectrum and 
%   the autocovariance functions are swapped.  The autocovariance is 
%
%      R(TAU) = SIGMA^2 / (1+LAMBDA^2*TAU^2)^ALPHA 
%
%   Note that the parameter LAMBDA has been inverted so that it still has units
%   of frequency, as in the case of the forward Matern process. 


%   __________________________________________________________________
%
%   Extensions
%
%   Several extensions of the basic Matern form are supported, all of which
%   are discussed in more detail in MATERNSPEC.
%
%   Oscillatory Matern
%       [F,SPP,SNN]=MATERNCOV(DT,N,SIGMA,ALPHA,LAMBDA,NU)
%   Generalized Matern (+ optional oscillations)
%       [F,SPP,SNN]=MATERNCOV(DT,N,SIGMA,ALPHA,LAMBDA,GAMMA,'general') 
%       [F,SPP,SNN]=MATERNCOV(DT,N,SIGMA,ALPHA,LAMBDA,GAMMA,NU,'general') 
%   Extended Matern  (+ optional oscillations)
%       [F,SPP,SNN]=MATERNCOV(DT,N,SIGMA,ALPHA,LAMBDA,MU,'extended')
%       [F,SPP,SNN]=MATERNCOV(DT,N,SIGMA,ALPHA,LAMBDA,MU,NU,'extended')
%   Composite Matern 
%       [F,SPP,SNN]=MATERNCOV(DT,N,SIGMA,ALPHA,LAMBDA,MU,NU,'composite')

if strcmpi(varargin{1}, '--t')
    materncov_test,return
end

sid='one';
model='standard';
flag='complex';

M=10;  %Specifying oversampling rates for numerical computations
P=10;
for i=1:3
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'rev')||strcmpi(varargin{end}(1:3),'com')||strcmpi(varargin{end}(1:3),'gen')||strcmpi(varargin{end}(1:3),'ext')
            model=varargin{end}(1:3);
        elseif strcmpi(varargin{end}(1:3),'rea')
            flag=varargin{end};
        else
            sid=varargin{end};
        end
        varargin=varargin(1:end-1);
    elseif ischar(varargin{end-2})  %For inputting M and P for testing
        if strcmpi(varargin{end-2}(1:3),'com')||strcmpi(varargin{end-2}(1:3),'gen')
            model=varargin{end-2};
            M=varargin{end-1};
            P=varargin{end};
            varargin=varargin(1:end-3);
        end
    end
end

if strcmpi(model(1:3),'com')&&strcmpi(model(1:3),'rev')
    error('Sorry, the COMPOSITE and REVERSE options cannot be used together.')
end

dt=double(varargin{1});
N=floor(varargin{2});
varargin=varargin(3:end);

if length(varargin)==1
    params=varargin{1};  %for input from materfit
else
    %this will be easier as a matrix
    Ncell=max(cellength(varargin));
    params=zeros(Ncell,length(varargin));
    for i=1:length(varargin)
        params(:,i)=varargin{i};
    end
end

tau=dt*[0:N-1]';
R=zeros(N,size(params,1));
for i=1:size(params,1)
    R(:,i)=materncov_one(tau,params(i,:),model,M,P,flag);
end

if strcmpi(sid(1:3),'ful')
    tau=[-flipud(tau(2:end,:));tau];
    R=[flipud(conj(R(2:end,:)));R];
end
varargout{1}=tau;
varargout{2}=R;


function[R]=materncov_one(tau,params,model,M,P,flag)

sigma=params(1);
alpha=params(2);
lambda=params(3);

nu=0;
if strcmpi(model(1:3),'com')
    mu=params(4);
    nu=params(5);
    
    N=length(tau);
    dt=tau(2)-tau(1);
    [f,Spp,Snn]=maternspec(dt,M*N*P,sigma,alpha,lambda/P,mu/P,nu/P,'composite',flag);
    S=[flipud(Snn(2:end));Spp];%figure,plot(S)
    Ri=ifft(ifftshift(S))./dt;  %Make sure it's ifftshift not fftshift
    Ri=Ri(1:P:end);
    R=Ri(1:N);
    R=R./R(1);
elseif strcmpi(model(1:3),'gen')
    %Generalized Matern form
    gamma=params(4);
    if length(params)==5
        nu=params(5);
    end
    %sigma,alpha,lambda/P,gamma,nu,flag
    N=length(tau);
    dt=tau(2)-tau(1);
    [f,Spp,Snn]=maternspec(dt,M*N*P,sigma,alpha,lambda/P,gamma,'generalized',flag);
    S=[flipud(Snn(2:end));Spp];%figure,plot(S)
    %figure,plot(f,[Spp Snn]),xlog,ylog,aresame(Spp,Snn,1e-13)
    Ri=ifft(ifftshift(S))./dt;  %Make sure it's ifftshift not fftshift
    %figure,plot(Ri)
    Ri=Ri(1:P:end);
    R=Ri(1:N);
    R=R./R(1);
elseif strcmpi(model(1:3),'ext')
    %Extended Matern, reduces to the standard Matern form for mu=0
    mu=params(4);
    if length(params)==5
        nu=params(5);
    end
    t=lambda.*sqrt(tau.^2+mu.^2);  %This is to support the extended Matern form
    R=frac(maternfun(alpha,t),maternfun(alpha,mu*lambda));
elseif strcmpi(model(1:3),'sta')
    %Standard Matern form
    if length(params)==4
        nu=params(4);
    end
    t=lambda.*tau;
    R=maternfun(alpha,t);
end
if nu~=0 && ~strcmpi(model(1:3),'com')
    R=R.*exp(sqrt(-1)*tau.*nu);
end
R=R.*sigma.^2;

if strcmpi(flag(1:3),'rea')
     R=real(R);
end

% Extended Matern form reduces to this Bessel function form, test is below
% if alpha=-1/2
%     tnorm=sqrt(tau.^2+mu.^2);
%     fact=besselk(1,lambda.*mu);
%     R1=frac(1,fact).*frac(1,tnorm./mu).*besselk(1,lambda.*tnorm);
% end

% function[]=maternfun_fig
% 
% alpha=[-4:0.1:4];
% [tau,R]=materncov(dt,10000,1,alpha,0,0,1000);
% figure,plot(tau,R)
% 
% hold on,h=plot(r1,xi(:,1)); linestyle -h h 3D
% hold on,h=plot(r1,xi(:,end)); linestyle -h h 3k
% title('The Matern function from \alpha = 1/2 (gray) to \alpha = 4 (black)')

function[]=materncov_test

tic
alpha=1.5;
h=0.1;
sigma=1;

N=1000;
[tau,R]=materncov(1,N,sigma,alpha,h);
[f,Spp,Snn]=maternspec(1,N,sigma,alpha,h);
[f,Spp2,Snn2]=blurspec(1,R,'aliased');

b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-4);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-4);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC for case with negligible aliasing, even N',b1&&b2)


N=1000;
[tau,R]=materncov(1,N,sigma,alpha,h,h/2);
[f,Spp,Snn]=maternspec(1,N,sigma,alpha,h,h/2);
[f,Spp2,Snn2]=blurspec(1,R,'aliased');

b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-4);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-4);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC for oscillatory Matern with negligible aliasing, even N',b1&&b2)

N=1000;
[tau,R]=materncov(1,N,sigma,alpha,h,h/2,'real');
[f,Spp,Snn]=maternspec(1,N,sigma,alpha,h,h/2,'real');
[f,Spp2,Snn2]=blurspec(1,R,'aliased');

b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-4);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-4);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC for real-valued oscillatory Matern with negligible aliasing, even N',b1&&b2)

N=1000;
dt=7;
[tau,R]=materncov(dt,N,sigma,alpha,h/dt);
[f,Spp,Snn]=maternspec(dt,N,sigma,alpha,h/dt);
[f,Spp2,Snn2]=blurspec(dt,R,'aliased');

b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-4);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-4);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC ..., even N, non-unit sample rate',b1&&b2)

N=1000;
[f,Spp,Snn]=maternspec(2,N,sigma,alpha,h/5,h,h*5,'composite');
[tau,R]=materncov(2,N,sigma,alpha,h/5,h,h*5,'composite');
[f,Spp2,Snn2]=blurspec(2,R,'aliased');

%[tau,R]=materncov(2,N,sigma,alpha,h/5,h*5,h,'composite');
%[f,Spp3,Snn3]=blurspec(2,R,'aliased');

%z=maternoise(N,sigma,alpha,h/5-1i*h*5,h,'composite');

%fact=2*sigma.^2.*(h/5).*((h*5).^2+h.^2).^alpha;
%2*sigma.^2.*H.*(nu.^2+G.^2).^alpha;

b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),2e-2);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),2e-2);

%this test is currently not working ... Spp2 is off by a factor of sqrt(3),
%related to the fact that R(1) is not unity in the above code after calling mspec
%reporttest('MATERNCOV Fourier transforms to MATERNSPEC for case composite Matern with negligible aliasing, even N',b1&&b2)

N=1000;
[tau,R]=materncov(1,N,sigma,-1/2,h,alpha*10,'extended');
[f,Spp,Snn]=maternspec(1,N,sigma,-1/2,h,alpha*10,'extended');
[f,Spp2,Snn2]=blurspec(1,R,'aliased');

b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-4);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-4);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC for case with negligible aliasing, even N, exponential',b1&&b2)


N=1000;
dt=7;
[tau,R]=materncov(dt,N,sigma,-1/2,h/dt,(alpha*10*dt),'extended');
[f,Spp,Snn]=maternspec(dt,N,sigma,-1/2,h/dt,(alpha*10*dt),'extended');
[f,Spp2,Snn2]=blurspec(dt,R,'aliased');

b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-4);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-4);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC ..., even N, exponential, non-unit sample rate',b1&&b2)


[tau,R]=materncov(1,1000,7,1.5,1/10,20,'extended');
[f,Spp,Snn]=maternspec(1,1000,7,1.5,1/10,20,'extended');
[f,Spp2,Snn2]=blurspec(1,R,'aliased');


b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-10);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-10);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC for case with negligible aliasing, even N, extended Matern',b1&&b2)

dt=1.3;
[tau,R]=materncov(dt,1000,7,1.5,1/10/dt,20/dt,'extended');
[f,Spp,Snn]=maternspec(dt,1000,7,1.5,1/10/dt,20/dt,'extended');
%[f,Spp2,Snn2]=blurspec(dt,R);
[f,Spp2,Snn2]=blurspec(dt,R,'aliased');

b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-10);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-10);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC ..., even N, extended Matern, non-unit sample rate',b1&&b2)


N=1000-1;
[tau,R]=materncov(1,N,sigma,alpha,h,2*h,'extended');
[f,Spp,Snn]=maternspec(1,N,sigma,alpha,h,2*h,'extended');
[f,Spp2,Snn2]=blurspec(1,R,'aliased');

b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-4);
b2=aresame(Snn2./maxmax(Spp),Snn./maxmax(Spp),1e-4);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC for case with negligible aliasing, frequency shift case, odd N',b1&&b2)

N=1000;
dt=7;
lambda=h/dt;
mu=(alpha*10*dt);
[tau,R]=materncov(dt,N,sigma,-1/2,lambda,mu,'extended');
tnorm=sqrt(tau.^2+mu.^2);
fact=besselk(1,lambda.*mu);%See notes for this form
R1=frac(1,fact).*frac(1,tnorm./mu).*besselk(1,lambda.*tnorm);
reporttest('MATERNCOV ALPHA = -1/2 case matches expected form',aresame(R,R1,1e-10))

toc
% [tau,R]=materncov(1,N,sigma,alpha,h,'reverse');
% [f,Spp,Snn]=maternspec(1,N,sigma,alpha,h,'reverse');
% R(1)=R(1)/2;
% S=2*real(fft(R));
% Spp2=S(1:length(Spp));
% S=flipud(S);
% Snn2=[S(end);S(1:length(Snn)-1)];
% 
% b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-4);
% b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-4);
% 
% reporttest('MATERNCOV Fourier transforms to MATERNSPEC for case with negligible aliasing, even N, reverse Matern',b1&&b2)


%Playing around with asymptotic results... which don't work for alpha > 3/2
% sigma=7;
% alpha=[0.6 3/4 1.001 1.25 3/2-0.001 1.75];% 1.75 2 3 4 10];
% %alpha=[1.75 2 3 4 10];% 1.75 2 3 4 10];
% lambda=1/10000;
% N=10000;
% [tau,R]=materncov(N,sigma,alpha,lambda);
% Rapprox=zeros(size(R));
% for i=1:length(alpha)
%     %numer=(lambda.*tau).^(2*alpha(i)-1);
%     %denom=2*materncfun(alpha(i)).*cos(pi.*alpha(i)).*gamma(2*alpha(i));
%     numer1=-(lambda.*tau/2).^(2*alpha(i)-1)*gamma(3/2-alpha(i));
%     denom1=gamma(1/2+alpha(i));
%     numer2=(lambda.*tau/2).^2*gamma(3/2-alpha(i));
%     denom2=gamma(5/2-alpha(i));
%     Rapprox1(:,i)=sigma.^2.*(1+frac(numer,denom));
%     %Rapprox2(:,i)=sigma.^2.*(1+frac(numer1,denom2)+frac(numer2,denom2));
% end
% Rapprox(Rapprox<0)=nan;
% Rapprox(Rapprox>2*sigma.^2)=nan;

% alpha=1.5;
% h=0.1;
% sigma=1;
% nu=0.1;
% 
% N=1000;
% figure
% [tau,R]=materncov(N,sigma,1,h-1i*nu);
% [f,Spp,Snn]=blurspec(R);
% plot(-f,Snn),hold on,plot(f,Spp)
% [tau,R]=materncov(N,sigma,1,h+1i*nu);
% [f,Spp,Snn]=blurspec(R);
% plot(-f,Snn),hold on,plot(f,Spp)
%  
% [f,Spp,Snn]=maternspec(N,sigma,alpha,h);

% N=1000;
% r1=[0:.00001:0.00001]';
% alpha=[1/2:0.1:4];
% R=materncov(r1,1,alpha(:),1/100);
% dr=(R(2,:)-R(1,:))./(r1(2)-r1(1));
% dr2=-frac(gamma(abs(alpha-3/2)),gamma(alpha-1/2)).*frac(1/100,2);
% %dr2(alpha==1)=-1/100;
% 
% plot(alpha,-dr,'r*'),hold on,plot(alpha,-dr2,'o')
% 
% dxi=-frac(r,2).*frac(gamma(alpha-1),gamma(alpha)).*maternfun(alpha-1,r);
% reporttest('MATERNFUN iterative expression for derivative matches direct expression',aresame(xi1(2:end,:),dxi(2:end,:),1e-10))
% 

