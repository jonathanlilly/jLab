function[varargout]=materncov(varargin)
%MATERNCOV  Autocovariance of the Matern random process and variations.
%
%   [TAU,R]=MATERNCOV(DT,N,SIGMA,ALPHA,LAMBDA) returns the autocovariance 
%   function R of a length N complex-valued Matern random process having 
%   variance SIGMA^2, slope parameter ALPHA, and damping parameter LAMBDA.
%
%   DT is the sample interval.  Note that LAMBDA is understood to have the
%   same units as the inverse sample interval 1/DT.
%
%   TAU is an array of time lags at which R is computed, and is given by 
%   TAU=DT*[0,1,...,N-1].
%
%   The Matern autocovariance is given by
%
%      R(TAU) = SIGMA^2 * (LAMBDA|TAU|)^(ALPHA-1/2) 
%                             * BESSELK(ALPHA-1/2,LAMBDA|TAU|) * C
%
%   where C is a normalizing constant dependent upon ALPHA and LAMBDA.
%
%   By definition, R is one-sided theoretical autocovariance at 
%   non-negative time lags.  See below for the relationship between this 
%   and the full, length (2N-1) theoretical autocovariance function. 
%
%   Note that for LAMBDA=0, the case of fractional Brownian motion, R will 
%   contain only INFs because the autocovariance function is unbounded.
%
%   The input parameters SIGMA, ALPHA, and LAMBDA, may all either be scalars 
%   of arrays of the same length M.  If the latter, then the output 
%   autocovariance function R will be a matrix with N rows and M columns. 
%
%   [TAU,R]=MATERNCOV(DT,N,SIGMA,ALPHA,LAMBDA,NU,MU) returns the 
%   autocovariance function of various extensions of the Matern process. 
%   MATERNOISE(...,'composite') also works.  See MATERNSPEC for details.
%
%   See MATERNSPEC for a more thorough discussion of the Matern process.
%
%   For further details, see Sykulski, Olhede, Lilly, and Danioux (2015),
%   "Lagrangian time series models for ocean surface drifter trajectories."
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
%   Composite Matern
%
%   [TAU,R]=MATERNCOV(DT,N,SIGMA,ALPHA,LAMBDA,NU,MU,'composite'), returns 
%   the autocovariance of the composite Matern process.  See MATERNSPEC for
%   more details.
%
%   This autocovariance does not have a simple analytic form, but is 
%   computed in approximate form from inverse Fourier tranforming the
%   spectrum with 10x oversampling and a length 100 x N. 
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
%          [tau,R]=materncov(dt,N,sigma,alpha,lambda,nu,mu);
%          [tau,R]=materncov(dt,N,sigma,alpha,lambda,nu,mu,'composite');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2015  J.M. Lilly --- type 'help jlab_license' for details

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


if strcmpi(varargin{1}, '--t')
    materncov_test,return
end

sid='one';
model='forward';

for i=1:3
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'rev')||strcmpi(varargin{end}(1:3),'com')
            model=varargin{end}(1:3);
        else
            sid=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end

if strcmpi(model(1:3),'com')&&strcmpi(model(1:3),'rev')
    error('Sorry, the COMPOSITE and REVERSE options cannot be used together.')
end

dt=double(varargin{1});
N=floor(varargin{2});
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

%dt,N,sigma,alpha,lambda,nu,mu

arrayify(sigma,alpha,lambda,nu,mu);
tau=dt*[0:N-1]';

R=zeros(N,length(lambda));

for i=1:length(lambda)
    if ~sigma(i)==0
        R(:,i)=materncov_one(tau,sigma(i),alpha(i),lambda(i),nu(i),mu(i),model);
    end
end

if strcmpi(sid(1:3),'ful')
    tau=[-flipud(tau(2:end,:));tau];
    R=[flipud(conj(R(2:end,:)));R];
end
varargout{1}=tau;
varargout{2}=R;


function[R]=materncov_one(tau,sigma,alpha,lambda,nu,mu,model)

%sigma,alpha,lambda,mu,model

if strcmpi(model(1:3),'com')
    M=10;
    P=10;
    N=length(tau);
    dt=tau(2)-tau(1);
    [f,Spp,Snn]=maternspec(dt,M*N*P,sigma,alpha,lambda/P,nu/P,mu/P,'composite');
    %[f,Spp,Snn]=maternspec(dt,M*N,sigma,alpha,lambda,nu,mu,'composite');
    S=[flipud(Snn(2:end));Spp];
    Ri=ifft(ifftshift(S))./dt;  %Make sure it's ifftshift not fftshift
    Ri=Ri(1:P:end);
    R=Ri(1:N);
else
    if mu==0
        fact=2*frac(1,gamma(alpha-1/2).*pow2(alpha-1/2));
        R=fact.*((lambda*tau).^(alpha-1/2)).*besselk(abs(alpha-1/2),lambda*tau);
    else
        if alpha==-1/2
            tnorm=sqrt(tau.^2+mu.^2);
            fact=besselk(1,mu.*lambda);
            R=frac(1,fact).*frac(mu,tnorm).*besselk(1,lambda.*tnorm);
        else
            tnorm=lambda.*sqrt(tau.^2+mu.^2);
            fact=1./(abs(mu.*lambda)).^(alpha-1/2)./besselk(abs(alpha-1/2),abs(mu.*lambda));
            %/maternfun(alpha(i),lambda(i).*sqrt(G(i).^2));
            R=fact*tnorm.^(alpha-1/2).*besselk(abs(alpha-1/2),tnorm);
        end
    end
    if nu~=0
        R=R.*exp(sqrt(-1)*tau*nu);
    end
    R(tau==0)=1;
    R=R.*sigma.^2;
end


%             M=10;
%             taui=[0:M*N-1]'./M;
%             taufull=[-flipud(taui(2:end,:));taui];
%             R1=maternfun(alpha(i)-1/2,G(i).*abs(taufull));
%             R2=maternfun(1/2,lambda(i).*abs(taufull)).*exp(sqrt(-1)*taufull*nu(i));
%             %figure,uvplot(R2),hold on,plot(R1)
%             Rfull=conv(R1,R2);
%             a=(length(Rfull)+1)/2-(length(taufull)-1)/2;
%             b=(length(Rfull)+1)/2+(length(taufull)-1)/2;
%             Rfull=Rfull(a:b);
%             Rfull=Rfull((end+1)/2:end);
%             Rfull=Rfull(1:M:end)./Rfull(1);
%             R(:,i)=sigma(i).^2.*Rfull;
% 
%         elseif strcmpi(model(1:3),'rev')
%             %disp('MATERNCOV computing reverse autocovariance function.')
%             R(:,i)=sigma(i).^2.*frac(1,1+squared(lambda(i).*abs(tau))).^alpha(i);
%        

function[]=materncov_test
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
dt=7;
[tau,R]=materncov(dt,N,sigma,alpha,h/dt);
[f,Spp,Snn]=maternspec(dt,N,sigma,alpha,h/dt);
[f,Spp2,Snn2]=blurspec(dt,R,'aliased');

b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-4);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-4);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC ..., even N, non-unit sample rate',b1&&b2)

N=1000;
[f,Spp,Snn]=maternspec(2,N,sigma,alpha,h/5,h*5,h,'composite');
[tau,R]=materncov(2,N,sigma,alpha,h/5,h*5,h,'composite');
[f,Spp2,Snn2]=blurspec(2,R,'aliased');

%[tau,R]=materncov(2,N,sigma,alpha,h/5,h*5,h,'composite');
%[f,Spp3,Snn3]=blurspec(2,R,'aliased');

%z=maternoise(N,sigma,alpha,h/5-1i*h*5,h,'composite');

%fact=2*sigma.^2.*(h/5).*((h*5).^2+h.^2).^alpha;
%2*sigma.^2.*H.*(nu.^2+G.^2).^alpha;

b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),2e-2);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),2e-2);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC for case composite Matern with negligible aliasing, even N',b1&&b2)

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

N=1000;
[tau,R]=materncov(1,N,sigma,-1/2,h,0,alpha*10);
[f,Spp,Snn]=maternspec(1,N,sigma,-1/2,h,0,alpha*10);
[f,Spp2,Snn2]=blurspec(1,R,'aliased');


b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-4);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-4);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC for case with negligible aliasing, even N, exponential',b1&&b2)


N=1000;
dt=7;
[tau,R]=materncov(dt,N,sigma,-1/2,h/dt,0,alpha*10*dt);
[f,Spp,Snn]=maternspec(dt,N,sigma,-1/2,h/dt,0,alpha*10*dt);
[f,Spp2,Snn2]=blurspec(dt,R,'aliased');


b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-4);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-4);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC ..., even N, exponential, non-unit sample rate',b1&&b2)



[tau,R]=materncov(1,1000,7,1.5,1/10,0,20);
[f,Spp,Snn]=maternspec(1,1000,7,1.5,1/10,0,20);
[f,Spp2,Snn2]=blurspec(1,R,'aliased');


b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-10);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-10);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC for case with negligible aliasing, even N, extended Matern',b1&&b2)

dt=7;
[tau,R]=materncov(dt,1000,7,1.5,1/10/dt,0,20*dt);
[f,Spp,Snn]=maternspec(dt,1000,7,1.5,1/10/dt,0,20*dt);
[f,Spp2,Snn2]=blurspec(dt,R,'aliased');


b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-10);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-10);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC ..., even N, extended Matern, non-unit sample rate',b1&&b2)


N=1000-1;
[tau,R]=materncov(1,N,sigma,alpha,h,2*h);
[f,Spp,Snn]=maternspec(1,N,sigma,alpha,h,2*h);
[f,Spp2,Snn2]=blurspec(1,R,'aliased');

b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-4);
b2=aresame(Snn2./maxmax(Spp),Snn./maxmax(Spp),1e-4);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC for case with negligible aliasing, frequency shift case, odd N',b1&&b2)


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

