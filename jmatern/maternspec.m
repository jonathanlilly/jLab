function[varargout]=maternspec(varargin)
%MATERNSPEC  Fourier spectrum of the Matern random process and variations.
% 
%   [F,S]=MATERNSPEC(DT,N,SIGMA,ALPHA,LAMBDA) returns the spectrum S of a  
%   length N complex-valued Matern random process having variance SIGMA^2, 
%   slope parameter ALPHA, and damping parameter LAMBDA.
%
%   DT is the sample interval.  Note that LAMBDA is understood to have the
%   same units as the inverse sample interval 1/DT.
%
%   F is an array of one-sided (positive) Fourier frequencies for a time
%   series of length N, F=FOURIER(N).  Note that F is a *radian* frequency. 
%
%   The lengths of the output variables F and S are N/2+1 for even N, and
%   (N+1)/2 for odd N.
%     
%   S is the postive rotary spectrum given by
%
%        S(F) = SIGMA^2 / (F^2+LAMBDA^2)^ALPHA * C
%
%   where C is a normalizing constant dependent upon ALPHA and LAMBDA.  The
%   negative rotary spectrum takes the same form.
%
%   For LAMBDA=0, the Matern spectrum reduces to the spectrum of fractional
%   Brownian motion.  
%
%   For further details, see Sykulski, Olhede, Lilly, and Danioux (2015),
%   "Lagrangian time series models for ocean surface drifter trajectories."
%   __________________________________________________________________
%
%   Matrix and cell array output
%
%   [F,S]=MATERNSPEC(DT,N,SIGMA,ALPHA,LAMBDA) where N is a scalar while the
%   other input arguments are all either scalars or arrays of the same 
%   length M, gives an output spectra S with LENGTH(F) rows and M columns. 
%
%   [F,S]=MATERNSPEC(DT,N,SIGMA,ALPHA,LAMBDA) where N is an array of M 
%   different lengths, returns F and S that are length M cell arrays.  Then 
%   SIGMA, ALPHA, and LAMBDA may all either be scalars or length M arrays.
%
%   This latter format is convenient for generating sets of spectra that 
%   do not all have the same size. 
%
%   When N is an array, MATERNSPEC(...,'parallel') parallelizes the 
%   computation of the various spectra using a PARFOR loop.  This option
%   requires that Matlab's Parallel Computing Toolbox be installed.
%
%   The matrix and cell array formats also work for the variations of the 
%   Matern process described below. 
%   __________________________________________________________________
%
%   Oscillatory Matern
%
%   [F,SPP,SNN]=MATERNSPEC(DT,N,SIGMA,ALPHA,LAMBDA,NU) with six input 
%   arguments modifies the spectrum to have a rotation frequency NU. 
%
%   This is accomplished by shifting the spectrum to be centered on F=NU 
%   rather than F=0.  SPP and SNN are now the postive rotary and negative
%   rotary spectra,with the spectrum for positive frequencies +F returned
%   in SPP, and for negative frequencies -F in SNN.  
%
%   With ALPHA=1, the oscillatory Matern becomes the complex Ornstein-
%   Uhlenbeck process.
%
%   Note that NU has units of radians per sample interval DT.
%   __________________________________________________________________
%
%   Extended Matern
%  
%   [F,S]=MATERNSPEC(DT,N,SIGMA,ALPHA,LAMBDA,0,MU) with seven arguments
%   returns the spectrum of the four-parameter "extended" Matern process:
%
%      S(F) = SIGMA^2 * BESSELK(ALPHA,SQRT(F^2+LAMBDA^2)/MU)
%                                    / (SQRT(F^2+LAMBDA^2)/MU)^ALPHA * C               
%
%   where C is a normalizing constant dependent upon ALPHA, LAMBDA, and MU. 
%   The additional parameter, MU, has units of frequency.
%
%   [F,SPP,SNN]=MATERNSPEC(DT,N,SIGMA,ALPHA,LAMBDA,NU,MU) shifts the 
%   extended Matern spectrum to be centered at F=NU rather than F=0.
%
%   With MU set to zero, this becomes the standard Matern spectrum.
%   __________________________________________________________________
%
%   Damped exponential 
%
%   [F,S]=MATERNSPEC(DT,N,SIGMA,-1/2,LAMBDA,0,MU) with ALPHA set to -1/2
%   returns the spectrum of the damped exponential process, having the form
%
%      S(F) = SIGMA^2 * EXP(-MU * SQRT(F^2+LAMBDA^2)) * C 
%                              
%   where C is again a normalizing constant, dependent upon LAMBDA and MU. 
%
%   This is a special case of the extended Matern process.  As for that 
%   process, setting NU to a nonzero value results in a shifted spectrum.
%   __________________________________________________________________
%
%   Composite Matern
%
%   [F,SPP,SNN]=MATERNSPEC(DT,N,SIGMA,ALPHA,LAMBDA,NU,MU,'composite') 
%   implements the "composite" Matern spectrum having the form
%
%       SPP(F) = B * SIGMA^2 / (F^2 + MU^2)^ALPHA / [(F-NU)^2 + LAMBDA^2] 
%       SNN(F) = B * SIGMA^2 / (F^2 + MU^2)^ALPHA / [(F+NU)^2 + LAMBDA^2] 
%
%   where B is a normalizing constant discussed shortly.  This consists of 
%   a Matern spectrum times an oscillatory Matern spectrum having ALPHA=1.
%
%   The quantity in square brackets is recognized as the transfer function
%   for a damped simple harmonic oscillator.  In oceanographic terms, this 
%   composite model gives the spectrum of a damped slab model of the 
%   surface mixed layer forced by winds having a Matern spectrum.
%
%   The interpretation of the variance SIGMA is different from the other 
%   cases in MATERNSPEC, because an analytic form of the total variance 
%   does not exist. Instead SIGMA^2 is an approximation to the variance
%   associated with the oscillatory peak at F=NU.  
%
%   The additional parameter here, MU, has units of *frequency* and is the
%   damping parameter associated with the background process, which in this
%   case reprents the structure of the wind spectrum.
%
%   Here B = 2 * LAMBDA * (NU^2 + MU^2) is a normalizing constant that lets
%   SIGMA^2 be interpreted as an approximation to the inertial variance. 
%   __________________________________________________________________
%
%   See also MATERNCOV, MATERNOISE, MATERNFIT, BLURSPEC.
%
%   'maternspec --f' generates some sample figures.
%   Tests for MATERNSPEC can be found in MATERNCOV.
%
%   Usage:  [f,s]=maternspec(dt,N,sigma,alpha,lambda);
%           [f,spp,snn]=maternspec(dt,N,sigma,alpha,lambda);
%           [f,spp,snn]=maternspec(dt,N,sigma,alpha,lambda,nu);
%           [f,spp,snn]=maternspec(dt,N,sigma,alpha,lambda,nu,mu);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2015 J.M. Lilly --- type 'help jlab_license' for details


%   No longer supported, sorry
%   Reverse Matern
%
%   [F,S]=MATERNSPEC(N,SIGMA,ALPHA,H,'reverse') returns the spectrum of the
%   reverse Matern process, in which the roles of the spectrum and the
%   autocovariance functions are swapped.  The spectrum is given by
%
%      S(F) = SIGMA^2 * C * (F/H)^(ALPHA-1/2) * BESSELK(ALPHA-1/2,F/H)
%
%   where C = 2 * SQRT(PI) / GAMMA(ALPHA) / 2^(ALPHA-1/2) / H is a 
%   normalizing constant.  Note that the parameter H has been inverted 
%   so that it still has units of frequency. 
%
%   [F,SPP,SNN]=MATERNSPEC(N,SIGMA,ALPHA,C,'reverse') for complex C also
%   works and gives the spectrum of an oscillatory reverse Matern.  
%   __________________________________________________________________


if strcmpi(varargin{1}, '--f')
    type makefigs_maternspec
    makefigs_maternspec;
    return
end

cores='serial';
model='forward';

for i=1:3
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'rev')||strcmpi(varargin{end}(1:3),'com')
            model=varargin{end}(1:3);
        else
            cores=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end


if strcmpi(cores(1:3),'par')
    if exist('parpool')~=2
        disp('Sorry, parallel algorithm requires the Parallel Computing Toolbox.')
        disp('Defaulting to the standard algorithm.')
        cores='serial';
    end
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

if length(N)==1
    [omega,Spp,Snn]=maternspec_one(dt,N,sigma,alpha,lambda,nu,mu,model);
else
    if strcmpi(cores(1:3),'ser')
        for i=1:length(N)
            [omega{i,1},Spp{i,1},Snn{i,1}]=maternspec_one(dt,N(i),sigma(i),alpha(i),lambda(i),nu(i),mu(i),model);
        end
    elseif strcmpi(cores(1:3),'par')
        parfor i=1:length(N)
            [omega{i,1},Spp{i,1},Snn{i,1}]=maternspec_one(dt,N(i),sigma(i),alpha(i),lambda(i),nu(i),mu(i),model);
        end
    end
end

varargout{1}=omega;
varargout{2}=Spp;
varargout{3}=Snn;

    
function[omega,Spp,Snn]=maternspec_one(dt,N,sigma,alpha,lambda,nu,mu,model)
%This loops and implements positive/negative rotary spectra
   
omega=fourier(dt,N);
Spp=zeros(length(omega),length(lambda));
Snn=zeros(length(omega),length(lambda));

for i=1:length(lambda)
    if sigma(i)~=0
        Spp(:,i)=maternspec_spec(omega,sigma(i),alpha(i),lambda(i),nu(i),mu(i),model);
        if nu(i)==0
            Snn(:,i)=Spp(:,i);
        else
            Snn(:,i)=maternspec_spec(omega,sigma(i),alpha(i),lambda(i),-nu(i),mu(i),model);
        end
    end
end


function[S]=maternspec_spec(omega,sigma,alpha,lambda,nu,mu,model)
%This implements the single/composite and forward/reverse/exponential options

if strcmpi(model(1:3),'com')
    fact=2*sigma.^2.*lambda.*(nu.^2+mu.^2).^alpha;
    S1=frac(fact,(omega.^2+mu.^2).^alpha);
    S2=frac(1,(omega-nu).^2+lambda.^2);
    S=S1.*S2;
else
    if mu==0
        d=frac(lambda.^(2*alpha-1),materncfun(alpha));
        S=frac(sigma.^2,((omega-nu).^2+lambda.^2).^alpha).*d;
    else
        if alpha==-1/2
            fact=lambda.*besselk(1,lambda./mu);
            S=frac(pi*sigma.^2,fact).*exp(-sqrt((omega-nu).^2+lambda.^2)./mu);
        else
            fact=sigma.^2.*frac(sqrt(2*pi),mu).*frac((lambda./mu).^(alpha-1/2),besselk(abs(alpha-1/2),abs(lambda./mu)));
            omnorm=sqrt((omega-nu).^2+lambda.^2)./mu;
            S=fact*omnorm.^(-alpha).*besselk(-alpha,omnorm);
        end
    end
end


%elseif strcmpi(model(1:3),'rev')
%    S=2*pi*materncfun(alpha).*frac(sigma.^2,lambda).*maternfun(alpha-1/2,(omega-nu)./lambda);
% elseif strcmpi(model(1:3),'com')
%     A2=2*sigma.^2.*lambda.*(nu.^2+G.^2).^alpha;
%     %if strcmpi(model(1:3),'for')
%         S1=frac(A2,(omega.^2+G.^2).^alpha);
%     %else
%         %Currently generating an error for this one
%         %This one doesn't really work yet... not sure it should
%         %S1=A2.*maternfun(alpha-1/2,(omega-nu)./lambda);
%     %end
%     S2=frac(1,(omega-nu).^2+lambda.^2);
%     S=S1.*S2;

%fact=frac(2*sigma.^2.*sqrt(pi),gamma(alpha).*2.^(alpha-1/2).*lambda);
%S=fact.*(abs(omega-nu)./lambda).^(alpha-1/2).*besselk(alpha-1/2,abs(omega-nu)./lambda);
   
 
%\*************************************************************************
    

%uv=fillbad(uv);
%uv=uv-mean(uv);

% tr=1/10000;
% %tr=1;
% [a,alpha,h] = maternfit(tr,uv);
% [psi,lambda] = sleptap(length(uv),3);
% [f,spp,snn,spn] = mspec(tr,uv,psi,lambda,'adaptive');
% [fm,sppm,snnm] = maternspec(tr,length(uv),a,alpha,h);
% 
% clf
% h=twospecplot(f,[spp sppm],[snn snnm]);
% axes(h(1)),xlin,axes(h(2)),xlin

