function[varargout]=maternspec(varargin)
%MATERNSPEC  Fourier spectrum of the Matern random process and variations.
% 
%   [F,S]=MATERNSPEC(DT,N,SIGMA,ALPHA,LAMBDA) returns the spectrum S of a  
%   length N real-valued or complex-valued Matern random process having 
%   variance SIGMA^2, slope parameter ALPHA, and damping parameter LAMBDA.
%
%   DT is the sample interval.  Note that LAMBDA is understood to have the
%   same units as the inverse sample interval 1/DT.
%
%   F is an array of one-sided (positive) Fourier frequencies for a time
%   series of length N, F=FOURIER(N), where F is a *radian* frequency. 
%
%   The lengths of the output variables F and S are N/2+1 for even N, and
%   (N+1)/2 for odd N.
%     
%   S is the postive or negative rotary spectrum given by
%
%        S(F) = SIGMA^2 / (F^2+LAMBDA^2)^ALPHA 
%                           * LAMBDA^(2*ALPHA-1) / C
%
%   where C is a normalizing constant dependent upon ALPHA.  Note that the
%   positive and negative spectra are identical for this process.
%
%   For LAMBDA=0, the Matern spectrum reduces to the spectrum of fractional
%   Brownian motion.  
%
%   For details on the Matern process and its spectrum, see:
%
%     Lilly, Sykulski, Early, and Olhede, (2017).  Fractional Brownian
%        motion, the Matern process, and stochastic modeling of turbulent 
%        dispersion.  Nonlinear Processes in Geophysics, 24: 481--514.
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
%   rotary spectra, with the spectrum for positive frequencies +F returned
%   in SPP, and for negative frequencies -F in SNN.  
%
%   With ALPHA=1, the oscillatory Matern becomes the complex Ornstein-
%   Uhlenbeck process.
%
%   Note that NU has units of radians per sample interval DT.
%
%   The oscillatory Matern is described in Lilly et al. (2017).
%   __________________________________________________________________
%
%   Experimental extensions
%
%   The remaining features are experimental extensions to the Matern 
%   process.  They are not yet documented in a publication, and should
%   be considered as 'beta features' that are to be used with caution.
%   __________________________________________________________________
%
%   Generalized Matern
%
%   [F,S]=MATERNSPEC(DT,N,SIGMA,ALPHA,LAMBDA,GAMMA,'general') returns the
%   spectrum of the generalized Matern, or generalized Whittle-Matern,
%   process, introduced by Lim and Teo (2009)a,b.
%
%   This spectrum is defined as
%
%    S(F) = SIGMA^2 / [F^(2*GAMMA)+LAMBDA^(2*GAMMA)]^(ALPHA/GAMMA)
%                   * LAMBDA^(2*ALPHA-1) / C
%
%   where C is a normalizing constant dependent upon ALPHA and GAMMA. 
%
%   Here we use slightly different choices of parameters from  Lim and Teo.
%   With our choices, GAMMA=1 corresponds to the usual Matern process.
%
%   [F,SPP,SNN]=MATERNSPEC(DT,N,SIGMA,ALPHA,LAMBDA,GAMMA,NU,'generalized') 
%   is an oscillatory version of the generalized Matern process, shifting 
%   the spectrum to be centered on F=NU rather than F=0.
%   __________________________________________________________________
%
%   Extended Matern
%  
%   [F,S]=MATERNSPEC(DT,N,SIGMA,ALPHA,LAMBDA,MU,'extended') returns the 
%   spectrum of what we will call the extended Matern process:
%
%      S(F) = SIGMA^2 * BESSELK(ALPHA,MU*SQRT(F^2+LAMBDA^2))
%                                    / (MU*SQRT(F^2+LAMBDA^2))^ALPHA / C              
%
%   where C is a normalizing constant dependent upon ALPHA, LAMBDA, and MU. 
%   The additional parameter, MU, has units of time.  Here ALPHA can 
%   take on any real value, unlike for the standard Matern case.
%
%   [F,SPP,SNN]=MATERNSPEC(DT,N,SIGMA,ALPHA,LAMBDA,MU,NU,'extended') is an
%   oscillatory version of the extended Matern process.
%
%   As MU increases with ALPHA>1/2, this becomes the standard Matern form.
%   __________________________________________________________________
%
%   Damped exponential 
%
%   A special case of the extended Matern is the damped exponential.
%
%   [F,S]=MATERNSPEC(DT,N,SIGMA,-1/2,LAMBDA,MU,'extended') with ALPHA=-1/2
%   returns the spectrum of the damped exponential process, with spectrum
%
%      S(F) = SIGMA^2 * EXP(-MU * SQRT(F^2+LAMBDA^2)) / C 
%                              
%   where C is again a normalizing constant, dependent upon LAMBDA and MU. 
%
%   This is a special case of the extended Matern process.  As for that 
%   process, setting NU to a nonzero value results in a shifted spectrum.
%   __________________________________________________________________
%
%   Composite Matern
%
%   [F,SPP,SNN]=MATERNSPEC(DT,N,SIGMA,ALPHA,LAMBDA,MU,NU,'composite') 
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
%   Real-valued processes
%
%   By default MATERSPEC returns the spectrum of a complex-valued process. 
%
%   MATERNSPEC(...,'real') instead returns the spectrum of a real-valued
%   process. This also works with any of extended versions described above.
%
%   In this case, the rotary spectra SPP and SNN be forced to be the same 
%   in models that return them.
%   __________________________________________________________________
%
%   See also MATERNCOV, MATERNIMP, MATERNOISE, MATERNFIT, BLURSPEC.
%
%   'maternspec --f' generates some sample figures.
%
%   Tests for MATERNSPEC can be found in MATERNCOV.
%
%   Usage:  [f,s]=maternspec(dt,N,sigma,alpha,lambda);
%           [f,spp,snn]=maternspec(dt,N,sigma,alpha,lambda);
%           [f,spp,snn]=maternspec(dt,N,sigma,alpha,lambda,nu);
%           [f,spp,snn]=maternspec(dt,N,sigma,alpha,lambda,mu,'extended');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2020 J.M. Lilly --- type 'help jlab_license' for details


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
model='standard';
flag='complex';

for i=1:3
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'sta')||strcmpi(varargin{end}(1:3),'com')||strcmpi(varargin{end}(1:3),'ext')||strcmpi(varargin{end}(1:3),'gen')
            model=varargin{end}(1:3);
        elseif strcmpi(varargin{end}(1:3),'rea')
            flag=varargin{end};
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

parcell=varargin(3:end);
%this will be easier as a matrix
M=max(cellength(parcell));
params=zeros(M,length(parcell));
for i=1:length(parcell)
    params(:,i)=parcell{i};
end

%dt,N,sigma,alpha,lambda,nu,mu
%arrayify(sigma,alpha,lambda,nu,mu);
    
if length(N)==1
    [omega,Spp,Snn]=maternspec_one(dt,N,params,model,flag);
else
    if strcmpi(cores(1:3),'ser')
        for i=1:length(N)
            [omega{i,1},Spp{i,1},Snn{i,1}]=maternspec_one(dt,N(i),params(i,:),model,flag);
        end
    elseif strcmpi(cores(1:3),'par')
        parfor i=1:length(N)
            [omega{i,1},Spp{i,1},Snn{i,1}]=maternspec_one(dt,N(i),params(i,:),model,flag);
        end
    end
end

varargout{1}=omega;
varargout{2}=Spp;
varargout{3}=Snn;

    
function[omega,Spp,Snn]=maternspec_one(dt,N,params,model,flag)
%This loops and implements positive/negative rotary spectra
   
omega=fourier(dt,N);
Spp=zeros(length(omega),size(params,1));
Snn=zeros(length(omega),size(params,1));

for i=1:size(params,1)
    Spp(:,i)=maternspec_spec(omega,params(i,:),model);
    Snn(:,i)=maternspec_spec(-omega,params(i,:),model);
    if strcmpi(flag(1:3),'rea')
        Spp(:,i)=frac(1,2)*(Spp(:,i)+Snn(:,i));
        Snn(:,i)=Spp(:,i);
    end
end

function[S]=maternspec_spec(omega,params,model)
%This implements the single/composite and extended/exponential options


sigma=params(1);
alpha=params(2);
lambda=params(3);

%switch to extended if alpha<1/2 
if strcmpi(model(1:3),'sta') && alpha<1/2  
    model='extended';
end
 
nu=0;
if strcmpi(model(1:3),'com')
    %Composite Matern form
    mu=params(4);
    nu=params(5);
    %[sigma,alpha,lambda,mu,nu]
    fact=2*sigma.^2.*lambda.*(nu.^2+mu.^2).^alpha;
    S1=frac(fact,(omega.^2+mu.^2).^alpha);
    S2=frac(1,(omega-nu).^2+lambda.^2);
    S=S1.*S2;
elseif strcmpi(model(1:3),'gen')
    %Generalized Matern form
    gamma=params(4);
    if length(params)==5
        nu=params(5);
    end
    omega=omega-nu;
 %   sigma,alpha,lambda,gamma,nu
    d=frac(lambda.^(2*alpha-1),materncfun(alpha,gamma));
    S=frac(sigma.^2,(abs(omega).^(2*gamma)+lambda.^(2*gamma)).^(alpha./gamma)).*d;
elseif strcmpi(model(1:3),'ext')
    %Extended Matern, reduces to the standard Matern form for mu=0
    mu=params(4);
    if length(params)==5
        nu=params(5);
    end
    omega=omega-nu;
    omnorm=mu*sqrt(omega.^2+lambda.^2);
    S1=frac(sigma.^2,lambda.*materncfun(alpha).*(omega.^2/lambda.^2+1).^alpha);
    S2=frac(maternfun(alpha+1/2,omnorm),maternfun(alpha,lambda*mu));
    S=S1.*S2;
elseif strcmpi(model(1:3),'sta') 
    %Standard Matern form
    if length(params)==4
        nu=params(4);
    end
    omega=omega-nu;
    d=frac(lambda.^(2*alpha-1),materncfun(alpha));
    S=frac(sigma.^2,(omega.^2+lambda.^2).^alpha).*d;
end

%Extended Matern form reduces to this exponential form, see notes
%if alpha==-1/2
%   fact=lambda.*besselk(1,lambda.*mu);
%   S=frac(pi*sigma.^2,fact).*exp(-sqrt(omega.^2+lambda.^2).*mu);
%end
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


 
%\*************************************************************************
    
% % 
% function[]=maternspec_figures;
% 
% %[f,s]=maternspec(dt,N,sigma,alpha,lambda);
% %[f,spp,snn]=maternspec(1,1000,1,1,1/10,4/10);
% %figure,plot(f,spp),xlog,ylog
% 
% [f,sppo,snno]=maternspec(1,1000,2,1,1/10,0,0);
% %[f,spp,snn]=maternspec(1,1000,2,1,1/10,0,10);
% [f,spp,snn]=maternspec(1,1000,2,1,1/10,0,10);
% figure,plot(f,[sppo spp]),xlog,ylog



