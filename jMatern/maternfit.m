function [structout] = maternfit(varargin)
%MATERNFIT  Parametric spectral fit to the Matern form. [with A. Sykulski]
%
%   MATERNFIT performs a parametric fit of the spectrum of a time series to
%   that expected for a Matern process plus optional spin.  
%
%   The Matern process plus spin has a spectrum given by, see MATERNSPEC, 
%
%        SPP(F) = SIGMA^2 / [(F-NU)^2/LAMBDA^2 + 1]^ALPHA 
%                                    / (LAMBDA * MATERNC(ALPHA))
%
%   where SIGMA is the standard deviation, NU is a frequency shift, LAMBDA 
%   is a damping coefficient, and ALPHA is 1/2 of the spectral slope.  The
%   coefficient on the second line lets us parameterize the spectrum in 
%   terms of the process variance SIGMA^2.
%
%   The optimal parameters are found using a frequency-domain maximum
%   likelihood method, accounting for both aliasing and spectral blurring.
%
%   For further details on the parameter inference method method, see
%
%     Sykulski, Olhede, Lilly, and Danioux (2016).  Lagrangian time series 
%        models for ocean surface drifter trajectories. Journal of the 
%        Royal Statisical Society, Series C. 65 (1): 29--50.
%
%     Sykulski, Olhede, and Lilly (2016). The de-biased Whittle likelihood 
%        for second-order stationary stochastic processes.  ArXiv preprint.
%
%   For details on the Matern process and its spectrum, see:
%
%     Lilly, Sykulski, Early, and Olhede, (2017).  Fractional Brownian
%        motion, the Matern process, and stochastic modeling of turbulent 
%        dispersion.  Nonlinear Processes in Geophysics, 24: 481--514.
%   __________________________________________________________________
%
%   Usage
%
%   FIT=MATERNFIT(DT,Z,FO) where Z is a complex-valued time series, returns
%   the result of a fitting the periodogram of Z to a background process
%   having a non-shifted Matern form, fit over all frequencies. The range 
%   of frequencies involved in the fit can be modified, as described below.
%
%   DT is the sample interval, which has units of days, while FO is a
%   signed reference frequency in units of radians per day.  FO is used in 
%   determining search ranges and initial guesses.  In oceanographic
%   applications, one would normally choose FO as the Coriolis frequency.
%
%   The output argument FIT is a structure which will be described later.
%
%   The default search ranges and initial guess values are as follows:
% 
%                    Low  Guess  High 
%        SIGMA   =   [0      1    100] * STD(Z)
%        ALPHA   =   [1/2    1     10] 
%        LAMBDA  =   [1e-3  2e-3   10] * ABS(FO)
%        NU      =   [0      0      0] * ABS(FO)
%
%   These can be modified, as described below. 
%   __________________________________________________________________
%
%   Output format
%
%   The following parameters are output as fields of the structure FIT.
%
%      SIGMA     Standard deviation of currents in cm/s
%      ALPHA     Slope parameter 
%      LAMBDA    Damping parameter in rad / day
%      NU        Oscillation frequency in rad /day
%      RANGE     Ranges of fit search for each parameter
%      PARAMS    Sub-structure with fields described below
%
%   The associated spectra can be then created from SIGMA, ALPHA, LAMBDA, 
%   and NU using MATERNSPEC. 
%
%   RANGE is substructure with fields SIGMA, ALPHA, LAMBDA, and NU, such 
%   that RANGE.SIGMA specifies the *dimensional* values of the associated
%   range and initial guess with format [MIN GUESS MAX], and so forth for 
%   all the other parameters. 
%
%   PARAMS is a substructure containing various parameter values
%   characterizing the fit itself: 
%
%      PARAMS.DT      Sample rate, as input to MATERNFIT
%      PARAMS.FO      Reference frequency, as input to MATERNFIT
%      PARAMS.A       Index into first frequency F(A) used in the fit
%      PARAMS.B       Index into last frequency F(B) used in the fit  
%      PARAMS.P       Number of free parameters used in the fit 
%      LIKE           Negative of the log-likelihood 
%      AICC           Akaike Information Criterion, corrected version 
%      ERR            Normalized error of fit to log spectra 
%      EXITFLAG       The exit flag from the optimization routine
%      ITER           The number of iterations in the optimization routine
%      PARAMS.SIDE    String specifying frequency side options
%      PARAMS.ALG     String specifying algorithm options
%      PARAMS.VER     String specifiying raw or difference option
%      PARAMS.CORES   String specifying series or parallel computation
%  
%   These parameters are described in more detail below.
%   __________________________________________________________________
%   
%   Multiple input time series
%
%   FIT=MATERNFIT(DT,Z,FO) may have Z being a matrix with N columns, or a 
%   cell array of N different time series.  In both of these cases, DT and 
%   FO may be scalars or length N arrays.
%
%   In these cases, the fields SIGMA, ALPHA, LAMBDA, and NU of FIT will 
%   also be arrays with N elements, as will LIKE, AICC, ERR, and EXITFLAG. 
%   The subfields of RANGE will all be N x 3 arrays.  In PARAMS, the fields
%   DT, FO, A, and B will all be length N.   
%   __________________________________________________________________
%   
%   Specifiying frequencies
%  
%   MATERNFIT(DT,Z,FO,RA,RB), where RA and RB are both real-valued scalars,
%   applies the fit to only frequencies in the range
%
%               ABS(FO)*RA < F < ABS(FO)*RB
%
%   with the default behavior corresponding to MATERNFIT(DT,Z,FO,0,INF).
%   Thus RA is the smallest permissible ratio of F to the reference 
%   frequency, while RB is similarly the largest permissible ratio.
%
%   MATERNFIT(DT,Z,FO,RA,[RB,RN]) also works, where the fifth argument is
%   an array of length two. In this case the fit is applied to the range
%   
%               ABS(FO)*RA < F < MIN( ABS(FO)*RB, PI/DT*RN )
% 
%   RB is the largest permissible ratio of F to the reference frequency, 
%   while RN is the largest permissible ratio of F to the Nyquist PI/DT.  
%   The fit extends to the smaller of these two frequencies.  
%
%   If RA and RB are imaginary numbers, rather than real numbers, then 
%   the fit is only applied to frequencies in the range
%
%               IMAG(RA) < F < IMAG(RB)
%
%   that is, the range is found without scaling by the reference frequency. 
%   
%   MATERNFIT can create the fit by utilizing both positive and negative 
%   frequency sides of the spectrum, the default behavior, or to only one 
%   side (plus the zero frequency). This is modified as follows:
% 
%      MATERNFIT(...,'both',...), the default, uses both sides.
%      MATERNFIT(...,'positive',...), uses the side where F/FO is positive.
%      MATERNFIT(...,'negative',...), uses the side where F/FO is negative.
%
%   Thus, changing the sign of the reference frequency also changes the 
%   side of the spectrum to be fit using the 'postive' or 'negative' flags.
%   __________________________________________________________________
%   
%   Parameter range specification
%   
%   The default search ranges and guess values can be modified.  As an 
%   example, to modify the default range and guess for the SIGMA parameter
%   corresponding to the background flow, use
%
%       MATERNFIT(...,'range.sigma',[MIN GUESS MAX],...)                   
%  
%   and so forth for other parameters. SIGMA ranges input to MATERNFIT 
%   represent *fractions* of the total signal standard deviation.  Thus 
%
%       MATERNFIT(...,'range.sigma',[0 1 100],...)
%
%   corresponds to the default setting.  The upper limit is set very high
%   because occasionally an optimum spectral fit is found that has much 
%   larger total variance than the signal. 
%
%   Similarly, ranges of LAMBDA and NU are nondimensional values representing
%   fractions of the magnitude of the reference frequency ABS(FO).  Thus 
%
%       MATERNFIT(...,'range.lambda',[1/1000 2/1000 1],...)
%
%   corresponds to the default range for LAMBDA. 
%
%   The slope parameter ALPHA can take on a value of no less than 1/2, so
%   the lower range for ALPHA cannot be less than 1/2.
%
%   This approach can be used to omit parameters from the fit, by setting
%   the MIN, GUESS, and MAX values to be identical.  For example,  
%
%      MATERNFIT(...,'range.alpha',[2 2 2]) 
%   
%   sets the ALPHA parameter to a value of 2.  Fixed parameters are then 
%   not included in the optimization, thus speeding up the fit.  
%   __________________________________________________________________
%   
%   Parameter value specification
%   
%   Parameters can also be set to particular values for each time series. 
%   This is done by setting the 'value' field as follows
%
%       MATERNFIT(...,'value.sigma',SIGMA,...)
%    
%   and so forth for the other parameters.  
%
%   Here SIGMA is an array of the same length as the number of time series
%   in Z.  Thus SIGMA is a scalar if Z is a single array, an array of
%   LENGTH(Z) is Z is a cell array, or of SIZE(Z,2) if Z is a matrix. 
%
%   Note that the values specified in this way are actual dimensional 
%   values, not nondimensional values as with setting the range. 
%
%   This approach works by internally setting the dimension MIN, GUESS, and
%   MAX to the value specified for each time series, overriding the default 
%   choices. This will be reflected in the output RANGE fields. 
%
%   If both ranges and values are specified for the same parameter, the
%   value settings take precedence.
%   __________________________________________________________________
%   
%   Numerical options
%
%   MATERNFIT has options for specifying numerical details of the fit and 
%   the optimization algorithm. 
%
%   MATERNFIT can use one of two different optimization algorithms.
%
%     MATERNFIT(...,'bnd',...), the default, uses FMINSEARCHBND by J. 
%       D'Errico included with JLAB in accordance with its license terms. 
%       This in turn calls Matlab's FMINSEARCH using Nelder-Mead.
%
%     MATERNFIT(...,'con',...) alternately uses FMINCON, using the default
%       interior-point algorithm.  This requires Matlab's Optimization 
%       Toolbox to be installed.  Again, this is mainly for testing.
%
%   In tests, FMINCON is generally faster, and most of the fits agree
%   closely with those using FMINSEARCHBND.  However, occasionally FMINCON
%   produces fits that are significantly worse than those obtained from 
%   FMINSEARCHBND, which is why the latter is preferred by default.
%
%   MATERNFIT can employ two different versions of the fit. 
%     
%     MATERNFIT(...,'difference',...) estimates the Matern parameters by 
%       fitting the first difference of the time series to the first 
%       difference of a Matern process.  This amounts to a form of pre-
%       whitening, and is the default behavior when a taper is not input.
%     MATERNFIT(...,'raw',...) fits the time series directly to a Matern. 
%       This is the default behavior when a taper *is* input.
%
%   These choices are reflected in the output fields PARAMS.ALG and 
%   PARAMS.VER, respectively.
%   __________________________________________________________________
%
%   Tapering
%
%   The default behavior of fitting the first difference of the spectrum 
%   is usually sufficient to account for spectral blurring.  For very steep 
%   spectra or those with a very large dynamic range, this is no longer the
%   case, because leakage from high-energy portions of the spectra will 
%   obscure the structure of low-energy portions. 
%  
%   To addess this, MATERNFIT can optionally perform the fit by fitting 
%   to the tapered and aliased spectrum, correctly accounting for the 
%   influence of tapering.  This is accomplished with 
%   
%      MATERNFIT(...,'taper',PSI,...)
%
%   where PSI is a data taper of the same length as Z.  If Z is a cell
%   array, then PSI is a cell array of data tapers having the same length
%   as the components of Z.  PSI is typically computed by SLEPTAP.
%
%   When a taper is input, the default behavior is *not* to difference 
%   the time series, corresponding to the 'raw' option.  If both a taper
%   and the 'difference' flag are both input, then the taper length must 
%   be one less than the data length.  
%   __________________________________________________________________
% 
%   Error
%
%   While the maximum likelikehood method is not about finding the best 
%   fit in a least squares sense, a measure of the mean squared error 
%   provides a useful measure of evaluating the degree of misfit.
%
%   The error ERR returned by MATERNFIT is the squared difference between 
%   the *natural log* of the periodogram and that of the fit spectrum, 
%   summed over all frequencies used in the fit, and divided by the sum of 
%   the squared natural log of the periodogram over the same frequencies.
%
%   The error is computed over the inertial side, anti-inertial side, or 
%   both sides, depending on which frequency range is used for the fit. 
%
%   When the 'difference' behavior is employed, as is the default, ERR will
%   reflect the error between the periodogram of the first difference of 
%   the time series and the spectrum of the first difference of the fit.
%   __________________________________________________________________
%   
%   Parallelization
%
%   With Matlab's Parallel Computing toolbox installed, when Z is a cell 
%   array or a matrix, MATERNFIT(DT,Z,'parallel') will loop over the 
%   elements of Z using a parallel "parfor" loop to speed things up.
%
%   This choice is reflected in the output field PARAMS.CORES, which 
%   will take on the values 'series' (the default) or 'parallel'.
%   __________________________________________________________________
%
%   'maternfit --t' runs a test.
%
%   Usage: fit=maternfit(dt,z);
%          fit=maternfit(dt,z,fc,a,b);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2017 J.M. Lilly and A.M. Sykulski
%                                --- type 'help jlab_license' for details

%Notes to self:
%Note, this is all set up to run the 5-parameter Matern process, I have
%just not checked the code.  Also make sure to comment, which I can grab 
%from SPECFIT, and also make sure to assign bgtype in output params.
%
%This is basically just a stripped-down version of SPECFIT.  Use Mac's
%Filemerge tool to compare them. 

if nargin>0
    if strcmpi(varargin{1}, '--f')
        maternfit_fig,return
    elseif strcmpi(varargin{1}, '--t')
        maternfit_test,return
    end
end
%--------------------------------------------------------------------------
%Set search ranges for parameters
range=struct;

%Ranges for the background process
range.background.sigma=[0 1 100];            %Range of standard deviation as fraction of total
range.background.alpha=[1/2 1 10];           %Range of slope parameter
range.background.lambda=[1/1000 2/1000 10];  %Range of decay parameter as multiple of Coriolis
range.background.nu=[0 0 0];                 %Frequency shift is fixed at NU=0
range.background.mu=[0 0 0];                 %Range of decay parameter as multiple of inverse Coriolis
%--------------------------------------------------------------------------
%Initialize value structure
names=fieldnames(range);
for i=1:length(names)
    subnames=fieldnames(eval(['range.' names{i}]));
    for j=1:length(subnames)
        eval(['value.' names{i} '.' subnames{j} '=[];'])
    end
end
%--------------------------------------------------------------------------
dt=double(varargin{1});
z=varargin{2};
varargin=varargin(3:end);
%--------------------------------------------------------------------------
%Sorting out input options 

opts.side='both';            %Determines side: opposite, same, or both
opts.alg='bnd';              %Algorithm: bnd or con
opts.ver='difference';       %Fit version: difference or raw
opts.cores='series';         %Series or parallel computation
opts.bgtype='matern';        %Type of model for background

psi=[];                      %The default data taper is the empty taper
verstrin=false;              %Flag for whether or not version string is input
for i=1:30
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'con')||strcmpi(varargin{end}(1:3),'bnd')
            opts.alg=varargin{end};
            if strcmpi(opts.alg(1:3),'con')
                if exist('fmincon.m')~=2
                    error('Sorry, MATERNFIT with ALG=''con'' requires the Optimization Toolbox')
                end
            end
        elseif strcmpi(varargin{end}(1:3),'raw')||strcmpi(varargin{end}(1:3),'dif')||strcmpi(varargin{end}(1:3),'sec')
            opts.ver=varargin{end};
            verstrin=true;
        elseif strcmpi(varargin{end}(1:3),'bot')||strcmpi(varargin{end}(1:3),'pos')||strcmpi(varargin{end}(1:3),'neg')
            opts.side=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'ser')||strcmpi(varargin{end}(1:3),'par')
            opts.cores=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'sta')||strcmpi(varargin{end}(1:3),'com')
            opts.model=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'mat')||strcmpi(varargin{end}(1:3),'exp')||strcmpi(varargin{end}(1:3),'ext')
            opts.bgtype=varargin{end};
        end
        varargin=varargin(1:end-1);
    elseif length(varargin)>1
        if ischar(varargin{end-1})&&~ischar(varargin{end})
            if strcmpi(varargin{end-1}(1:3),'tap')
                psi=varargin{end};
            elseif strcmpi(varargin{end-1}(1:3),'ran')
                name=varargin{end-1}(7:end);
                range.background=setfield(range.background,name,varargin{end});
            elseif strcmpi(varargin{end-1}(1:3),'val')
                name=varargin{end-1}(7:end);
                value.background=setfield(value.background,name,varargin{end});
            end
            varargin=varargin(1:end-2);
        end
    end
end
%--------------------------------------------------------------------------
%Convert value structure to a cell array
if iscell(z)
    N=length(z);
else
    N=size(z,2);
end
names=fieldnames(range);
for i=1:length(names)
    subnames=fieldnames(eval(['range.' names{i}]));
    for j=1:length(subnames)
        subvalue=eval(['value.' names{i} '.' subnames{j}]);
        if ~isempty(subvalue)
            values(i,j)=subvalue(k);
        else
            values(i,j)=nan;
        end
    end
end
%Have to convert this to numeric values because I can't pass a structure
%through to the workers, for some unknown reason
%--------------------------------------------------------------------------
%Make sure ranges and flags are set correctly for specified background type 
if strcmpi(opts.bgtype(1:3),'exp')
    range.background.alpha=[-1/2 -1/2 -1/2];
    range.background.mu=[1/100 10 100];  
    %range.background.lambda=[1/100000 2/1000 10];  %Range of decay parameter as multiple of Coriolis
elseif strcmpi(opts.bgtype(1:3),'ext')
    range.background.alpha=[-1/2 1 10];
    range.background.mu=[1/100 10 100];  
elseif strcmp(opts.bgtype(1:3),'mat')
    if range.background.alpha(1)<=1/2
        %Correction for alpha being less than permitted value
        range.background.alpha(1)=1/2+1e-10;
    end
end
%--------------------------------------------------------------------------
%Frequency range
flow=0;
fhigh=inf;
fc=varargin{1};
varargin=varargin(2:end);

if length(varargin)==2
    flow=varargin{1};
    fhigh=varargin{2};
end
if length(fhigh)==1
    fhigh=[fhigh inf];
end
%--------------------------------------------------------------------------
%Re-sizing of dt and fc 
dt=dt(:);

if length(dt)==1
    if iscell(z)
        dt=dt+zeros(length(z),1);
    else
        dt=dt+zeros(size(z,2),1);
    end
end

fc=fc(:);
if length(fc)==1
    if iscell(z)
        fc=fc+zeros(length(z),1);
    else
        fc=fc+zeros(size(z,2),1);
    end
end
%--------------------------------------------------------------------------
%Difference z if requested, and compute standard deviation

%Default to raw algorithm if taper is input
if ~isempty(psi)&&~verstrin
    opts.ver='raw';
end

if ~iscell(z)
    stdz=std(z);
    if strcmpi(opts.ver(1:3),'dif')
        z=diff(z);
    elseif strcmpi(opts.ver(1:3),'sec')
        z=diff(diff(z));
    end
else
    for i=1:length(z)
        stdz{i}=std(z{i});
        if strcmpi(opts.ver(1:3),'dif')
            z{i}=diff(z{i});
        elseif strcmpi(opts.ver(1:3),'sec')
            z{i}=diff(diff(z{i}));
        end
    end
end
%--------------------------------------------------------------------------
%Compute data window from taper or periodogram
if ~isempty(psi)
    if iscell(psi)
        for i=1:length(psi)
            if size(psi{i},1)~=size(z{i},1)
                error('Length of taper does not match length of data.')
            end
            win{i}=conv(psi{i}(:),psi{i}(:));
            win{i}=win{i}(end-length(z{i})+1:end);
        end
    else
        if size(psi,1)~=size(z,1)
            error('Length of taper does not match length of data.')
        end
        win=conv(psi(:),psi(:));
        win=win(end-size(z,1)+1:end);
    end
else
    if iscell(z)
        for i=1:length(z)
            N=length(z{i});
            win{i}=[N:-1:1]'./N;
            psi{i}=[];
        end
    else
        N=size(z,1);
        win=[N:-1:1]'./N;
        psi=[];
    end
end
%--------------------------------------------------------------------------
if ~iscell(z)&&size(z,2)==1
    %Single time series input, no loop
    [x,xa,xb,xo,like,aicc,P,a,b,err,exitflag,iters]=...
        maternfit_one(dt,z,flow,fhigh,range,opts,fc,psi,win,stdz,values);
    %The rest of this simply implements different types of loops
elseif ~iscell(z)
    %Loop over matrix columns
    N=size(z,2);[x,xa,xb,xo]=vzeros(N,5);
    [like,aicc,P,a,b,err,exitflag,iters]=vzeros(N,1);
    if strcmpi(opts.cores(1:3),'ser')
        %For series loop over matrix columns
        for i=1:size(z,2)
            disp(['MATERNFIT processing series ' int2str(i) ' of ' int2str(size(z,2)) '.'])
            [x(i,:),xa(i,:),xb(i,:),xo(i,:),like(i),aicc(i),P(i),a(i),b(i),err(i),exitflag(i),iters(i)]=...
                maternfit_one(dt(i),z(:,i),flow,fhigh,range,opts,fc(i),psi,win,stdz(i),values);
        end
    elseif strcmpi(opts.cores(1:3),'par')
        %For parallel loop over matrix columns
        disp('MATERNFIT using parallel for loop...')
        parfor i=1:size(z,2)
            %dt(i),flow,fhigh,range,opts,fc(i),fit
            disp(['MATERNFIT processing series ' int2str(i) ' of ' int2str(size(z,2)) '.'])
            [x(i,:),xa(i,:),xb(i,:),xo(i,:),like(i),aicc(i),P(i),a(i),b(i),err(i),exitflag(i),iters(i)]=...
                maternfit_one(dt(i),z(:,i),flow,fhigh,range,opts,fc(i),psi,win,stdz(i),values);
        end
    end
elseif iscell(z)
    %Loop over cell array
    N=length(z);[x,xa,xb,xo]=vzeros(N,5)
    [like,aicc,P,a,b,err,exitflag,iters]=vzeros(N,1);
    if strcmpi(opts.cores(1:3),'ser')
        %For series loop over cell array
        for i=1:length(z)
            disp(['MATERNFIT processing series ' int2str(i) ' of ' int2str(length(z)) '.'])
            [x(i,:),xa(i,:),xb(i,:),xo(i,:),like(i),aicc(i),P(i),a(i),b(i),err(i),exitflag(i),iters(i)]=...
                maternfit_one(dt(i),z{i},flow,fhigh,range,opts,fc(i),psi{i},win{i},stdz{i},values);
        end
    elseif strcmpi(opts.cores(1:3),'par')
        %For parallel loop over cell array
        disp('MATERNFIT using parallel for loop...')
        parfor i=1:length(z) 
            disp(['MATERNFIT processing series ' int2str(i) ' of ' int2str(length(z)) '.'])
            [x(i,:),xa(i,:),xb(i,:),xo(i,:),like(i),aicc(i),P(i),a(i),b(i),err(i),exitflag(i),iters(i)]=...
                maternfit_one(dt(i),z{i},flow,fhigh,range,opts,fc(i),psi{i},win{i},stdz{i},values);
        end
    end
end

sigma=x(:,1);
alpha=x(:,2);
lambda=x(:,3);
nu=x(:,4);
mu=x(:,5);

clear range
range.sigma= [xa(:,1)  xo(:,1)  xb(:,1)];
range.alpha= [xa(:,2)  xo(:,2)  xb(:,2)];
range.lambda=[xa(:,3)  xo(:,3)  xb(:,3)];
range.nu=    [xa(:,4)  xo(:,4)  xb(:,4)];
range.mu=    [xa(:,5)  xo(:,5)  xb(:,5)];
    
structout=[];

use opts

make params dt fc a b P like aicc err exitflag iters side alg ver cores 
%make params dt fc a b P like aicc err exitflag iters side alg ver bgtype cores 
make structout sigma alpha lambda nu mu range params

if ~(strcmpi(opts.bgtype(1:3),'exp')||strcmpi(opts.bgtype(1:3),'ext'))
    structout=rmfield(structout,'mu');
    structout.range=rmfield(structout.range,'mu');
end

%--------------------------------------------------------------------------
function[x,exitflag,iters]=maternfit_optimize(dt,N,om,spp,snn,index,xa,xo,xb,opts)
%Optimization for all parameters using both sides
%Optimizations like to work with doubles... see comment below
xa=double(xa);
xo=double(xo);
xb=double(xb);
%xa,xo,xb

xanorm=(xa./xo);
xbnorm=(xb./xo);
xanorm(~isfinite(xanorm))=0;
xbnorm(~isfinite(xbnorm))=0;

guess=ones(size(xo));
guess((xanorm==0)&(xbnorm==0))=0;

tol=1e-6;%tol=1e-3;
if strcmpi(opts.alg(1:3),'bnd')
    options=optimset('GradObj','on','MaxFunEvals',10000,'MaxIter',10000,'TolFun',tol,'TolX',tol);
    [xf,fval1,exitflag,struct]=fminsearchbnd(@(z) specmodel(dt,z.*xo,spp,snn,index,N,opts),...
        guess,xanorm,xbnorm,options);
elseif strcmpi(opts.alg(1:3),'con')
    options=optimoptions('fmincon','Algorithm','interior-point','MaxFunEvals',10000,'MaxIter',10000,'TolFun',tol,'TolX',tol,'Display','off');
    [xf,fval1,exitflag,struct]=fmincon(@(z) specmodel(dt,z.*xo,spp,snn,index,N,opts),...
        guess,[],[],[],[],xanorm,xbnorm,[],options);
end
iters=struct.iterations;

x=xf.*xo;% scale back to correct units
[like,Spp,Snn,f]=specmodel(dt,x,spp,snn,index,N,opts);
%--------------------------------------------------------------------------
function [like,Spp,Snn,f,indexn]=specmodel(dt,x,spphat,snnhat,index,N,opts)
% Computes the value of the Whittle likelihood for parameter set X, 
% given periodogram estimates SPPHAT and SNNHAT of positive and negative
% sides of the spectrum, respectively.  Only frequencies in the locations 
% given by INDEX are used to compute the likelihood.  N is the data length.

%Optimizations like to work with doubles... this prevents a strange error 
%where, inside the optimization loop, the spectrum takes on negative values
x=double(x);

if strcmpi(opts.ver(1:3),'dif')
    N=N+1;
elseif strcmpi(opts.ver(1:3),'sec')
    N=N+2;
end

%Compute the blurred spectrum by first computing the covariance.  
R=maternfit_materncov(dt,N,x(1,:));    
[f,Spp,Snn]=maternfit_blurspec(dt,R,opts.ver,opts.win);
    
%[f,Spp2,Snn2]=blurspec(1,R,opts.ver,'window',opts.win);aresame(Spp,Spp2)
%Same answers but slower

indexn=fixindex(N,index);        %This deals with zero and Nyquist for negative frequencies
if strcmpi(opts.side(1:3),'pos')
    %Include only positive frequencies
    like=sum(log(Spp(index))+spphat(index)./Spp(index));
elseif strcmpi(opts.side(1:3),'neg')
    %Include only negative frequencies
    like=sum(log(Snn(index))+snnhat(index)./Snn(index));
else
    % vsize(Spp,spphat,index)
    likepp=sum(log(Spp(index))+spphat(index)./Spp(index));
    likenn=sum(log(Snn(indexn))+snnhat(indexn)./Snn(indexn));
    like=likepp+likenn;
end

%FMINCON likes to work with doubles
like=double(like);

%hold on,plot(f,spphat),plot(-f, Spp)
%--------------------------------------------------------------------------
%function[x,xar,xbr,like,aicc,P,f,spp,snn,Spp,Snn,a,b,err]=maternfit_one(dt,z,flow,fhigh,side,alg,ver,range,fc)
function[x,xa,xb,xo,like,aicc,P,a,b,err,exitflag,iters]=maternfit_one(dt,z,flow,fhigh,range,opts,fc,psi,win,stdz,values)

%This is pretty dumb, but Matlab doesn't want me to pass a structure 
%through parfor.  Thus I have to convert to a matrix and then back again.
value.background.sigma=values(1,1);
value.background.alpha=values(1,2);
value.background.lambda=values(1,3);
value.background.nu=values(1,4);
value.background.mu=values(1,5);

[x,xa,xb,xo]=vzeros(1,5,nan);
[like,aicc,P,a,b,err,exitflag,iters]=vzeros(1,1,nan);
if anyany(~isfinite(z))
    return
end

N=length(z);

%Put window into options to simplify argument passing
opts.win=win;

use range
%Multiply standard deviation ranges by sample standard deviation 
background.sigma  = background.sigma*stdz;

%Multiply damping and frequency ranges by Coriolis frequency
background.lambda  = background.lambda*abs(fc);     
background.mu      = background.mu.*abs(fc);   %Is this a time or frequency scale?

%Determine if any explicit values are input
names=fieldnames(value);
for i=1:length(names)
    subnames=fieldnames(eval(['value.' names{i}]));
    for j=1:length(subnames)
        subvalue=eval(['value.' names{i} '.' subnames{j}]);
        if ~isempty(subvalue)&&~isnan(subvalue)
            %Overwrite parameter ranges with explicit values
            eval([names{i} '.' subnames{j} '(1)=subvalue;'])
            eval([names{i} '.' subnames{j} '(2)=subvalue;'])
            eval([names{i} '.' subnames{j} '(3)=subvalue;'])
        end
    end
end

xa=[background.sigma(1)  background.alpha(1)   background.lambda(1)  background.nu(1)  background.mu(1)];
xo=[background.sigma(2)  background.alpha(2)   background.lambda(2)  background.nu(2)  background.mu(2)];
xb=[background.sigma(3)  background.alpha(3)   background.lambda(3)  background.nu(3)  background.mu(3)];

%Periodogram from mspec  ... use dt = 1
%vsize(dt,z,psi),length(find(isfinite(z))),isreal(z)
[om,spp,snn]=mspec(dt,z,psi);

if ~isnan(fc)
    if isreal(flow(1))
        a=find(om>=flow*abs(fc),1,'first');   %Look up multiple of Corilois frequency
    else
        a=find(om>=imag(flow),1,'first');   
    end
    if isreal(fhigh(1))
        b=find(om<min(fhigh(1)*abs(fc),fhigh(2)*pi/dt),1,'last');    %Look up multiple of Corilios frequency
    else
        b=find(om<min(imag(fhigh(1)),real(fhigh(2))*pi/dt),1,'last');  
    end
else
    a=1;
    b=length(om);
end

%figure,plot(om,spp),hold on,plot(-om,snn),ylog,vlines([om(a) om(b) -om(a) -om(b)])

%[like,Spp,Snn,f,indexn]=specmodel(dt,x(:,4),spp,snn,a:b,N,opts);
%figure,subplot(1,2,1),plot(f,[Spp spp]),ylog,vlines(x(9,2),'r')
%subplot(1,2,2),plot(f,[Snn snn]),ylog,
    
%Fit to background on requested side
%xa,xo,xb
[x,exitflag,iters]=maternfit_optimize(dt,N,om,spp,snn,a:b,xa,xo,xb,opts);
%x
 
%Final value of likelihood and spectra
[like,Spp,Snn,f,indexn]=specmodel(dt,x,spp,snn,a:b,N,opts);

% Compute the number of free parameters
P=0;
for i=1:size(x,1)
    if ~isnan(xo(i,1))
        P=P+length(find(xa(i,:)~=xb(i,:)));
    end
end

% Compute the value of the AICC information criterion
if (strcmpi(opts.side(1:3),'opp')||strcmpi(opts.side(1:3),'sam'))
    M=length(a:b);
else
    M=length(a:b)+length(indexn);
end
aicc=2*like + frac(2*M*P,M-P-1);
    
%figure,plot(f,[Spp spp]),hold on,plot(-f,[Snn snn]),ylog
%vlines(f(a)),vlines(f(b))

%Compute an error measure
if strcmpi(opts.side(1:3),'sam')
    numer=sum(squared(log(spp(a:b))-log(Spp(a:b))));
    denom=sum(squared(log(spp(a:b))));
    err=frac(numer,denom);
elseif strcmpi(opts.side(1:3),'opp')
    numer=sum(squared(log(snn(a:b))-log(Snn(a:b))));
    denom=sum(squared(log(snn(a:b))));
    err=frac(numer,denom);
else
    numerp=sum(squared(log(spp(a:b))-log(Spp(a:b))));
    denomp=sum(squared(log(spp(a:b))));
    numern=sum(squared(log(snn(a:b))-log(Snn(a:b))));
    denomn=sum(squared(log(snn(a:b))));
    err=frac(numerp+numern,denomp+denomn);
end

vtranspose(x,xa,xb,xo);
vcolon(x,xa,xb,xo);
vtranspose(x,xa,xb,xo);
%--------------------------------------------------------------------------
function[index]=fixindex(N,index)
%This modifies the frequency index appropriate for negative frequencies, 
%due to the fact that zero is repeated for both even and odd, while the  
%Nyquist is repeated for odd. See comments at MSPEC.
if iseven(N)
    if index(end)==N
        index=index(1:end-1);
    end
end
if index(1)==1
    index=index(2:end);
end
%--------------------------------------------------------------------------
function[omega,Spp,Snn]=maternfit_blurspec(dt,R,ver,win)
%This is a version of BLURSPEC, stripped down for speed.

%I don't take dt into account during difference
if strcmpi(ver(1:3),'dif')    
    R=2*R(1:end-1,:)-R(2:end,:)-[conj(R(2,:));R(1:end-2,:)];
elseif strcmpi(ver(1:3),'sec')
    R=2*R(1:end-1,:)-R(2:end,:)-[conj(R(2,:));R(1:end-2,:)];
    R=2*R(1:end-1,:)-R(2:end,:)-[conj(R(2,:));R(1:end-2,:)];
end

R=R.*win;
R(1,:)=R(1,:)./2;  %Don't forget to divide first element by two
S=dt*2*real(fft(R)); %But I do take it into account in spectrum

S=abs(S);  %Sometimes there are small negative parts after blurring
N=size(R,1);
omega=frac(1,dt)*2*pi*(0:floor(N/2))'./N;
%omega=fourier(N);

Spp=S(1:length(omega),:);
Snn=[S(1,:);S(end:-1:end-length(omega)+2,:)];
%--------------------------------------------------------------------------
function[R]=maternfit_materncov(dt,N,x)
sigma=x(1);
alpha=x(2);
lambda=x(3);
nu=x(4);
mu=x(5);

tau=dt*[0:N-1]';

if isnan(sigma)
    R=[];
else
    %This is just copied from MATERNCOV, but it is faster to have it internal
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
            R=fact*tnorm.^(alpha-1/2).*besselk(abs(alpha-1/2),tnorm);
        end
    end
    if nu~=0
        R=R.*exp(sqrt(-1)*tau*nu);
    end
    R(tau==0)=1;
    R=R.*sigma.^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of function body; begin tests and figures
function[]=maternfit_test

sigo=17;
alpha=1.5;
h=1/10;

rng(0);
dt=1;
z=maternoise(dt,1000,sigo,alpha,h);
fit=maternfit(dt,z,frac(1,2)*pi);
bool=aresame(abs(fit.sigma./sigo-1),0,0.2)&&aresame(fit.alpha,alpha,0.06)&&aresame(fit.lambda,h,0.002);
reporttest('MATERNFIT recovers Matern parameters with unit sample rate',allall(bool))

fit2=maternfit(dt,z,frac(1,2)*pi,'range.sigma',[0 17 100],'range.alpha',[0 1.5 100],'range.lambda',[1/1000 1/10 10]);

clear bool
bool(1)=aresame(fit.sigma,fit2.sigma,1e-2);
bool(2)=aresame(fit.alpha,fit2.alpha,1e-4);
bool(3)=aresame(fit.lambda,fit2.lambda,1e-4);

reporttest('MATERNFIT is independent of initial guess',allall(bool))

% use fit
% [fi,spp,snn]=maternspec(dt,length(z),sigma,alpha,lambda);
% [f,Spp,Snn]=mspec(dt,z,[]);
% figure,plot(f,[Spp spp]),hold on,plot(-f,[Snn snn]),ylog
     
rng(0);
dt=3600;
z=maternoise(dt,1000,sigo,alpha,h./dt);
fit=maternfit(dt,z,frac(1,2)*pi./dt);
bool=aresame(abs(fit.sigma./sigo-1),0,0.2)&&aresame(fit.alpha,alpha,0.06)&&aresame(fit.lambda,h./dt,0.002);
reporttest('MATERNFIT recovers Matern parameters with non-unit sample rate',allall(bool))

h=1;
rng(0);
dt=3600;

z=maternoise(dt,1000,sigo,alpha,h./dt);
psi=sleptap(size(z,1),3,1);
fit=maternfit(dt,z,frac(1,2)*pi./dt,'tapered',psi);
bool=aresame(abs(fit.sigma./sigo-1),0,0.2)&&aresame(fit.alpha,alpha,0.06)&&aresame(fit.lambda,h./dt,0.002);
reporttest('MATERNFIT recovers Matern parameters with non-unit sample rate, tapered version',allall(bool))

psi=sleptap(size(z,1)-1,3,1);
fit=maternfit(dt,z,frac(1,2)*pi./dt,'tapered',psi,'difference');
bool=aresame(abs(fit.sigma./sigo-1),0,0.2)&&aresame(fit.alpha,alpha,0.06)&&aresame(fit.lambda,h./dt,0.002);
reporttest('MATERNFIT recovers Matern parameters with non-unit sample rate, differenced tapered version',allall(bool))




