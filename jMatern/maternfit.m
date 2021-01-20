function [structout] = maternfit(varargin)
%MATERNFIT  Parametric spectral fit to the Matern form. [with A. Sykulski]
%
%   MATERNFIT performs a parametric fit of the spectrum of a time series to
%   that expected for a Matern process plus optional spin.  
%
%   The time series may either be real-valued or complex-valued.
%
%   The Matern process plus spin has a spectrum given by, see MATERNSPEC, 
%
%        SPP(F) = SIGMA^2 / [(F-NU)^2/LAMBDA^2 + 1]^ALPHA 
%                           / (LAMBDA * MATERNC(ALPHA)) + 2*pi * EPSILON^2
%
%   where SIGMA is the standard deviation, NU is a frequency shift, LAMBDA 
%   is a damping coefficient, ALPHA is 1/2 of the spectral slope, and 
%   EPSILON^2 is the variance of an optional additive noise component.
%
%   The coefficient LAMBDA * MATERNC(ALPHA) in the above form lets us 
%   parameterize the spectrum in terms of the process variance SIGMA^2.
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
%     Sykulski, Olhede, Guillaumin, Lilly, and Early (2019). The de-biased
%        Whittle likelihood. Biometrika, 106 (2): 251--266.
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
%   FIT=MATERNFIT(DT,Z,FO) where Z is a times series oriented as column 
%   vector, returns the result of a fitting the periodogram of Z to a
%   spectrum having a Matern form, fit over all frequencies. The range 
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
%                     Low  Guess  High 
%        SIGMA    =  [0      1    100] * STD(Z)
%        ALPHA    =  [1/2    1     10] 
%        LAMBDA   =  [1e-3  2e-3   10] * ABS(FO)
%        EPSILON  =  [0      0      0] * STD(Z)
%        NU       =  [0      0      0] * ABS(FO)
%
%   Note that by default, neither the EPSILON nor NU parameters are used. 
%
%   The search ranges and guesses can all be modified, as described below. 
%
%   A common special case is the inclusion of an additive noise component. 
%   Calling flag MATERNFIT(...'noisy'), which sets
% 
%        EPSILON  =  [0     0.01     2] * STD(Z)
%
%   for the range of the EPSILON parameter.  According to the spectral
%   normalizations used here, a white noice time series with standard 
%   deviation EPSILON will have a spectral  
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
%   Extensions
%
%   Two extensions of the basic Matern form are supported, discussed in 
%   more detail in MATERNSPEC.
% 
%   An oscillatory version, with nonzero MU, and a noisy version, with 
%   nonzero EPSILON, are available as options with these two extensions.
%
%   Generalized Matern
%
%   MATERFIT(...,'general') fits to the spectrum of generalized Matern
%   process. A parameter GAMMA is used with the default range
%
%        GAMMA  =  [1/2     1    20] 
%
%   such that the initial guess is the standard Matern process.  
%
%   The difference between the generalized Matern and standard Matern is
%   limited to frequency band in the vicinity of the falloff frequency, so 
%   unless you are averaging over many time series (decribed below), you
%   can expect a lot of variability in the fit value of GAMMA.
%
%   Extended Matern
%
%   MATERFIT(...,'extended') fits to the spectrum of an extended Matern
%   process. A parameter MU is used with the default range
%
%            MU  =  [1/100  10  100] * ABS(FO)
%
%   and in this case, the ALPHA parameter has the default range
%
%       ALPHA    =  [-1/2    1    10]
%
%   which has a different lower bound than for the standard Matern process.
%   __________________________________________________________________
%   
%   Options for multiple input time series
%
%   FIT=MATERNFIT(DT,Z,FO) may have Z being a matrix with N columns, or a 
%   cell array of N different time series.  In both of these cases, DT and 
%   FO may be scalars or length N arrays.
%
%   In these cases, the fields SIGMA, ALPHA, LAMBDA, and NU of FIT will 
%   also be arrays with N elements, as all non-string fields of PARAMS. The 
%   fields of RANGE will then all be N x 3 arrays.  
%
%   Alternatively, MATERNFIT(...,'average') with Z a matrix or a cell
%   array of matrices causes the columns of each matrix to be interpreted 
%   as members of an ensemble, averaging over columns to create an average 
%   spectrum. One fit per matrix is returned, rather than one per column. 
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
%   The ranges of LAMBDA and NU are nondimensional values representing
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
%     MATERNFIT(...,'con',...) alternately uses FMINCON with the default
%       interior-point algorithm.  This requires Matlab's Optimization 
%       Toolbox to be installed.  This is mainly for testing.
%
%     MATERNFIT(...,'nlo',...) uses the Nelder-Mead algorithm from the
%       NLopt toolbox at https://nlopt.readthedocs.io.  This requires NLopt
%       to be installed.  Again, this is mainly for testing at the moment.
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
%   array or a matrix, MATERNFIT(...,'parallel') will loop over the 
%   elements of Z using a parfor loop to speed things up.
%
%   This choice is reflected in the output field PARAMS.CORES, which 
%   will take on the values 'series' (the default) or 'parallel'.
%   __________________________________________________________________
%
%   'maternfit --t' runs a test.
%
%   Usage: fit=maternfit(dt,z,fo);
%          fit=maternfit(dt,z,fo,a,b);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2020 J.M. Lilly and A.M. Sykulski
%                                --- type 'help jlab_license' for details


%   Composite Matern
%
%   MATERFIT(...,'composite') fits to the spectrum of an composite Matern
%   process. A parameter MU is used with the default range

%Notes to self:
%Note, this is all set up to run the 5-parameter Matern process, I have
%just not checked the code.  Also make sure to comment, which I can grab 
%from SPECFIT, and also make sure to assign bgtype in output params.
%
%This is basically just a stripped-down version of SPECFIT.  Use Mac's
%Filemerge tool to compare them. 

if nargin>0
%    if strcmpi(varargin{1}, '--f')
 %       maternfit_fig,return
    if strcmpi(varargin{1}, '--t')
        maternfit_test,return
    end
end
%--------------------------------------------------------------------------
%Set search ranges for parameters
range=struct;

%Ranges for the background process
range.sigma=[0 1 100];            %Range of standard deviation as fraction of total
range.alpha=[1/2 1 10];           %Range of slope parameter
range.lambda=[1/1000 2/1000 10];  %Range of decay parameter as multiple of Coriolis
range.mu=[0 0 0];                 %Range of damping parameter as multiple of inverse Coriolis
range.nu=[0 0 0];                 %Frequency shift is fixed at NU=0
range.epsilon=[0 0.01 2];         %Noise level
%--------------------------------------------------------------------------
%Initialize value structure
names=fieldnames(range);
for i=1:length(names)
    eval(['value.' names{i} '=[];'])
end
%--------------------------------------------------------------------------
dt=double(varargin{1});
z=varargin{2};
varargin=varargin(3:end);
%--------------------------------------------------------------------------
%Sorting out input options 

opts.side='both';            %Determines side: opposite, same, or both
opts.alg='bnd';              %Algorithm: bnd, con, or nlopt
opts.ver='difference';       %Fit version: difference or raw
opts.cores='serial';         %Serial or parallel computation
opts.bgtype='mat';      %Type of model for background
opts.noise='clean';          %Noisy or clean
opts.cols='fit';             %Fit or average

psi=[];                      %The default data taper is the empty taper
verstrin=false;              %Flag for whether or not version string is input
for i=1:30
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'con')||strcmpi(varargin{end}(1:3),'bnd')||strcmpi(varargin{end}(1:3),'nlo')
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
        elseif strcmpi(varargin{end}(1:3),'mat')||strcmpi(varargin{end}(1:3),'exp')||strcmpi(varargin{end}(1:3),'ext')||strcmpi(varargin{end}(1:3),'gen')
            opts.bgtype=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'noi')||strcmpi(varargin{end}(1:3),'cle')
            opts.noise=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'fit')||strcmpi(varargin{end}(1:3),'ave')
            opts.cols=varargin{end};
        end
        varargin=varargin(1:end-1);
    elseif length(varargin)>1
        if ischar(varargin{end-1})&&~ischar(varargin{end})
            if strcmpi(varargin{end-1}(1:3),'tap')
                psi=varargin{end};
            elseif strcmpi(varargin{end-1}(1:3),'ran')
                name=varargin{end-1}(7:end);
                range=setfield(range,name,varargin{end});
            elseif strcmpi(varargin{end-1}(1:3),'val')
                name=varargin{end-1}(7:end);
                value=setfield(value,name,varargin{end});
            end
            varargin=varargin(1:end-2);
        end
    end
end
%--------------------------------------------------------------------------
%Make sure ranges and flags are set correctly for specified model type
if strcmpi(opts.bgtype(1:3),'exp')
    range.alpha=[-1/2 -1/2 -1/2];
    range.mu=[1/100 10 100];  
    %range.lambda=[1/100000 2/1000 10];  %Range of decay parameter as multiple of Coriolis
elseif strcmpi(opts.bgtype(1:3),'ext')
    range.alpha=[-1/2 1 10];
    range.mu=[1/100 10 100];  
elseif strcmpi(opts.bgtype(1:3),'gen')
    range.mu=[1/2 1 20];  
elseif strcmp(opts.bgtype(1:3),'mat')
    if range.alpha(1)<=1/2
        %Correction for alpha being less than permitted value
        range.alpha(1)=1/2+1e-10;
    end
end
if strcmpi(opts.noise(1:3),'cle')
    range.epsilon=[0 0 0];
end
if isfield(range,'gamma')
    range.mu=range.gamma;
    range=rmfield(range,'gamma');
end
%--------------------------------------------------------------------------
%Convert value structure to a numeric array
if iscell(z)
    N=length(z);
else
    N=size(z,2);
end
values=nan*zeros(6,3);
values(1,:)=range.sigma;
values(2,:)=range.alpha;
values(3,:)=range.lambda;
values(4,:)=range.mu;
values(5,:)=range.nu;
values(6,:)=range.epsilon;
if ~isempty(value.sigma),   values(1,:)=value.sigma;end
if ~isempty(value.alpha),   values(1,:)=value.alpha;end
if ~isempty(value.lambda),  values(1,:)=value.lambda;end
if ~isempty(value.mu),      values(1,:)=value.mu;end
if ~isempty(value.nu),      values(1,:)=value.nu;end
if ~isempty(value.epsilon), values(1,:)=value.epsilon;end
%Have to convert this to numeric values because I can't pass a structure
%through to the workers, for some unknown reason
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
     elseif ~strcmpi(opts.cols(1:3),'ave')
        dt=dt+zeros(size(z,2),1);
    end
end

fc=fc(:);
if length(fc)==1
    if iscell(z)
        fc=fc+zeros(length(z),1);
    elseif ~strcmpi(opts.cols(1:3),'ave')
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
if ~iscell(z)&&( size(z,2)==1 || strcmpi(opts.cols(1:3),'ave'))
    %Single time series input, no loop
    [x,xa,xb,xo,like,aicc,P,a,b,err,exitflag,iters]=...
        maternfit_one(dt,z,flow,fhigh,opts,fc,psi,win,stdz,values);
    %The rest of this simply implements different types of loops
elseif ~iscell(z)
    %Loop over matrix columns
    N=size(z,2);[x,xa,xb,xo]=vzeros(N,6);
    [like,aicc,P,a,b,err,exitflag,iters]=vzeros(N,1);
    if strcmpi(opts.cores(1:3),'ser')
        %For series loop over matrix columns
        for i=1:size(z,2)
            disp(['MATERNFIT processing series ' int2str(i) ' of ' int2str(size(z,2)) '.'])
            [x(i,:),xa(i,:),xb(i,:),xo(i,:),like(i),aicc(i),P(i),a(i),b(i),err(i),exitflag(i),iters(i)]=...
                maternfit_one(dt(i),z(:,i),flow,fhigh,opts,fc(i),psi,win,stdz(i),values);
        end
    elseif strcmpi(opts.cores(1:3),'par')
        %For parallel loop over matrix columns
        disp('MATERNFIT using parallel for loop...')
        parfor i=1:size(z,2)
            %dt(i),flow,fhigh,opts,fc(i),fit
            disp(['MATERNFIT processing series ' int2str(i) ' of ' int2str(size(z,2)) '.'])
            [x(i,:),xa(i,:),xb(i,:),xo(i,:),like(i),aicc(i),P(i),a(i),b(i),err(i),exitflag(i),iters(i)]=...
                maternfit_one(dt(i),z(:,i),flow,fhigh,opts,fc(i),psi,win,stdz(i),values);
        end
    end
elseif iscell(z)
    %Loop over cell array
    N=length(z);[x,xa,xb,xo]=vzeros(N,6);
    [like,aicc,P,a,b,err,exitflag,iters]=vzeros(N,1);
    if strcmpi(opts.cores(1:3),'ser')
        %For series loop over cell array
        for i=1:length(z)
            disp(['MATERNFIT processing series ' int2str(i) ' of ' int2str(length(z)) '.'])
            [x(i,:),xa(i,:),xb(i,:),xo(i,:),like(i),aicc(i),P(i),a(i),b(i),err(i),exitflag(i),iters(i)]=...
                maternfit_one(dt(i),z{i},flow,fhigh,opts,fc(i),psi{i},win{i},stdz{i},values);
        end
    elseif strcmpi(opts.cores(1:3),'par')
        %For parallel loop over cell array
        disp('MATERNFIT using parallel for loop...')
        parfor i=1:length(z) 
            disp(['MATERNFIT processing series ' int2str(i) ' of ' int2str(length(z)) '.'])
            [x(i,:),xa(i,:),xb(i,:),xo(i,:),like(i),aicc(i),P(i),a(i),b(i),err(i),exitflag(i),iters(i)]=...
                maternfit_one(dt(i),z{i},flow,fhigh,opts,fc(i),psi{i},win{i},stdz{i},values);
        end
    end
end

sigma=x(:,1);
alpha=x(:,2);
lambda=x(:,3);
mu=x(:,4);
nu=x(:,5);
epsilon=x(:,6);

clear range
range.sigma=  [xa(:,1)  xo(:,1)  xb(:,1)];
range.alpha=  [xa(:,2)  xo(:,2)  xb(:,2)];
range.lambda= [xa(:,3)  xo(:,3)  xb(:,3)];
if strcmpi(opts.bgtype(1:3),'gen')
    range.gamma=[xa(:,4)  xo(:,4)  xb(:,4)];
else
    range.mu=[xa(:,4)  xo(:,4)  xb(:,4)];
end
range.nu=     [xa(:,5)  xo(:,5)  xb(:,5)];
range.epsilon=[xa(:,6)  xo(:,6)  xb(:,6)];

structout=[];

use opts

make params dt fc a b P like aicc err exitflag iters side alg ver cores 
%make params dt fc a b P like aicc err exitflag iters side alg ver bgtype cores 

if strcmpi(opts.bgtype(1:3),'gen')
    gamma=mu;
    make structout sigma alpha lambda gamma nu epsilon range params
else
    make structout sigma alpha lambda mu nu epsilon range params
end

%remove some fields I don't use
if strcmpi(opts.bgtype(1:3),'mat')
    structout=rmfield(structout,'mu');
    structout.range=rmfield(structout.range,'mu');
end
if strcmpi(opts.noise(1:3),'cle')
    structout=rmfield(structout,'epsilon');
    structout.range=rmfield(structout.range,'epsilon');
end
if all(range.nu==0)
    structout=rmfield(structout,'nu');
    structout.range=rmfield(structout.range,'nu');
end



%--------------------------------------------------------------------------
function[x,exitflag,iters]=maternfit_optimize(dt,N,om,spp,snn,index,xa,xo,xb,opts,realflag)
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
    options=optimset('MaxFunEvals',10000,'MaxIter',10000,'TolFun',tol,'TolX',tol);
    [xf,fval1,exitflag,struct]=fminsearchbnd(@(z) specmodel(dt,z.*xo,spp,snn,index,N,opts,realflag),...
        guess,xanorm,xbnorm,options);
    iters=struct.iterations;
elseif strcmpi(opts.alg(1:3),'con')
    options=optimoptions('fmincon','Algorithm','interior-point','MaxFunEvals',10000,'MaxIter',10000,'TolFun',tol,'TolX',tol,'Display','off');
    [xf,fval1,exitflag,struct]=fmincon(@(z) specmodel(dt,z.*xo,spp,snn,index,N,opts,realflag),...
        guess,[],[],[],[],xanorm,xbnorm,[],options);
    iters=struct.iterations;
elseif strcmpi(opts.alg(1:3),'nlo')
    %opt.algorithm = NLOPT_GN_DIRECT_L;
    opt.algorithm = NLOPT_LN_NELDERMEAD;  %Works, not too much slower than Matlab
    %opt.algorithm = NLOPT_LN_PRAXIS;   %Works,but slower than NM
    %opt.algorithm = NLOPT_GN_CRS2_LM;
    %opt.algorithm = NLOPT_GD_STOGO;
    %opt.algorithm = NLOPT_GN_ISRES;
    %opt.algorithm = NLOPT_GN_ESCH;
    %opt.algorithm = NLOPT_LN_COBYLA;
    %opt.algorithm =NLOPT_LN_SBPLX;  %Works,but slower than NM
    %opt.algorithm = NLOPT_LN_BOBYQA;  %Somewhat slower than NM, and different results
    %opt.algorithm = NLOPT_AUGLAG;    opt.lower_bounds = xanorm;
    opt.lower_bounds = xanorm;
    opt.upper_bounds = xbnorm;
    opt.min_objective = @(z) specmodel(dt,z.*xo,spp,snn,index,N,opts,realflag);
    opt.fc_tol = [tol tol tol tol];
    opt.xtol_rel = tol;
    [xf, fmin, exitflag] = nlopt_optimize(opt, guess);
    iters=nan;  %Not sure how to return this
end

x=xf.*xo;% scale back to correct units
[like,Spp,Snn,f]=specmodel(dt,x,spp,snn,index,N,opts,realflag);
%--------------------------------------------------------------------------
function [like,Spp,Snn,f,indexn]=specmodel(dt,x,spphat,snnhat,index,N,opts,realflag)
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
R=maternfit_materncov(dt,N,x,realflag,opts.bgtype);    
[f,Spp,Snn]=maternfit_blurspec(dt,R,opts.ver,opts.win);

%length(R)
%[f,Spp2,Snn2]=blurspec(1,R,opts.ver,'window',opts.win);aresame(Spp,Spp2)
%Same answers but slower

%length(index),N,length(f),length(Spp)
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
function[x,xa,xb,xo,like,aicc,P,a,b,err,exitflag,iters]=maternfit_one(dt,z,flow,fhigh,opts,fc,psi,win,stdz,values)

%This is pretty dumb, but Matlab doesn't want me to pass a structure 
%through parfor.  Thus I have to convert to a matrix and then back again.
sigma=values(1,:);
alpha=values(2,:);
lambda=values(3,:);
mu=values(4,:);
nu=values(5,:);
epsilon=values(6,:);

[x,xa,xb,xo]=vzeros(1,6,nan);
[like,aicc,P,a,b,err,exitflag,iters]=vzeros(1,1,nan);
if anyany(~isfinite(z))
    return
end

N=size(z,1);
realflag=isreal(z);

%Put window into options to simplify argument passing
opts.win=win;

if strcmpi(opts.cols(1:3),'ave')
    stdz=sqrt(vmean(squared(stdz),2));
end

%Multiply standard deviation ranges by sample standard deviation 
sigma    = sigma*stdz;
epsilon  = epsilon*stdz;

%Multiply damping and frequency ranges by Coriolis frequency
lambda  = lambda*abs(fc);   
if ~strcmpi(opts.bgtype(1:3),'gen')
    mu      = mu.*abs(fc);   %Is this a time or frequency scale?
end
xa=[sigma(1)  alpha(1)   lambda(1)  mu(1)  nu(1)  epsilon(1)];
xo=[sigma(2)  alpha(2)   lambda(2)  mu(2)  nu(2)  epsilon(2)];
xb=[sigma(3)  alpha(3)   lambda(3)  mu(3)  nu(3)  epsilon(3)];

%Periodogram from mspec  ... use dt = 1
%vsize(dt,z,psi),length(find(isfinite(z))),isreal(z)
[om,spp,snn]=mspec(dt,z,psi);
if isreal(z)
    snn=spp;
end

if strcmpi(opts.cols(1:3),'ave')
    spp=vmean(spp,2);
    snn=vmean(snn,2);
end

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
%vsize(dt,N,om,spp,snn,a:b,xa,xo,xb,opts,realflag)
%[xa;xo;xb]
[x,exitflag,iters]=maternfit_optimize(dt,N,om,spp,snn,a:b,xa,xo,xb,opts,realflag);

%Final value of likelihood and spectra
%vsize(dt,x,spp,snn,a:b,N,opts,realflag)
[like,Spp,Snn,f,indexn]=specmodel(dt,x,spp,snn,a:b,N,opts,realflag);

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
aicc=2*like + frac(4*M*P,2*M-P-1);
    
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
%Nyquist is also repeated for even. See comments at MSPEC.

if iseven(N)
    if index(end)==((N/2)+1)
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
R(1,:)=R(1,:)./2;    %Don't forget to divide first element by two
S=dt*2*real(fft(R)); %But I do take it into account in spectrum

S=abs(S);  %Sometimes there are small negative parts after blurring
N=size(R,1);
omega=frac(1,dt)*2*pi*(0:floor(N/2))'./N;
%omega=fourier(N);

Spp=S(1:length(omega),:);
Snn=[S(1,:);S(end:-1:end-length(omega)+2,:)];
%--------------------------------------------------------------------------
function[R]=maternfit_materncov(dt,N,x,realflag,model)
sigma=x(1);
alpha=x(2);
lambda=x(3);
mu=x(4);
nu=x(5);
epsilon=x(6);

tau=dt*[0:N-1]';

if isnan(sigma)
    R=[];
else
    %This is just copied from MATERNCOV, but it is faster to have it internal
    if strcmpi(model(1:3),'mat')
        fact=2*frac(1,gamma(alpha-1/2).*pow2(alpha-1/2));
        R=fact.*((lambda*tau).^(alpha-1/2)).*besselk(abs(alpha-1/2),lambda*tau);
        R(1)=1;%because of being undefiend there
    elseif strcmpi(model(1:3),'ext')
        if alpha==-1/2
            tnorm=sqrt(tau.^2+mu.^2);
            fact=besselk(1,mu.*lambda);
            R=frac(1,fact).*frac(mu,tnorm).*besselk(1,lambda.*tnorm);
        else
            tnorm=lambda.*sqrt(tau.^2+mu.^2);
            fact=1./(abs(mu.*lambda)).^(alpha-1/2)./besselk(abs(alpha-1/2),abs(mu.*lambda));
            R=fact*tnorm.^(alpha-1/2).*besselk(abs(alpha-1/2),tnorm);
        end
    elseif strcmpi(model(1:3),'gen')
        M=10;P=10;  %Specifying oversampling rates for numerical computations
        [f,Spp,Snn]=maternspec(dt,M*N*P,sigma,alpha,lambda/P,mu,'generalized');
        S=[flipud(Snn(2:end));Spp];%plot(Spp),hold on
        Ri=ifft(ifftshift(S))./dt;  %Make sure it's ifftshift not fftshift
        Ri=Ri(1:P:end);
        R=Ri(1:N);%plot(Ri),hold on
        R=R./R(1);
    elseif strcmpi(model(1:3),'com')
        M=10;P=10;  %Specifying oversampling rates for numerical computations
        [f,Spp,Snn]=maternspec(dt,M*N*P,sigma,alpha,lambda/P,mu/P,nu/P,'composite');
        S=[flipud(Snn(2:end));Spp];%figure,plot(S)
        Ri=ifft(ifftshift(S))./dt;  %Make sure it's ifftshift not fftshift
        Ri=Ri(1:P:end);
        R=Ri(1:N);
        R=R./R(1);
    end
    if (nu~=0)&&~strcmpi(model(1:3),'com')
        R=R.*exp(sqrt(-1)*tau*nu);
    end
    R=R.*sigma.^2;
    R(1)=R(1)+squared(epsilon);
end

if realflag
     R=real(R);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of function body; begin tests and figures
function[]=maternfit_test
%--------------------------------------------------------------------------
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

rng(0);
epsilono=2;
zn=epsilono.*(randn(length(z),1)+1i.*randn(length(z),1))./sqrt(2);
fit=maternfit(dt,zn,frac(1,2)*pi,'noisy','value.sigma',0);
bool=(abs(fit.sigma)<1e-7)&&aresame(abs(fit.epsilon./epsilono-1),0,0.05);
reporttest('MATERNFIT recovers Matern parameters with unit sample rate, noise only',allall(bool))

rng(0);
fit=maternfit(dt,z+zn,frac(1,2)*pi,'noisy');
bool=aresame(abs(fit.sigma./sigo-1),0,0.2)&&aresame(fit.alpha,alpha,0.5)&&aresame(fit.lambda,h,0.1)&&aresame(abs(fit.epsilon./epsilono-1),0,0.3);
reporttest('MATERNFIT recovers Matern parameters with unit sample rate, noisy version',allall(bool))
%--------------------------------------------------------------------------
sigo=17;
alpha=1.5;
h=1/10;

rng(0);
dt=1;
z=maternoise(dt,1000,sigo,alpha,h,'real');
fit=maternfit(dt,z,frac(1,2)*pi);
bool=aresame(abs(fit.sigma./sigo-1),0,0.2)&&aresame(fit.alpha,alpha,0.06)&&aresame(fit.lambda,h,0.0025);
reporttest('MATERNFIT recovers Matern parameters with unit sample rate and real signal',allall(bool))

fit2=maternfit(dt,z,frac(1,2)*pi,'range.sigma',[0 17 100],'range.alpha',[0 1.5 100],'range.lambda',[1/1000 1/10 10]);

clear bool
bool(1)=aresame(fit.sigma,fit2.sigma,1e-2);
bool(2)=aresame(fit.alpha,fit2.alpha,1e-4);
bool(3)=aresame(fit.lambda,fit2.lambda,1e-4);

reporttest('MATERNFIT is independent of initial guess for real signal',allall(bool))

rng(0);
epsilono=2;
zn=epsilono.*randn(length(z),1);
fit=maternfit(dt,zn,frac(1,2)*pi,'noisy','value.sigma',0);
bool=(abs(fit.sigma)<1e-7)&&aresame(abs(fit.epsilon./epsilono-1),0,0.1);
reporttest('MATERNFIT recovers Matern parameters with unit sample rate and real signal, noise only',allall(bool))

rng(0);
fit=maternfit(dt,z+zn,frac(1,2)*pi,'noisy');
bool=aresame(abs(fit.sigma./sigo-1),0,0.2)&&aresame(fit.alpha,alpha,0.5)&&aresame(fit.lambda,h,0.1)&&aresame(abs(fit.epsilon./epsilono-1),0,0.3);
reporttest('MATERNFIT recovers Matern parameters with unit sample rate and real signal, noisy version',allall(bool))
%--------------------------------------------------------------------------

% rng(0);
% dt=3600;
% z=maternoise(dt,1000,sigo,alpha,h./dt);
% fit=maternfit(dt,real(z),frac(1,2)*pi./dt);
% bool=aresame(abs(fit.sigma./sigo-1),0,0.2)&&aresame(fit.alpha,alpha,0.06)&&aresame(fit.lambda,h./dt,0.002);
% reporttest('MATERNFIT recovers Matern parameters with non-unit sample rate',allall(bool))

rng(0);
dt=3600;
z=maternoise(dt,1000,sigo,alpha,h./dt);
fit=maternfit(dt,z,frac(1,2)*pi./dt);
bool=aresame(abs(fit.sigma./sigo-1),0,0.2)&&aresame(fit.alpha,alpha,0.06)&&aresame(fit.lambda,h./dt,0.002);
reporttest('MATERNFIT recovers Matern parameters with non-unit sample rate',allall(bool))

rng(0);
epsilono=2;
zn=epsilono.*(randn(length(z),1)+1i.*randn(length(z),1))./sqrt(2);
fit=maternfit(dt,zn,frac(1,2)*pi,'noisy','value.sigma',0);
bool=(abs(fit.sigma)<1e-7)&&aresame(abs(fit.epsilon./epsilono-1),0,0.1);
reporttest('MATERNFIT recovers Matern parameters with non-unit sample rate, noise only',allall(bool))

h=1;
rng(0);
dt=3600;

z=maternoise(dt,1000,sigo,alpha,h./dt);
psi=sleptap(size(z,1),3,1);
fit=maternfit(dt,z,frac(1,2)*pi./dt,'tapered',psi);
bool=aresame(abs(fit.sigma./sigo-1),0,0.2)&&aresame(fit.alpha,alpha,0.06)&&aresame(fit.lambda,h./dt,0.002);
reporttest('MATERNFIT recovers Matern parameters with non-unit sample rate, tapered version',allall(bool))

psi=sleptap(size(z,1)-1,3,1);
tic;
fit=maternfit(dt,z,frac(1,2)*pi./dt,'tapered',psi,'difference');
bool=aresame(abs(fit.sigma./sigo-1),0,0.2)&&aresame(fit.alpha,alpha,0.06)&&aresame(fit.lambda,h./dt,0.002);
reporttest('MATERNFIT recovers Matern parameters with non-unit sample rate, differenced tapered version',allall(bool))
etime1=toc;

sigo=17;
alpha=1.5;
h=1/10;
rng(0);

dt=3600;
gamma=5;
N=1000;
% clf
% z=maternoise(dt,[N,2000],sigo,alpha,h./dt,1,'generalized');
% plot(vmean(abs(fft(z)),2)),xlog,ylog
% z=maternoise(dt,[N,2000],sigo,alpha,h./dt,'chol');
% hold on,plot(vmean(abs(fft(z)),2)),xlog,ylog
z=maternoise(dt,[N,5000],sigo,alpha,h./dt,5,'generalized');
zn=epsilono.*randn(size(z,1),5000);
%hold on,plot(vmean(abs(fft(z)),2))

% [tau,tpz1]=materncov(dt,N,sigo,alpha,h./dt);
% [tau,tpz2]=materncov(dt,N,sigo,alpha,h./dt,1,'generalized');
% [tau,tpz3]=materncov(dt,N,sigo,alpha,h./dt,2,'generalized');
% [tau,tpz4]=materncov(dt,N,sigo,alpha,h./dt,5,'generalized');
% figure,plot(tau,[tpz1 tpz2 tpz3 tpz4])

%[f,s]=maternspec(dt,N,sigo,alpha,h./dt,5,'generalized');plot(f,s./maxmax(s))
%[f,s]=maternspec(dt,N,sigo,alpha,h./dt,1,'generalized');hold on,plot(f,s./maxmax(s))
%[f,s]=maternspec(dt,N,sigo,alpha,h./dt,1/2,'generalized');hold on,plot(f,s./maxmax(s))
%[f,s]=maternspec(dt,N,sigo,alpha,h./dt,1/10,'generalized');hold on,plot(f,s./maxmax(s))

fit=maternfit(dt,z,frac(1,2)*pi./dt,'generalized','average');
%fit=maternfit(dt,z,frac(1,2)*pi./dt,'generalized');
bool=aresame(abs(fit.sigma./sigo-1),0,0.2)&&aresame(fit.alpha,alpha,0.06)&&aresame(fit.lambda,h./dt,0.002)&&aresame(fit.gamma,5,0.25);
reporttest('MATERNFIT recovers generalized Matern parameters with non-unit sample rate',allall(bool))

fit=maternfit(dt,z+zn,frac(1,2)*pi./dt,'generalized','average','noisy');
%fit=maternfit(dt,z,frac(1,2)*pi./dt,'generalized');
bool=aresame(abs(fit.sigma./sigo-1),0,0.2)&&aresame(fit.alpha,alpha,0.06)&&aresame(fit.lambda,h./dt,0.002)&&aresame(fit.gamma,5,0.25)&&aresame(abs(fit.epsilon./epsilono-1),0,0.1);
reporttest('MATERNFIT recovers generalized Matern parameters with non-unit sample rate, noisy version',allall(bool))

% if exist('nlopt_optimize')==3
%     psi=sleptap(size(z,1)-1,3,1);
%     tic;
%     fit=maternfit(dt,z,frac(1,2)*pi./dt,'tapered',psi,'difference','nlopt');
%     bool=aresame(abs(fit.sigma./sigo-1),0,0.2)&&aresame(fit.alpha,alpha,0.06)&&aresame(fit.lambda,h./dt,0.002);
%     reporttest('MATERNFIT using NLopt version of previous',allall(bool))
%     etime2=toc;
%     disp(['MATERNFIT NLopt version took ' num2str(etime2./etime1) ' as much time as FMINSEARCH.'])
% end


% sigo=17;
% alpha=1.5;
% h=1/10;
% 
% rng(0);
% dt=1;
% z=maternoise(dt,1000,sigo,alpha,h);
% fit=maternfit(dt,[z z],frac(1,2)*pi);
