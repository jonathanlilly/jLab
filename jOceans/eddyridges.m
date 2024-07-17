function[struct]=eddyridges(varargin)
%EDDYRIDGES  Coherent eddy ridges from Lagrangian trajectories.
%
%   EDDYRIDGES extracts eddy-like displacement signals from a float or 
%   drifter position record using wavelet ridge analysis. 
%
%   By default, EDDYRIDGES analyzes Lagrangian trajectories on the earth, 
%   or from an earth-like model.  For analsis of trajectories in Cartesian
%   model domains, see the 'Cartesian domains' section below.
%
%   STRUCT=EDDYRIDGES(NUM,LAT,LON,FMAX,FMIN,P,M,RHO) runs the ellipse
%   extraction algorithm and returns the results in a data structure,
%   described below.
%
%   The input arguments are as follows:
%
%      NUM     Date in Matlab's DATENUM format 
%      LAT     Latitude record of float
%      LON     Longitude record of float
%      FMAX    Maximum ratio of frequency to Coriolis frequency
%      FMIN    Minimum ratio of frequency to Coriolis frequency
%      P       Wavelet duration, P>=1 
%      M       Ridge length cutoff in units of wavelet length 2P/pi
%      RHO     Ridge trimming factor, removing RHO*P/pi from each end
%
%   NUM, LAT, and LON are all column vectors of the same length, or cell 
%   arrays of column vectors, while the other arguments are all scalars.  
%
%   FMAX, FMIN, P, M, and RHO are parameters controlling the ridge
%   analysis, described in more detail next.
%   _______________________________________________________________________
%
%   Analysis parameters
%
%   EDDYRIDGES works by calling RIDGEWALK to determine time/frequency
%   curves called "ridges" that are composed of special points, "ridge
%   points", that identify modulated oscillations such as the signal of a
%   coherent eddy.  RIDGEWALK in turn is based on a wavelet transform. 
%
%   FMAX and FMIN define the frequency band for the ridge analysis, P 
%   controls the wavelet used, and L and R are length and amplitude cutoffs
%   applied during the ridge analysis to minimize spurious features.
%
%   Frequency range: FMAX and FMIN
%
%   FMAX and FMIN define a moving band of frequencies with respect to the
%   local Coriolis frequency.  Both of these are positive numbers, with 1
%   being the Corilios frequency.  Ridge point lying outside this band are
%   rejected.  A typical choice of eddy band would be FMIN=1/32, FMAX=1/2.
%
%   Wavelet duration: P
%
%   The wavelet duration P controls the degree of time localization.  The 
%   number of oscillations contained within the wavelet is roughly 2P/pi,
%   so frequency resolution increases as P increases.  P should be greater
%   than about SQRT(3). 
%
%   Ridge length: M
%
%   M sets the minimum length for a ridge, in units of the approximate 
%   number of oscillations spanned by the wavelet, 2P/pi.  Thus the ridge
%   must contain about M times as many oscillations as the wavelet.  The 
%   ridges become more statistically significant as M increases.
%
%   Ridge trimming: RHO
%
%   After the ridge length criterion is applied, RHO*(2/pi) oscillations,
%   corresponding to one wavelet half-width, are removed from each end of
%   the ridge, as these are generally contaminated by edge effects.  The 
%   choice RHO=1/2 is recommended.  Note the minium ridge length will then 
%   be (M-RHO)*2P/pi. See RIDGETRIM for details.  
%   _______________________________________________________________________
%
%   Output 
%  
%   The following parameter are output as fields of the structure STRUCT. 
%   Type 'use struct' (using the actual name of the structure) to put these
%   fields into named variables in the workspace. 
%
%      -------------------------------------------------------------------     
%      Signal fields:
%
%      NUM       Date in Matlab's DATENUM format along each ridge
%      LAT       Latitude along each ridge
%      LON       Longitude along each ridge
%      LATRES    Latitude residual after subtracting eddy signal
%      LONRES    Longitude residual after subtracting eddy signal
%      ZHAT      Estimated eddy displacement XHAT+iYHAT, in kilometers     
%      -------------------------------------------------------------------
%      Ellipse fields:
%
%      KAPPA     Ellipse amplitude SQRT((A^2+B^2)/2), in kilometers 
%      XI        Ellipse circularity 2AB/(A^2+B^2), nondimensional
%      THETA     Ellipse orientation, in radians
%      PHI       Ellipse orbital phase, in radians
%      OMEGA     Ellipse instantaneous frequency, in radians/day 
%      UPSILON   Ellipse instantaneous bandwidth, in radians/day
%      CHI       Ellipse normalized instantaneous curvature, nondimensional
%      R         Ellipse geometric mean radius, in kilometers
%      V         Ellipse kinetic energy velocity, in cm/s
%      RO        Ellipse Rossby number, nondimensional
%      -------------------------------------------------------------------
%      Ridge fields:
%
%      IR      Indices into rows of wavelet transform (time) 
%      JR      Indices into columns of wavelet transform (scale) 
%      KR      Indices into data columns or cells (time series number)
%      LEN      Ridge length in number of complete periods
%      -------------------------------------------------------------------
%      PARAMS  Sub-structure of internal parameters used in the analysis
%  
%   PARAMS is a structure that contains the parameter values of FMAX, FMIN,
%   P, M, RHO, GAMMA, BETA, and FS used in the analysis.  BETA is an array
%   of length LENGTH(NUM), FS is a cell array of that length, and the other
%   elements of PARAMS are all scalars. 
%
%   Apart from PARAMS, all other output fields are cell arrays of column 
%   with one ridge per cell.  
%
%   Note that XI and V are signed quantities that are positive for
%   counterclockwise rotation and negative for counterclockwise rotation, 
%   while RO is a signed quantity that is positive for cyclonic and 
%   negative for anticyclonic rotations.
% 
%   If more than one time series is input, the ridge field KR provides an
%   index into the original column. To find all cells associated with the
%   5th time series, for instance, use the index FIND(CELLFIRST(KR)==5).
%
%   If you want the ridges to be separated for each input signal, call 
%   EDDYRIDGES with an external loop.
%   _______________________________________________________________________
%
%   Cartesian domains
%
%   EDDYRIDGES can also be used to analyze trajectories from idealized 
%   model domains, which are better represented in Cartesian coordinates
%   rather than latitude and longitude. 
%
%   STRUCT=EDDYRIDGES(F,NUM,Z,FMAX,FMIN,P,M,RHO) performs the eddy 
%   extraction on data with positions given in Cartesian coordinates.
%
%   The first three input arguments have changed from the planetary case:
%
%      F      Coriolis frequency of the model domain, in radians per second
%      NUM    Time in days along each float trajectory
%      Z      Complex-valued float location X+iY, in kilometers
%
%   F must be a scalar, with LENGTH(F)=1, in order for EDDYRIDGES to know
%   that the Cartesian algorithm is intended.  F is used to calculate the
%   frequency range for the ridge analysis, and for the eddy Rossby number.  
%   
%   Note that for definiteness, EDDYRIDGES requires that the model fields 
%   are expressed in particular physical units.
%
%   EDDYRIDGES(F,FBETA,NUM,...) accounts for variations of the Coriolis
%   frequency due to the beta effect in calculating the ridge frequency
%   range.  FBETA is a scalar with units of radians per meter per second.
%
%   The output arguments are then all the same, *except* for the signal 
%   fields, which become
%
%      -------------------------------------------------------------------     
%      Signal fields:
%      NUM     Time in days along each ridge
%      Z       Complex-valued float location X+iY, in kilometers
%      ZRES    Residual after subtracting eddy XRES+iYRES, in kilometers
%      ZHAT    Estimated eddy displacement XHAT+iYHAT, in kilometers     
%      -------------------------------------------------------------------
% 
%   Only the signal fields changed.  All of the ellipse fields and ridge
%   fields are the same as in the planetary case.  
%  
%   The values of F and BETA used are included in the sub-structure PARAMS.
%   _______________________________________________________________________
%
%   Algorithm variations
%   
%   By default, EDDYRIDGES uses a split-sides approach to finding the
%   ridges, in which ridges dominated by positive rotations (XI>1) and
%   those dominated by negative rotations (XI<1) are found separately
%   and then combined; this is done using RIDGEWALK's 'mask' functionality.
%
%   Under this approach, ridges are not permitted to change their sign of
%   rotation.  This is physically realistic as particle paths in eddies are
%   expected to be close to circular, with some distortion expected due
%   ambient strain and other processes.  Eddies do not suddenly change
%   their rotation sense from one sign to the other.  
%
%   This algorithm has the effect of attenuating the detection of events 
%   with wandering polarizations that are observed to occur regularly in 
%   noise, but that are unlikely to correspond to eddies.
%
%   For details, see 
%
%       Lilly and Perez-Brunius (2021b), Extracting statistically 
%           significant eddy signals from large Lagrangian datasets using
%           wavelet ridge analysis, with application to the Gulf of Mexico.
% 
%   EDDYRIDGES(...,'nosplit') suppresses the use of this split-sides 
%   algorithm, such that any polarization senses are permitted.
%   _______________________________________________________________________
%
%   Outputting the transforms
%
%   Often it can be useful to look at the actual wavelet transforms, even
%   if one does not want to output these when analyzing large datasets.  
%
%   Therefore, if the input fields NUM, LAT, and LON contain only a single
%   time series, the rotary transforms WP and WN computed internally by
%   EDDYRIDGES are returned as the last two sub-fields of PARAMS.
%   _______________________________________________________________________
%
%   Parallelization
%
%   EDDYRIDGES(...,'parallel') parellelizes the computation by using a 
%   PARFOR when looping over multiple trajectories.
%   _______________________________________________________________________
%
%   'eddyridges --t' runs a test.
%
%   See also RIDGEWALK.
%
%   Usage: struct=eddyridges(num,lat,lon,1/2,1/64,3,0,1);
%          struct=eddyridges(num,lat,lon,fmax,fmin,P,M,rho);
%          struct=eddyridges(f,num,z,fmax,fmin,P,M,rho);
%          struct=eddyridges(f,beta,num,z,fmax,fmin,P,M,rho);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2023 J.M. Lilly --- type 'help jlab_license' for details

if strcmp(varargin{1}, '--t')
    eddyridges_test,return
end


%Sort out input arguments
ga=3;
D=16;
alpha=0.05;
alg='sph';
fo=[];
beta=0;

parstr='series';
splitstr='split';
for i=1:2
    if isstr(varargin{end})
        if strcmpi(varargin{end}(1:3),'ser')||strcmpi(varargin{end}(1:3),'par')
            parstr=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'spl')||strcmpi(varargin{end}(1:3),'nos')
            splitstr=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end


if ~iscell(varargin{1})&&length(varargin{1})==1
    fo=varargin{1};
    alg='car';
    varargin=varargin(2:end);
    if ~iscell(varargin{1})&&length(varargin{1})==1
        beta=varargin{1};
        varargin=varargin(2:end);
    end
end

num=varargin{1};
if strcmp(alg,'car')  %Swap lat and lon, since lat=x and lon=y
    lat=imag(varargin{2});  %This is Y in km
    lon=real(varargin{2});  %This is X in km
    varargin=varargin(3:end);
else
    lat=varargin{2};
    lon=varargin{3};
    varargin=varargin(4:end);
end
fmax=varargin{1};
fmin=varargin{2};
P=varargin{3};
M=varargin{4};
rho=varargin{5};

if strcmpi(parstr(1:3),'ser')
    struct=eddyridges_one_series(ga,D,num,lat,lon,fmax,fmin,P,M,rho,alpha,fo,beta,alg,splitstr);
elseif strcmpi(parstr(1:3),'par')
    struct=eddyridges_one_parallel(ga,D,num,lat,lon,fmax,fmin,P,M,rho,alpha,fo,beta,alg,splitstr);
end

% if length(struct.num)==1
%     names=fieldnames(struct);
%     for i=1:length(names)-1
%         eval(['struct.' names{i} '=struct.' names{i} '{1};'])
%     end
% end


%lofotenridges=catstruct(struct);


%--------------------------------------------------------------------------
%Series version
function [struct]=eddyridges_one_series(ga,D,numo,lato,lono,fmax,fmin,P,M,rho,alpha,fo,beta,alg,splitstr)

%If it isn't in cell form, put in cell form
if ~iscell(lato)
    if size(numo,2)==1
        numo=vrep(numo,size(lato,2),2);
    end
    for i=1:size(numo,2)
        numcell{i,1}=numo(:,i);
        latcell{i,1}=lato(:,i);
        loncell{i,1}=lono(:,i);
    end
    numo=numcell;lato=latcell;lono=loncell;
    clear numcell latcell loncell
end


[latrefo,ir,jr,kr,xr,yr,chi,num,lat,lon,latref,...
    kappa,lambda,theta,phi,len,V,R,a,b,c,fs,...
    omega,upsilon,zhat,Ro,z,zres,latres,lonres]=initializecell(length(numo),1);

be=zeros(length(numo),1);

for i=1:length(numo)
    disp(['EDDYRIDGES working on time series ' int2str(i) ' of ' int2str(length(numo)) '.'])
    
    if strcmp(alg,'car')
        %This is the Cartesian case
        fcor=abs(fo+beta*1000*lato{i})*3600*24;  %lato is just y in km here ... rad per day
        latrefo{i}=asind(fcor./24./3600./(2*7.292e-5));  %For future reference
    else
        %This is the planetary case
        fcor=abs(corfreq(lato{i}))*24;  %Local Coriolis frequency in radians per day
        latrefo{i}=lato{i};  %For future reference
    end
    
    %vsize(numo{i},lato{i},lono{i},alpha,P(i),L(i),fmax,fmin,ga,D,C(i),fcor)
    if length(numo)==1  %get wavelet transforms as output
        [ir{i},jr{i},xr{i},yr{i},omega{i},chi{i},be(i),fs{i},V{i},R{i},kappa{i},lambda{i},theta{i},phi{i},zhat{i},wp,wn]=...
            eddyridges_loop(numo{i},lato{i},lono{i},alpha,P,M,rho,fmax,fmin,ga,D,fcor,alg,splitstr);
    else
        [ir{i},jr{i},xr{i},yr{i},omega{i},chi{i},be(i),fs{i},V{i},R{i},kappa{i},lambda{i},theta{i},phi{i},zhat{i}]=...
            eddyridges_loop(numo{i},lato{i},lono{i},alpha,P,M,rho,fmax,fmin,ga,D,fcor,alg,splitstr);
    end
    if ~isempty(ir{i})
        dt=numo{i}(2)-numo{i}(1);
        bool=(~isnan(ir{i}));
        
        num{i}=ir{i};
        lat{i}=ir{i};
        lon{i}=ir{i};
        latref{i}=ir{i};
        kr{i}=i+0*ir{i};
        
        num{i}(bool)=numo{i}(ir{i}(bool));
        lat{i}(bool)=lato{i}(ir{i}(bool));
        lon{i}(bool)=lono{i}(ir{i}(bool));
        latref{i}(bool)=latrefo{i}(ir{i}(bool));
   
        len{i}=ridgelen(dt,omega{i});

        %V{i}=ellvel(dt,kappa{i},lambda{i},theta{i},phi{i},1e5/24/3600,'kin');
        %R{i}=ellrad(kappa{i},lambda{i});
        [a{i},b{i},c{i},upsilon{i}]=ellband(dt,kappa{i},lambda{i},theta{i},phi{i});
        
        %Note I'm using the reference latitude here, computed above
        Ro{i}=ellrossby(latref{i},sign(lambda{i}).*sqrt(1-squared(lambda{i})),omega{i});
        
        if strcmp(alg,'car')
            z{i}=lon{i}+1i*lat{i};
            zres{i}=z{i}-zhat{i};
        else
            [latres{i},lonres{i}]=sig2latlon(xr{i},yr{i},lat{i},lon{i});
        end
    end
end

disp('EDDYRIDGES done with ridge calculation, rearranging ridges...')

make params fmax fmin P M rho
params.gamma=ga;
params.beta=be;
params.fs=fs;

%convert linearity to circularity
xi=col2cell(sign(cell2col(lambda)).*sqrt(1-squared(cell2col(lambda))));

if strcmp(alg,'car')
    %the third line enforces the shortest length, after multiplicity adjustment
    cell2col(num,z,zres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi,'nonans');
    col2cell(num,z,zres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi);
    %vindex(num,z,zres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi,cellfirst(len)>(M-rho)*2*P/pi,1);
    cellprune(num,z,zres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi);

    fbeta=beta;
    make params fo fbeta
    %make struct num z zres zhat kappa xi theta phi omega R V Ro upsilon ir jr kr chi len params
    make struct num z zres zhat kappa xi theta phi omega upsilon chi R V Ro ir jr kr len params
else  
    %the third line enforces the shortest length, after multiplicity adjustment
    cell2col(num,lat,lon,latres,lonres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi,'nonans');
    col2cell(num,lat,lon,latres,lonres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi);    
    %vindex(num,lat,lon,latres,lonres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi,cellfirst(len)>(M-rho)*2*P/pi,1);
    cellprune(num,lat,lon,latres,lonres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi);
    
%    make struct num lat lon latres lonres zhat kappa xi theta phi omega R V Ro upsilon ir jr kr chi len params
    make struct num lat lon latres lonres zhat kappa xi theta phi omega upsilon chi R V Ro ir jr kr len params
end

if length(numo)==1
    fs=fs{1};
    make params fs wp wn
    make struct params
end

disp('EDDYRIDGES finished.')

%--------------------------------------------------------------------------
%Parallel version... exactly repeated code except for the parfor
function [struct]=eddyridges_one_parallel(ga,D,numo,lato,lono,fmax,fmin,P,M,rho,alpha,fo,beta,alg,splitstr)

%If it isn't in cell form, put in cell form
if ~iscell(lato)
    if size(numo,2)==1
        numo=vrep(numo,size(lato,2),2);
    end
    for i=1:size(numo,2)
        numcell{i,1}=numo(:,i);
        latcell{i,1}=lato(:,i);
        loncell{i,1}=lono(:,i);
    end
    numo=numcell;lato=latcell;lono=loncell;
    clear numcell latcell loncell
end

N=cell(length(numo),1);

[latrefo,ir,jr,kr,xr,yr,chi,num,lat,lon,latref,...
    kappa,lambda,theta,phi,len,V,R,a,b,c,fs,...
    omega,upsilon,zhat,Ro,z,zres,latres,lonres]=initializecell(length(numo),1);

be=zeros(length(numo),1);

parfor i=1:length(numo)
     disp(['EDDYRIDGES working on time series ' int2str(i) ' of ' int2str(length(numo)) '.'])
    
    if strcmp(alg(1:3),'car')
        %This is the Cartesian case
        fcor=abs(fo+beta*1000*lato{i})*3600*24;  %lato is just y in km here ... rad per day
        latrefo{i}=asind(fcor./24./3600./(2*7.292e-5));  %For future reference
    else
        %This is the planetary case
        fcor=abs(corfreq(lato{i}))*24;  %Local Coriolis frequency in radians per day
        latrefo{i}=lato{i};  %For future reference
    end
    
    %vsize(numo{i},lato{i},lono{i},alpha,P(i),L(i),fmax,fmin,ga,D,C(i),fcor)
    
    if length(numo)==1  %get wavelet transforms as output
        [ir{i},jr{i},xr{i},yr{i},omega{i},chi{i},be(i),fs{i},V{i},R{i},kappa{i},lambda{i},theta{i},phi{i},zhat{i},wp,wn]=...
            eddyridges_loop(numo{i},lato{i},lono{i},alpha,P,M,rho,fmax,fmin,ga,D,fcor,alg,splitstr);
    else
        [ir{i},jr{i},xr{i},yr{i},omega{i},chi{i},be(i),fs{i},V{i},R{i},kappa{i},lambda{i},theta{i},phi{i},zhat{i}]=...
            eddyridges_loop(numo{i},lato{i},lono{i},alpha,P,M,rho,fmax,fmin,ga,D,fcor,alg,splitstr);
    end
    
    if ~isempty(ir{i})
        dt=numo{i}(2)-numo{i}(1);
        bool=(~isnan(ir{i}));
        
        num{i}=ir{i};
        lat{i}=ir{i};
        lon{i}=ir{i};
        latref{i}=ir{i};
        kr{i}=i+0*ir{i};
        
        num{i}(bool)=numo{i}(ir{i}(bool));
        lat{i}(bool)=lato{i}(ir{i}(bool));
        lon{i}(bool)=lono{i}(ir{i}(bool));
        latref{i}(bool)=latrefo{i}(ir{i}(bool));
      
        %[kappa{i},lambda{i},theta{i},phi{i}]=ellparams(xr{i},yr{i});
        len{i}=ridgelen(dt,omega{i});
        %V{i}=ellvel(dt,kappa{i},lambda{i},theta{i},phi{i},1e5/24/3600,'kin');
        %R{i}=ellrad(kappa{i},lambda{i});
        [a{i},b{i},c{i},upsilon{i}]=ellband(dt,kappa{i},lambda{i},theta{i},phi{i});
        
        %Note I'm using the reference latitude here, computed above
        Ro{i}=ellrossby(latref{i},sign(lambda{i}).*sqrt(1-squared(lambda{i})),omega{i});
        
        if strcmp(alg,'car')
            z{i}=lon{i}+1i*lat{i};
            zres{i}=z{i}-zhat{i};
        else
            [latres{i},lonres{i}]=sig2latlon(xr{i},yr{i},lat{i},lon{i});
        end
    end
end

%vsize(kappa,nc,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi)

disp('EDDYRIDGES done with ridge calculation, rearranging ridges...')

make params fmax fmin P M rho
params.gamma=ga;
params.beta=be;
params.fs=fs;

%convert linearity to circularity
xi=col2cell(sign(cell2col(lambda)).*sqrt(1-squared(cell2col(lambda))));

if strcmp(alg,'car')
    %the third line enforces the shortest length, after multiplicity adjustment
    %vsize(num,z,zres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi)
    %num,z,zres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi
    cell2col(num,z,zres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi,'nonans');
    col2cell(num,z,zres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi);
    %vindex(num,z,zres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi,cellfirst(len)>(M-rho)*2*P/pi,1);
    cellprune(num,z,zres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi);

    make params fo beta
    make struct num z zres zhat kappa xi theta phi omega upsilon chi R V Ro ir jr kr len params
else
    %the third line enforces the shortest length, after multiplicity adjustment
    cell2col(num,lat,lon,latres,lonres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi,'nonans');
    col2cell(num,lat,lon,latres,lonres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi);    
    %vindex(num,lat,lon,latres,lonres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi,cellfirst(len)>(M-rho)*2*P/pi,1);
    cellprune(num,lat,lon,latres,lonres,zhat,kappa,xi,theta,phi,omega,R,V,Ro,upsilon,ir,jr,kr,len,chi);

    make struct num lat lon latres lonres zhat kappa xi theta phi omega upsilon chi R V Ro ir jr kr len params
end


if length(numo)==1
    fs=fs{1};
    make params fs wp wn
    make struct params
end

disp('EDDYRIDGES finished.')

function[ir,jr,xr,yr,omega,chi,be,fs,Vr,Rr,kappar,lambdar,thetar,phir,zhat,wp,wn]=...
    eddyridges_loop(num,lat,lon,alpha,P,M,rho,rmax,rmin,ga,D,fcor,alg,splitstr)

[ir,jr,xr,yr,omega,chi,fs,Vr,Rr,kappar,lambdar,thetar,phir,zhat]=vempty;
be=squared(P)./ga;

if length(num)>=1
    dt=num(2)-num(1);        %Days per sample point
    fmax=maxmax(fcor)*rmax;  %Radians per day
    fmin=minmin(fcor)*rmin;  %Radians per day
    
    fs=morsespace(ga,be,{alpha,fmax*dt},max(fmin*dt,frac(P,2).*2*pi./length(num)),D);  %Radians per sample point
    %[psi,psif]=morsewave(length(num),ga,be,fs);
    %figure,uvplot(psi(:,1))
    
    if  ~isempty(fs)
        if strcmpi(alg,'car')
            %Note that lon = x and lat = y herein
            [wx,wy]=wavetrans(lon,lat,{ga,be,fs,'bandpass'},'mirror');
        else
            [wx,wy]=spheretrans(lat,lon,{ga,be,fs,'bandpass'},'mirror');
        end
        
        wp=wx+sqrt(-1)*wy;%Positive rotations
        wn=wx-sqrt(-1)*wy;%Negative rotations
            
        %w(:,:,1)=wx;
        %w(:,:,2)=wy;
        %figure,jpcolor(log(vmean(squared(w),3))'),hold on
        %bool=isridgepoint(w,fs,0,'amplitude',[],[],[]);
        %[ii,jj]=find(bool);
        
        %plot(ii,jj,'m.')
        %figure,plot(fs)1
        %figure,plot([lon lat])
        %figure,jpcolor(num,fs,(sqrt(squared(wx')+squared(wy')))),shading interp
        
        %Not using this anymore
        %snro=log10(frac(squared(wp),squared(wn)));
        %[wx,wy,fs]=wavetransderiv(ga,be,wx,wy,fs);
        %figure,subplot(2,1,1),jpcolor(abs(wx)'),flipy,subplot(2,1,2),jpcolor(abs(wy)'),flipy
        if strcmpi(splitstr(1:3),'nos')
            %Version not splitting ridges
            [xr,yr,ir,jr,omega,chi]=ridgewalk(dt,wx,wy,fs,P,M,rho,[rmax*fcor,rmin*fcor]);
        else            
            %Splitting positive and negative ridges ... really reduces artifacts
            %         %[P,M,rho,rmax,rmin]
            [xrp,yrp,irp,jrp,omegap,chip]=ridgewalk(dt,wx,wy,fs,P,M,rho,[rmax*fcor,rmin*fcor],'mask',abs(wp)>abs(wn));
            [xrn,yrn,irn,jrn,omegan,chin]=ridgewalk(dt,wx,wy,fs,P,M,rho,[rmax*fcor,rmin*fcor],'mask',abs(wp)<abs(wn));
            
            %figure,plot(irp,jrp,'bo'),hold on,plot(irn,jrn,'r')
            if isempty(xrn)
                xr=xrp;yr=yrp;ir=irp;jr=jrp;omega=omegap;chi=chip;
            elseif isempty(xrp)&&~isempty(xrn)
                xr=xrn;yr=yrn;ir=irn;jr=jrn;omega=omegan;chi=chin;
            elseif ~isempty(xrp)&&~isempty(xrn)
                xr=[xrp;nan+1i*nan;xrn];
                yr=[yrp;nan+1i*nan;yrn];
                ir=[irp;nan;irn];
                jr=[jrp;nan;jrn];
                omega=[omegap;nan;omegan];
                chi=[chip;nan;chin];
            end
        end
        %------------------------------------------------------------------
        %This section computes radius and velocity with interpolation,
        %removing some of the small-scale jaggedness of velocity
        clear w
        w(:,:,1)=wx;
        w(:,:,2)=wy;
        [~,rq]=isridgepoint(w,fs,0,'amp',[],[],[]);
        [kappa,lambda,theta,phi]=ellparams(wx,wy);
        V=ellvel(dt,kappa,lambda,theta,phi,1e5/24/3600,'kin');
        Vr=vzeros(length(ir),1,'nan');
        bool=isfinite(ir);
        Vr(bool)=ridgeinterp(rq,ir(bool),jr(bool),V);
        [kappar,lambdar,thetar,phir]=ellparams(xr,yr);
        Rr=ellrad(kappar,lambdar);
        %------------------------------------------------------------------
        zhat=real(xr)+1i*real(yr);
%         %sort out multiplicity...
%         %keep whichever ridge explains more velocity variance over the overlap
%         iro=min(ir);
%         mult=vzeros(max(ir)-iro+1,1);%length of number of index points
%         irc=col2cell(ir);
%         for i=1:length(irc)
%             mult(irc{i}-iro+1)=mult(irc{i}-iro+1)+1;%adding one
%         end
%         
%         multr=ir;
%         multr(~isnan(ir))=mult(ir(~isnan(ir))-iro+1);%putting back along ridge
%         
%         
%         %figure,plot(ir,jr),hold on,plot(ir,multr)
%         
%         if anyany(mult>1)
%             omegamat=nan*zeros(max(ir)-iro+1,length(find(isnan(ir))));
%             %cell2col(ir,jr,xr,yr,omega,chi,Vr,Rr,kappar,lambdar,thetar,phir,multr,zhat,len);
%             omegac=col2cell(omega);
%             for i=1:length(irc)
%                 omegamat(irc{i}-iro+1,i)=omegac{i};
%             end
%             %figure,plot(omegamat)
%             %figure,plot(ir,omega,'.')
%             
%             omegamaxmat=vrep(max(omegamat,[],2),size(omegamat,2),2);
%             boolmat=(omegamat==omegamaxmat);
%             %figure,plot(boolmat)
%             boolrc=irc;
%             for i=1:size(boolmat,2)
%                 boolrc{i}=boolmat(irc{i}-iro+1,i);
%             end
%             bool=cell2col(boolrc);
%             %vsize(ir,bool)
%             %aresame(find(isnan(bool)),find(isnan(ir)))
%             
%             % vsize(ir,jr,xr,yr,omega,chi,Vr,Rr,kappar,lambdar,thetar,phir,multr)
%              
%             ir(bool==0)=nan; 
% %             len=ridgelen(dt,omega);
% %             
% %             lenr=nan*zeros(length(mult),length(find(isnan(len))));
% %             [irc,lenc]=col2cell(ir,len);
% %             for i=1:length(irc)
% %                   lenr(irc{i}-iro+1,i)=lenc{i};
% %             end
% %             figure,plot(lenr)
% %             
% %             col2mat(ir,jr,xr,yr,omega,chi,Vr,Rr,kappar,lambdar,thetar,phir,multr,zhat,len);
% %             lenr=ir;
% %             lenr(~isnan(ir))=len(ir(~isnan(ir)));%put 
% %             
% %             
% %             [~,ii]=max(len(1,:));
% %             %figure,plot(len)
%                       
%          
% %older verison, using eke criterion
% if 0
%     col2mat(ir,jr,xr,yr,omega,chi,Vr,Rr,kappar,lambdar,thetar,phir,multr,zhat);
%     
%     if strcmp(alg,'car')
%         cv=vdiff(lon+1i*lat,1)./dt;
%     else
%         cv=latlon2uv(num,lat,lon);
%     end
%     cvr=ir;
%     cvr(~isnan(ir))=cv(ir(~isnan(ir)));
%     
%     %vsize(ir,cvr,irc,numer,multr)
%     numer=squared(cvr-vdiff(zhat,1)./dt);
%     denom=squared(cvr);
%     
%     %figure,plot(numer./denom)
%     numer(multr<=1)=nan;
%     denom(multr<=1)=nan;
%     
%     err=vsum(numer,1)./vsum(denom,1);
%     [~,ii]=min(err);
%     
%     
%     %set ridges to infs during periods when multiplicity exceeds one
%     %except for that ridge which minimizes the error
%     %(or has the longest length)
%     for i=1:size(ir,2)
%         if max(multr(:,i))>1
%             if i~=ii
%                 ir(multr(:,i)>1,i)=inf;
%             end
%         end
%     end
% end
% 
% %             remove columns with all nans
%             %bool=sum(isnan(ir),1)~=size(ir,1)
%             %size(ir)
%             %vindex(ir,jr,xr,yr,omega,chi,Vr,Rr,kappar,lambdar,thetar,phir,zhat,bool,2);
%             %size(ir)
%             %       find(isnan(ir(1,:)))
%             %find(isnan(jr(1,:)))
%             %mat2col(ir,jr,xr,yr,omega,chi,Vr,Rr,kappar,lambdar,thetar,phir,zhat);
%             col2cell(ir,jr,xr,yr,omega,chi,Vr,Rr,kappar,lambdar,thetar,phir,zhat);
%             cellpack(ir,jr,xr,yr,omega,chi,Vr,Rr,kappar,lambdar,thetar,phir,zhat);
%             cell2col(ir,jr,xr,yr,omega,chi,Vr,Rr,kappar,lambdar,thetar,phir,zhat);
%             %find(isnan(ir(1,:)))
%             %find(isnan(jr(1,:)))
%             %figure,plot(ir,jr)
%         end
        
        %figure,plot(mr)
        %figure,plot(ir,mr),hold on,plot(ir,jr)
        
        %[a{i},b{i},c{i},upsilon{i}]=ellband(dt,kappa{i},lambda{i},theta{i},phi{i});
        %zhat{i}=real(xr{i})+1i*real(yr{i});
        
        %Note I'm using the reference latitude here, computed above
        %Ro{i}=ellrossby(latref{i},lambda{i},omega{i});
        
        %vsize(xr,yr,ir,jr,omega,err)
        %length(find(isnan(xr)))
        
        %[xr(end,:) yr(end,:) ir(end,:) jr(end,:) omega(end,:) err(end,:)]
        
        %Determine multiplicity
        %[~,multo]=ridgemap(length(num),xr,ir);
        %bool=(~isnan(ir));
        %mult=ir;
        %mult(bool)=multo(ir(bool));
    end
end

function[varargout]=initializecell(M,N)

for i=1:nargout
    varargout{i}=cell(M,N);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin tests

function[]=eddyridges_test

eddyridges_test1;


function[]=eddyridges_test1

N=1000;
num=[1:N]'/24;
phi=4*num;

[x,y,z]=ellsig(100,0,0,phi,0,0,'real');
z=sqrt(radearth^2-x.^2-y.^2);
[x,y,z]=vectmult(jmat3(20*2*pi/360,1),x,y,z);

[lat,lon]=xyz2latlon(x,y,z);
struct=eddyridges(num,lat,lon,1/2,1/64,3,2,1);

ii=25:925;
use struct
cell2col(lat,lon,latres,lonres,kappa);
vindex(lat,lon,latres,lonres,kappa,ii,1);

tol=5e-2;
bool=aresame(median(latres),70,tol)&aresame(median(lonres),-90,tol)&aresame(median(kappa),100,1/2);
reporttest('EDDYRIDGES latres and lonres match expected value for 1000 km radius circle',allall(bool))


%testing Cartesian vs. spherical algorithm
rng(0);
dt=0.25;
N=1000; 
num=[0:N-1]'*dt;
z=randn(N,1)+1i*randn(N,1);
z=z-mean(z);z=z./vstd(z,1)*(1/100000);%make very small amplitude
pos=cumsum(z).*(dt*24*3600)/1000/100;
tic;testeddies1=eddyridges(corfreq(30)/3600,num,pos,2,1/64,sqrt(6),1,0);time1=toc;
[lat,lon]=uv2latlon(num,real(z),imag(z),30,0);
%plot(pos),hold on,plot(latlon2xy(lat,lon,30,0))
tic;testeddies2=eddyridges(num,lat,lon,2,1/64,sqrt(6),1,0);time2=toc;
bool=aresame(length(testeddies1.num),length(testeddies2.num));
bool=bool&maxmax(abs(cellength(testeddies1.num)-cellength(testeddies2.num)))<=1;
reporttest('EDDYRIDGES spherical and Cartesian algorithms match',bool)



% N=1000;
% num=[1:N]'/24;
% phi=4*num;
% 
% [x,y,z]=ellsig(100,0,0,phi,0,0,'real');
% z=sqrt(radearth^2-x.^2-y.^2);
% [x,y,z]=vectmult(jmat3(20*2*pi/360,1),x,y,z);
% 
% [lat,lon]=xyz2latlon(x,y,z);
% [lat,lon]=xyz2latlon(x,y,z);
% lat(400:600)=0;lon(400:600)=0;
% %lat=lat+randn(size(lat));
% %lon=lon+randn(size(lon));
% struct=eddyridges([num num],[lat lat],[lon lon],1/2,1/64,3,2,1);



% function[]=eddyridges_test2
% 
% load drifterbetty
% use drifterbetty
% 
% struct1=eddyridges(num,lat,lon,1/2,1/64,3,log10(3),4);
% 
% 
% 
% 
% cx=latlon2xy(lat,lon,35,-70);
% struct2=eddyridges(abs(corfreq(35))/3600,num,cx,1/2,1/64,3,log10(3),4);
% 
% %Okay!  These are very close...  just write the test code
% %Other tests... Test works fine for negative latitudes
% use drifterbetty
% struct3=eddyridges(num,-lat,lon,1/2,1/64,3,log10(3),4);
% 
% cx=latlon2xy(-lat,lon,-35,-70);
% struct4=eddyridges(abs(corfreq(35))/3600,num,cx,1/2,1/64,3,log10(3),4);
% 
% 
% 
% function[]=eddyridges_test6
% load ebasnfloats
% use ebasnfloats
% struct=eddyridges(num,lat,lon,1/2,1/64,3,2,1);
% 
% 
% 
% 
% struct=eddyridges(num,lat,lon,1/2,1/64,3,log10(3),4);
% 
% 
% 
% function[]=eddyridges_test7
% load npg2006
% use npg2006
% 
% struct=eddyridges(num,lat,lon,1/2,1/64,3,log10(3),4);
% %Looks great!!  
% 
% struct2=eddyridges(num,lat,lon,3,1,1/2,1/64,'2d');
% struct3=eddyridges(num,lat,lon,3,1,1/2,1/64,'3d');
% 
% use struct2.ridges
% struct2.ridges.xhat=ellsig(kappa,lambda,theta,phi);
% use struct3.ridges
% struct3.ridges.xhat=ellsig(kappa,lambda,theta,phi);
% 
% 
% clear bool
% bool(1)=aresame(struct2.ridges.kappa,struct3.ridges.kappa,0.08);
% bool(2)=aresame(struct2.ridges.lambda,struct3.ridges.lambda,1e-2);
% bool(3)=allall(vmedian(abs(rot(2*struct2.ridges.theta)-rot(2*struct3.ridges.theta)),1)<1e-2);
% bool(4)=aresame(struct2.ridges.xhat,struct3.ridges.xhat,.2);
% bool(5)=aresame(struct2.data.latres,struct3.data.latres,1e-3);
% bool(6)=aresame(struct2.data.lonres,struct3.data.lonres,2e-3);
% reporttest('EDDYRIDGES 3D vs. 2D algorithm for npg2006 data',allall(bool))
% 
% function[]=eddyridges_test3
% 
% load sitelfloats
% use sitelfloats
% 
% num=num{4};
% lat=lat{4};
% lon=lon{4};
% 
% struct=eddyridges(num,lat,lon,1/2,1/64,3,log10(3),2);
% 
% %does it help if I double length?
% num2=[num(1):1/2:num(end)-1/2]';
% 
% struct=eddyridges(num2,doublen(lat(1:end-1)),doublen(lon(1:end-1)),1/2,1/64,3,log10(3),2);
% struct=eddyridges(num2,doublen(lat(1:end-1)),doublen(lon(1:end-1)),1/2,1/64,sqrt(3),log10(3),2);
% 
% %sqrt(3) seems better here
% 
% 
% %Decide on frequencies
% fs=morsespace(3,3,{.2,pi},2*pi/500,8);
% 
% cx=latlon2xy(lat,lon,36,-62.5);
% %Compute wavelet transforms using generalized Morse wavelets
% [wx,wy]=wavetrans(real(cx),imag(cx),{1,3,3,fs,'bandpass'},'mirror');
% [ir,jr,wr,fr]=ridgewalk(wx,wy,fs,{1.5,0,'amp'});
% [ir2,jr2,wr2,fr2]=ridgewalk(wx,wy,fs,{1.5,0,'amp'},{ga,be,1});
% 
% 
% 
% 
% 
% figure,jimage(rms(wx,wy)'),flipy
% plot(ir,jr,'m')
% 
% 
% cv=latlon2uv(num,lat,lon);
% %Compute wavelet transforms using generalized Morse wavelets
% [wu,wv]=wavetrans(real(cv),imag(cv),{1,3,3,fs,'bandpass'},'mirror');
% [iru,jru,wru,fru]=ridgewalk(wu,wv,fs,{1.5,0,'amp'});
% figure,jimage(rms(wu,wv)'),flipy
% plot(iru,jru,'m'),plot(ir,jr,'k')
% 
% 
% [iru,jru,wru,fru]=ridgewalk(wu,wv,fs,{5,1,'amp'});
% 
% [cvhatu,f]=ridgemap(length(cv),wru,fru,iru,'collapse');
% [cxhat,f]=ridgemap(length(cv),wr,fr,ir,'collapse');
% 
% 
% 
% %Testing wavetransderiv
% ga=3;be=3;
% fact=-sqrt(-1)*frac(morsea(ga,be+1),morsea(ga,be))*frac(morsefreq(ga,be),fs');
% fact=vrep(fact,size(wu,1),1);
% wx2=fact.*wu;
% wy2=fact.*wv;
% 
% 
% w=wu;
% w(:,:,2)=wv;
% rq=ridgequantity(w,fs,'amplitude');
% [wxr2,wyr2]=ridgeinterp(fs,rq,ir,jr,wx2,wy2);
% wr2=wxr2;
% wr2(:,2)=wyr2;
% [cxhat2,f]=ridgemap(length(cv),wr2,fr,ir,'collapse');
% 
% 
% figure,jimage(num,fs,(abs(fact).*rms(wu,wv))')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Just playing around here
% load ebasnfloats
% use ebasnfloats
% jj=36;num=num{jj};lat=lat{jj};lon=lon{jj};
% figure,plot(lon,lat);
% ga=3;be=3;
% fs=morsespace(ga,be,1/10,1/1000,4);  %Radians per sample point
% [wx,wy]=spheretrans(lat,lon,{ga,be,fs,'bandpass'},'mirror');
% 
% 
% %Question... is polarization distribution of a ridge the same as that of
% %noise?
% N=100000;
% z=randn(N,1)+1i*randn(N,1);
% 
% fs=morsespace(3,3,1/10,1/1000,4);
% [wp,wn]=wavetrans(z,conj(z),{3,3,fs,'amplitude'});
% 
% [ir,jr,wpr,wnr,fpr,fnr]=ridgewalk(1,wp,wn,fs,{1,0});
% 
% fbar=fpr;
% for i=1:length(fpr)
%     fbar{i}=frac(squared(wpr{i}).*fpr{i}+squared(wnr{i}).*fnr{i},squared(wpr{i})+squared(wnr{i}));
% end
% 
% 
% len1=ridgelen(dt,fbar);
% len=zeros(length(ir),1);
% for i=1:length(len)
%     len(i)=len1{i}(1);
% end
% 
% figure,
% [n,x]=hist(len,[0:1:100]/5);
% plot(x,n),hold on
% 
% %[n,x]=hist(len,100);
% %plot(x,n),hold on
% 
% 
% %Does red noise make longer ridges?
% N=100000;
% z=cumsum(randn(N,1)+1i*randn(N,1));
% 
% fs=morsespace(3,3,1/10,1/1000,4);
% [wp,wn]=wavetrans(z,conj(z),{3,3,fs,'amplitude'});
% 
% [ir,jr,wpr,wnr,fpr,fnr]=ridgewalk(1,wp,wn,fs,{1,0});
% 
% fbar=fpr;
% for i=1:length(fpr)
%     fbar{i}=frac(squared(wpr{i}).*fpr{i}+squared(wnr{i}).*fnr{i},squared(wpr{i})+squared(wnr{i}));
% end
% 
% len1=ridgelen(dt,fbar);
% len=zeros(length(ir),1);
% for i=1:length(len)
%     len(i)=len1{i}(1);
% end
% 
% [n,x]=hist(len,[0:1:100]/5);
% plot(x,n),hold on
% 
% 
% %Yes, but not greatly.  Not an order of magnitude
% %So, ridge length is strongly dependent upon slope.  This makes sense.
% 
% %What is distribution of mean power like?
% 
% [pbar,nbar]=vzeros(length(ir),1);
% for i=1:length(ir)
%     pbar(i)=sqrt(mean(squared(wpr{i})));
%     nbar(i)=sqrt(mean(squared(wnr{i})));
% end
% chi=log(squared(frac(pbar,nbar)));
% 
% 
% %Very interesting, longer ridges are less likely to have a particular
% %polarization value... longer ridges tend to be more isotropic
% plot(chi,len,'.')
% 
% %So, the double space of length vs. polarization appears the most
% %promising. Amplitude does not matter
% 
% %Could do: 
% %1) Matern fit to each time series.... 
% %2) Generate say 100 or 1000 fake time series of this length
% %3) Solve, time series by time series, for chi and len cutoffs....
% %  ... or just len?   Well, try this with a 
% 
% 
% %Does red noise make longer ridges?
% N=1000000;
% z=(randn(N,1)+1i*randn(N,1));
% [zp,zn]=anatrans(z,conj(z));
% chi=log(frac(squared(zp),squared(zn)));
% chi=log(frac(vfilt(squared(zp),[1 1 1 1 1]),vfilt(squared(zn),[1 1 1 1 1])));
% [n,x]=hist(chi,1000);
% plot(x,n)
% 
% 
% chi=log(frac(vfilt(squared(zp),[1 1 1 1 1]),vfilt(squared(zn),[1 1 1 1 1])));
% hold on 
% 
% [n,x]=hist(chi,1000);
% plot(x,n)
% %How can I do this without anything fancy
% %1) assume that +/- are independent Gaussians, average n squares
% %2) ratio of squares ==> chi^2_n/chi_^2_n = F_n,n
% %3) So now I have distribution of polarization *given* length ...
% %4) This avoids asking what the chance of something a certain length occuring is
% %So, can easily code up the f-dist and solve for confidence intervals...
% %5)  *Test* that noise ridges follow this distribution... 
% 
% This is very easy
% 
% 
%  


