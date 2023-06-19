function[zhat,beta,aux,res]=polymap(varargin)
%POLYMAP  Mapping using local polynomial fitting, also known as loess.  
%
%   POLYMAP generates a map from scattered data in two dimensions using
%   a locally weighted least squares fit to a polynomial.
%
%   This method is variously known as local polynomial fitting, local
%   polynomial smoothing, multivariate locally weighted least squares 
%   regression, lowess (originally for LOcally WEighted Scatterplot
%   Smoothing), and loess.  All of these are essentially synonyms.
%
%   POLYMAP has the support for all of the following options:
%
%       --- Cartesian or spherical geometry
%       --- a constant, linear, or quadratic fit in space
%       --- an additional linear or quadratic fit in time
%       --- fixed bandwidth or fixed population algorithms
%       --- a broad family of possible choices for the kernel
%       --- spatially-varying bandwidth, population, or kernel
%       --- additional datapoint weighting factors, e.g. by confidence
%       --- parallelization using Matlab's Parallel Computing toolbox 
%       --- the median-based robustification of Cleveland (1979) XXX
%
%   POLYMAP is implemented using a numerically efficient algorithm that
%   avoids an explicit loop over individual grid points.  Instead, there is
%   a parallelizable loop over the y-axis, or latitude axis, of the grid.
% 
%   For problems involving repeated mapping on the same measurement grid, 
%   the sub-components of POLYMAP can be called individually for a great 
%   computational savings, as described in "One grid, many fields" below.
%
%   For algorithm details, see XXX.
%   __________________________________________________________________
%
%   Local polynomial fitting on the plane
% 
%   Let's say we have an array Z of observations at locations X,Y. X,Y,
%   and Z can be arrays of any size provided they are all the same size.
%
%   The problem is to obtain a mapped field ZHAT on some regular grid 
%   specified by the vectors XO and YO.  These are arrays of length M and L
%   giving the bin center locations of an L x M matrix of grid points, i.e.
%
%      XO= [XO_1 XO_2 ... XO_M]      YO =  [YO_1;    
%                                           YO_2; 
%                                            ...
%                                           YO_L].
%
%   ZHAT=POLYMAP(X,Y,[],[],Z,XO,YO,{P,H,ALPHA,BETA}) performs a Pth order 
%   local polynomial fit to the data at each grid point, with the weighting
%   controlled by H, ALPHA, and BETA. ZHAT will be a matrix of size L x M.
% 
%   The fit is found by minimizing the weighted mean squared error between 
%   the fitted surface and the observed values.  
%
%   The data locations (X,Y) and grid point locations (X0,Y0) should have
%   the same units as the bandwidth H (e.g., kilometers).  Note that any
%   non-finite values in X, Y, or Z are ignored.
% 
%   The fit order P may be chosen as P=0 (fit to a constant), P=1 (fit to a
%   plane), or else P=2 (fit to a parabolic surface). 
%
%   The bandwidth H sets the distance outside at which the weight assigned
%   to the data will vanish.  It may either be a scalar or an L x M matrix.
%
%   ALPHA and BETA set the choice of weighting function, known as the 
%   kernel, as described next. 
%   __________________________________________________________________
%
%   Choice of kernel
%
%   Let R be the radial distance between an observation point and a grid 
%   point.  POLYMAP uses a kernel of the form
%
%        K = KAPPA * (1-(R/H)^ALPHA)^BETA
%
%   which is a general type of kernel called the generalized beta kernel.
%   KAPPA is a normalization coefficient that depends on ALPHA and BETA.
%
%   This form includes many kernels commonly used in the literature:
% 
%        ALPHA=1, BETA=0 --- uniform kernel, K=constant
%        ALPHA=2, BETA=1 --- Epanechnikov (parabolic) kernel, K~1-(DS/H)^2
%        ALPHA=2, BETA=2 --- bisquare kernel, K~(1-(DS/H)^2)^2
%        ALPHA=3, BETA=3 --- tricubic kernel, K~(1-(DS/H)^3)^3
%
%   Note that the kernel is defined to vanish for distances R>=H.
%
%   Like the bandwidth H, ALPHA and BETA may either be scalars, or matrices 
%   of size L x M to specify a spatially-varying kernel. 
%   __________________________________________________________________
%
%   Minimum population
%
%   POLYMAP(...,'minpop',NMIN) uses the input bandwidth H, but also
%   guarantees that the population of data points---that is, the number of 
%   data points subject to a nonzero weight---is at least NMIN at each
%   grid point. This is done by expanding the bandwidth H as necessary.  
%   __________________________________________________________________
%
%   Masking out mapping points
%
%   ZHAT=POLYMAP(...,'mask',BOOL), where BOOL is a boolean L x M matrix,
%   suppresses the computation of the polynomial fit at locations where 
%   BOOL is false.  At these locations, ZHAT is set to NaN.
% 
%   This is useful in accelerating POLYMAP by avoiding calculations where
%   no fit is desired.  For example, for mapping ocean quantities BOOL
%   could be set to false over land.
%   __________________________________________________________________
%
%   Fixed population
%
%   ZHAT=POLYMAP(...,{P,HO,ALPHA,BETA},'population',N) specifies a fixed 
%   population, or variable bandwidth, fit with population N and maximum 
%   bandwidth HO.  N can be a scalar or an L x M matrix.
% 
%   Here the bandwidth is varied spatially to be just large enough such 
%   that N points have nonzero weights at every grid point, if possible.
% 
%   The second argument of the cell array, HO, is now the maximum
%   bandwidth over which distances are computed. 
% 
%   The number of data points employed at each mapping point may be smaller
%   than the target population N, if additional points are required to 
%   reach N that lie outside the maximum bandwidth radius HO.
%  
%   To check the population and bandwidth actually used in the fix, see the 
%   first two pages of the auxilliary output variable AUX, described below.
%
%   The fixed population algorithm can give excellent results when the data
%   spacing is uneven, particularly when used with a higher-order fit.
%   __________________________________________________________________
%
%   Smoothing on the surface of the Earth
%
%   ZHAT=POLYMAP(LON,LAT,[],[],Z,LONO,LATO,{P,H,ALPHA,BETA},'sphere')
%   performs the local polynomial fit on the surface of a sphere having the 
%   radius of the Earth.
%
%   The bandwidth H in this case should have units of kilometers. The mean
%   radius of the Earth of 6371 kilometers is used.
%
%   The modification of local polynomial fitting to the sphere is described
%   in XXX.
%   __________________________________________________________________
%
%   Inclusion of temporal variability
%
%   It may be that one wishes to take temporal variability of the 
%   underlying field into account in constructing the estimate ZHAT.
%
%   In this case, data points with values Z are taken at locations X,Y and
%   also at times T. POLYMAP is then called as follows:
%
%       ZHAT=POLYMAP(X,Y,T,[],Z,XO,YO,{P,H,ALPHA,BETA},{MU,TAU,TAL,TBE})
%
%   which will fit the data to the sum of spatial polynomial of bandwidth H
%   and order P, and a temporal polynomial of bandwidth TAU and order MU.
%
%   Here TAL and TBE are the alpha and beta parameters controlling the 
%   kernel K to be used along the time axis, with K~(1-|T/TAU|^TAL)^TBE.
%
%   TAU, TAL, and TBE may all either be scalars or matrices of size L x M.
%   The units of TAU should be the same as those of the times T.  
%
%   The times T should be given relative to the center of the time window,
%   that is, time T=0 should correspond to the time at which you wish to 
%   construct the map. 
%
%   When the 'population' flag is input, described above, it only affects
%   the spatial kernel. The temporal kernel remains specified in terms of 
%   the temporal bandwidth TAU. 
%   __________________________________________________________________
%
%   Weighted data points
%
%   POLYMAP can incorporate an additional weighting factor on the data
%   points. Let W be an array of positive values the same size as the data 
%   array Z.  One may form a map incorporating these weights as follows:
%
%       ZHAT=POLYMAP(X,Y,[],W,Z,XO,YO,{P,H,ALPHA,BETA})
%
%   The weights W could represent the confidence in the measurement values,
%   a weighting used in robustification, or an aggregation of invididual
%   measurements into clusters.  The latter approach may be used to 
%   condense large datasets to a managable size.
%
%   This also works together with the temporal fitting described above.
%
%   When used with the fixed population or minimum population options, each
%   data point contributes an amount to the population specified by W. 
%   __________________________________________________________________
%
%   Additional output arguments
%
%   [ZHAT,BETA,AUX]=POLYMAP(...) returns two additional arguments.
%
%   BETA contains the estimated field, the same as ZHAT, together with 
%   estimates of all spatial derivatives of the field up to Pth order:
% 
%        BETA(:,:,1) = ZHAT     --- The field z
%        BETA(:,:,2) = ZXHAT    --- The first derivative dz/dx
%        BETA(:,:,3) = ZYHAT    --- The first derivative dz/dy
%        BETA(:,:,4) = ZXXHAT   --- The second derivative d^2z/dx^2
%        BETA(:,:,5) = ZXYHAT   --- The second derivative d^2z/dxdy 
%        BETA(:,:,6) = ZYYHAT   --- The second derivative d^2z/dy^2
%
%   The length of the third dimension of BETA is set by the total number of
%   derivatives of order P or less.  This number, denoted Q, is given by 
%   Q = 1, 3, and 6 for P = 0, 1, and 2 respectively. 
%
%   After these, in the case that TS is input, the estimated time 
%   derivatives are returned up to order MU:
%
%       BETA(:,:,Q+1) = ZT   --- The first time derivative dz/dt
%       BETA(:,:,Q+2) = ZTT  --- The second time derivative d^2z/dt^2
%   
%   AUX is an M x N x 5 array of auxiliary fields associated with the fit.
%
%        AUX(:,:,1) = N  --- The number of data points with nonzero weight
%        AUX(:,:,2) = H  --- The bandwidth used at each point
%        AUX(:,:,3) = W  --- The total weight 
%        AUX(:,:,4) = C  --- The matrix condition number
%        AUX(:,:,5) = E  --- The RMS error between the data and the fit
%
%   N, called the population, is the total number of data points less than
%   one bandwidth distance H from each of the (L,M) grid points. 
%
%   The condition number C arises because a matrix must be inverted at each
%   (L,M) location for the P=1 or P=2 fits. C is equal to 1 for P=0. 
%
%   C is computed by COND.  At (L,M) points where C is large, the least 
%   squares solution is unstable, and one should consider using a lower-
%   order fit P or a larger value of the bandwidth H.
%
%   The root-mean-squared error E is computed using the same weighted 
%   kernel that is applied to the data.
%   __________________________________________________________________
%
%   Robustification 
%
%   ZHAT=POLYMAP(...,'robust',NI) implements the median-based iterative 
%   robust algorithm of Cleveland (1979), p 830--831, using NI iterations.  
%
%   This can be useful when outliers are present in the data.  Typically
%   a single iteration is sufficient to remove most outliers. 
%
%   ZHAT will in this case have NI+1 entries along its third dimension.
%   The iterative estimates are stored in reverse order, with the last 
%   iteration in ZHAT(:,:,1) and the original estimate in ZHAT(:,:,NI+1).
%   __________________________________________________________________
%
%   Parallelization
%
%   POLYMAP(...,'parallel') parallelizes the computation using a PARFOR
%   loop, by operating on each latitude (or matrix row) separately.   This 
%   requires that Matlab's Parallel Computing toolbox be installed.  
%
%   POLYMAP will then using an existing parallel pool, or if one does 
%   not exist, a pool will be created using all availabale workers.
%
%   POLYMAP(...'parallel',Nworkers) alternately specifies the number of
%   workers to use. If you run into memory constraints, reduce Nworkers.
%
%   If you are working on multiple maps simultaneously, depending on the 
%   size of your problem, it may be faster to use an exterior PARFOR loop,
%   rather than calling POLYMAP with the 'parallel' flag. 
%   __________________________________________________________________
%
%   Verbose option
%
%   POLYMAP(...,'verbose') displays a status message saying what row it
%   is working on.  POLYMAP(...,'quiet') is the default.
%   _________________________________________________________________
%
%   One grid, many fields
%
%   It is often the case that the field to be mapped, Z, consists of many 
%   repeated sets of observations at the same (X,Y) points. 
%
%   For example, X and Y could be mark the locations of measurements that
%   are repeated at different times (as in satellite altimetry), or else
%   there could be multiple fields Z that are measured simultaneously. 
%   Alternatively, one may simply wish to try out different parameters.
%
%   In such cases, there is a vast computational savings to be gained by
%   only computing quantities which depend on the grid a single time.  
%  
%   The standard spatial mapping of POLYMAP is then performed as follows. 
%
%    (1)  [DS,XS,YS,INDEXS]=PM_SORT(X,Y,XO,YO,HO); 
%    (2)  [AMAT,XMAT,WMAT,H,C]=PM_WEIGHT(DS,XS,YS,[],[],{P,HO,ALPHA,BETA});
%    (3)  Z1S=PM_INDEX(INDEXS,Z1); 
%    (4)  [Z1HAT,BETA,AUX]=PM_APPLY(Z1S,AMAT,XMAT,WMAT,H,C);
%
%   The first call (1) sorts the data points x,y relative to the grid
%   points xo,yo, returning their distances DS up to the bandwidth HO,
%   the x- and y-deviations X-XO and Y-YO in XS and YS respectively, and
%   an index INDEXS.  See PM_SORT for the structure of the output fields.
% 
%   Then (2) creates certain important matrix quantities AMAT, XMAT, and
%   WMAT that depend only on the observations and the grid, and not on the
%   observed values Z.  The output H is an L x M matrix of bandwidths that
%   is modified from HO if any population flags are passed to PM_WEIGHT.
%
%   Given an array of observed values Z1 of the same size as X and Y, (3)
%   sorts those values into Z1S, which is the same size as DS, XS, etc.
%
%   Finally, (4) uses the quantities created in step (2) together with the
%   sorted observational values from step (3) to create the final map.
%
%   Steps (3) and (4) can then be iterated for Z2, Z3, Z4, etc, provided
%   these are all on the same grid as Z1. In typical problems steps (1)
%   and (2) take a least an order of magnitude more time that (3) and (4).
% 
%   Moreover, because the least squares problem is essentially already 
%   solved in (2), steps (3) and (4) do not increase in computational cost
%   as the fit order increases.
%
%   A subtlety is the case in which the grid is the same, but Z1, Z2, etc.
%   may have missing data at different points. This will leads to NaNs in
%   DS etc. that propagate through to the final map.  This cannot be
%   solved by e.g. swapping NaNs for 0s in Z1S, Z2S, etc. because the grid
%   itself has in fact changed, and the output of (2) should reflect this.
% 
%   If not too much of the grid is affected, one workaround is to first 
%   compute (1)--(4) for an entire grid with no missing data, and then 
%   re-do all of steps (1)--(4) as needed for each Z1, Z2, etc., only for
%   those grid points for which the original fit yeilded NaNs.
%
%   Some of the other options for POLYMAP are implemented as follows.
%
%   For mapping on the sphere, only (1) needs to be changed:
%
%    (1) [DS,XS,YS,INDEXS]=PM_SORT(LON,LAT,LONO,LATO,HO,'sphere');  
%
%   For weighted data points, fixed population, or minimum population, 
%   only (2) needs to be changed:
%
%    (2)  [AMAT,XMAT,WMAT,H,C]=PM_WEIGHT(DS,XS,YS,[],WS,{P,HO,ALPHA,BETA});
%    (2)  [AMAT,XMAT,WMAT,H,C]=...
%            PM_WEIGHT(DS,XS,YS,[],[],{P,HO,ALPHA,BETA},'population',N);
%    (2)  [AMAT,XMAT,WMAT,H,C]=...
%            PM_WEIGHT(DS,XS,YS,[],[],{P,HO,ALPHA,BETA},'minpop',NMIN);
%
%   For temporal fitting, one should add a step (0) before (1) 
% 
%    (0)  [X,Y,T,W,Z1,Z2,...,ZK]=PM_WINDOW(X,Y,T,W,Z1,Z2,...,ZK,TAU);
% 
%   which remove out-of-range time points and then change (2) to be
%
%    (2)  [AMAT,XMAT,WMAT,H,C]=...
%            PM_WEIGHT(DS,XS,YS,TS,[],{P,HO,ALPHA,BETA},{MU,TAU,TAL,TBE});
%
%   For parallelization, each of (1)--(4) can be called with the 'parallel'
%   flag. Each can also be called with the 'verbose' flag. 
%   __________________________________________________________________
%
%   'polymap --t' runs some tests.
%   'polymap --f' generates a sample figure.
%
%   Usage: zhat=polymap(x,y,[],[],z,xo,yo,{P,H,alpha,beta});
%          zhat=polymap(x,y,t,[],z,xo,yo,{P,H,alpha,beta},{mu,tau,tal,tbe});
%          zhat=polymap(lon,lat,[],[],z,lono,lato,{P,H,alpha,beta},'sphere');
%          [zhat,beta,aux]=polymap(x,y,[],[],z,xo,yo,{P,H,alpha,beta});
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2022 J.M. Lilly --- type 'help jlab_license' for details
 
%   'polymap --f2' with jData installed generates the figure shown
%       above, which may require a relatively powerful computer.

%        AUX(:,:,5) = R  --- The weighted mean distance to the data points
%        AUX(:,:,6) = V  --- The weighted standard deviation of data values

%   The derivatives appearing in BETA are now the derivatives given in a
%   local tangent plane.  These can be converted into derivatives in terms
%   of latitude and longitude following Lilly and Lagerloef (2018). XXX?

if strcmpi(varargin{1}, '--t')
    polymap_test,return
elseif strcmpi(varargin{1}, '--f')
    type makefigs_polymap
    makefigs_polymap;
    return
elseif strcmpi(varargin{1}, '--f2')
    type makefigs_polymap2
    makefigs_polymap2;
    return
end

parstr='series';
robstr='non';
verbstr='quiet';
varstr='bandwidth';
geostr='cartesian';
iters=0;
Nworkers=[];
Pmin=[];
bool=[];
Npop=nan;

%First parse the string arguments
for i=1:7
    if ischar(varargin{end})
        tempstr=varargin{end};
        if strcmpi(tempstr(1:3),'ver')||strcmpi(tempstr(1:3),'qui')
            verbstr=tempstr;
        elseif strcmpi(tempstr(1:3),'ser')||strcmpi(tempstr(1:3),'par')
            parstr=tempstr;
        elseif strcmpi(tempstr(1:3),'car')||strcmpi(tempstr(1:3),'sph')
            geostr=tempstr;
        end
        varargin=varargin(1:end-1);
    elseif ~ischar(varargin{end})&&ischar(varargin{end-1})
        tempstr=varargin{end-1};
        if strcmpi(tempstr(1:3),'mas')
            bool=varargin{end};
        elseif strcmpi(tempstr(1:3),'rob')||strcmpi(tempstr(1:3),'non')
            robstr=tempstr;
            iters=varargin{end};
        elseif strcmpi(tempstr(1:3),'pop')
            varstr=tempstr;
            Npop=varargin{end};
        elseif strcmpi(tempstr(1:3),'par')
            parstr=tempstr;
            Nworkers=varargin{end};
        end
        varargin=varargin(1:end-2);
    end
end

pool = gcp('nocreate');
if strcmpi(parstr(1:3),'par')
    if isempty(Nworkers)
        if isempty(pool)
            parpool('local');
        end
    else
        if ~isempty(pool)
            parpool('local',Nworkers);
        elseif pool.NumWorkders~=Nworkers
            parpool('local',Nworkers);
        end
    end
end

%--------------------------------------------------------------------------
%sort out rest of input arguments
xdata=varargin{1};
ydata=varargin{2};
xo=varargin{6}(:)';
yo=varargin{7}(:);

tdata=varargin{3};
wdata=varargin{4};
zdata=varargin{5};
sarg=varargin{8};

if length(nargin)==9
    targ=varargin{9};
else
    targ=[];
end

if isempty(bool)
    bool=true(length(yo),length(xo));
end

%--------------------------------------------------------------------------
%remove non-finite data points
bool_finite=isfinite(xdata)&isfinite(ydata)&isfinite(zdata);
xdata=xdata(bool_finite);
ydata=ydata(bool_finite);
zdata=zdata(bool_finite);
if ~isempty(tdata)
    tdata=tdata(bool_finite);
end   
if ~isempty(wdata)
    wdata=wdata(bool_finite);
end  

%--------------------------------------------------------------------------
%initialize output arrays to correct sizes
L=length(yo);
M=length(xo);

if ~isempty(targ)
    mu=targ{1};
else
    mu=0;
end

Q=sum(0:sarg{1}+1)+mu;

if isempty(bool)
    bool=true(L,M);
end

if numel(Npop)==1
    Npop=Npop+zeros(L,M);
end

%--------------------------------------------------------------------------
%Rearrange information in kernel cells
for i=2:4
    if length(sarg{i})==1
        sarg{i}=sarg{i}+zeros(L,M);
    end
end
scell=cell(L,1);
for i=1:L
    scell{i}={sarg{1},sarg{2}(i,:),sarg{3}(i,:),sarg{4}(i,:)};
end
tcell=cell(L,1);

if ~isempty(targ)
    for i=2:4
        if length(targ{i})==1
            targ{i}=targ{i}+zeros(L,M);
        end
    end
    for i=1:L
        tcell{i}={targ{1},targ{2}(i,:),targ{3}(i,:),targ{4}(i,:)};
    end
end

[zhat,beta,aux,res]=polymap_loop(L,M,Q,xdata,ydata,tdata,wdata,zdata,...
    xo,yo,scell,tcell,bool,geostr,varstr,Npop,verbstr,parstr,iters);


function[zhat,beta,aux,res]=polymap_loop(L,M,Q,xdata,ydata,tdata,wdata,zdata,xo,yo,scell,tcell,bool,geostr,varstr,Npop,verbstr,parstr,iters)

zhat=zeros(L,M,iters+1);
beta=zeros(L,M,Q);
aux=zeros(L,M,5);
%res=zdata;

%--------------------------------------------------------------------------
%execute loop
if ~isempty(xdata)
    if  strcmpi(parstr(1:3),'ser')
        for i=1:L
            if strcmpi(verbstr(1:3),'ver')
                disp(['POLYMAP computing fit for row ' int2str(i) ' of ' int2str(L) '.'])
            end
            [zhat(i,:,:),beta(i,:,:),aux(i,:,:),res{i}]=...
                polymap_one(xdata,ydata,tdata,wdata,zdata,xo,yo(i),scell{i},tcell{i},geostr,varstr,Npop(i,:),bool(i,:),iters);
        end
    elseif strcmpi(parstr(1:3),'par')
        parfor i=1:L
            if strcmpi(verbstr(1:3),'ver')
                disp(['POLYMAP computing fit for row ' int2str(i) ' of ' int2str(L) '.'])
            end
            [zhat(i,:,:),beta(i,:,:),aux(i,:,:),res{i}]=...
                polymap_one(xdata,ydata,tdata,wdata,zdata,xo,yo(i),scell{i},tcell{i},geostr,varstr,Npop(i,:),bool(i,:),iters);
        end
    end
    for i=1:L
        res{i}=res{i}{1};
    end
end

function[zhat,beta,aux,res]=polymap_one(xdata,ydata,tdata,wdata,zdata,xo,yo,scell,tcell,geostr,varstr,Npop,bool,iters)

if ~isempty(tcell)
    [xdata,ydata,tdata,wdata,zdata]=pm_window(xdata,ydata,tdata,wdata,zdata,tcell{2});
end

[ds,xs,ys,zs,ts,ws]=pm_sort(xdata,ydata,zdata,tdata,wdata,xo,yo,scell{2},geostr,'mask',bool,varstr,Npop); 
%Note that PM_WEIGHT can be called with ZS in the place of INDEXS
[zs,amat,xmat,wmat,H,C]=pm_weight(ds,xs,ys,zs,ts,ws,scell,tcell,varstr,Npop);
zhat=zeros(length(ds),size(ds{1},2),iters+1);
[zhat(:,:,iters+1),beta,aux,res]=pm_apply(zs,amat,xmat,wmat,H,C);

%boolinf=isinf(zhat(:,:,iters+1));

%robustification, see Cleveland (1979), p 830--831
while iters>0
    %update weight with result of PM_ROBUST
    [~,amat,xmat,wmat,H,C]=pm_weight(ds,xs,ys,zs,ts,pm_robust(ws,res),scell,tcell,varstr,Npop);
    %length(find(isinf(zs{1}(:))))
    [zhat(:,:,iters),beta,aux,res]=pm_apply(zs,amat,xmat,wmat,H,C);
    iters=iters-1;
end

function[]=polymap_test

polymap_test_tangentplane;
polymap_test_peaks_cartesian;
polymap_test_peaks_spherical;
polymap_test_peaks_cartesian_robust;

function[]=polymap_test_tangentplane

%Testing tangent plane equations 
load goldsnapshot
use goldsnapshot

phip=vshift(lat,-1,1);
phin=vshift(lat,1,1);
thetap=vshift(lon,-1,1);
thetan=vshift(lon,1,1);

[long,latg]=meshgrid(lon,lat);
[thetapg,phipg]=meshgrid(thetap,phip);
[thetang,phing]=meshgrid(thetan,phin);

xp=radearth*cosd(latg).*sind(thetapg-long);
xn=radearth*cosd(latg).*sind(thetang-long);
yp=radearth*(cosd(latg).*sind(phipg)-sind(latg).*cosd(phipg));%.*cosd(thetapg-long));
yn=radearth*(cosd(latg).*sind(phing)-sind(latg).*cosd(phing));%.*cosd(thetang-long));

for i=1:4
    if i==1
        ssh=latg;
    elseif i==2
        ssh=squared(latg);
    elseif i==3
        ssh=squared(long);
    elseif i==4
        ssh=goldsnapshot.ssh;
    end
    
    sshp=vshift(ssh,-1,2);
    sshn=vshift(ssh,1,2);
    dZdx=frac(sshn-sshp,xn-xp);
    
    sshp=vshift(ssh,-1,1);
    sshn=vshift(ssh,1,1);
    dZdy=frac(sshn-sshp,yn-yp);
    
    [fx,fy]=spheregrad(lat,lon,ssh);
    fx=fx*1000;fy=fy*1000;
 
    %figure,
    %subplot(1,3,1),jpcolor(fx),caxis([-1 1]/5)
    %subplot(1,3,2),jpcolor(dZdx),caxis([-1 1]/5)
    %subplot(1,3,3),jpcolor(fx-dZdx),caxis([-1 1]/1000/1000/5)
    
    %figure
    %subplot(1,3,1),jpcolor(fy),caxis([-1 1]/5)
    %subplot(1,3,2),jpcolor(dZdy),caxis([-1 1]/5)
    %subplot(1,3,3),jpcolor(fy-dZdy),caxis([-1 1]/1000/1000/5)
    
    %Gradients agree to like one part in one million
    
    del2=spherelap(lat,lon,ssh)*1000*1000;
    
    sshp=vshift(ssh,-1,2);
    sshn=vshift(ssh,1,2);
    %d2Zdx2=frac(2,xn-xp).*(frac(sshn-ssh,xn)-frac(ssh-sshp,-xp));
    d2Zdx2=frac(1,squared((xn-xp)/2)).*(sshn+sshp-2*ssh);
    
    sshp=vshift(ssh,-1,1);
    sshn=vshift(ssh,1,1);
    %d2Zdy2=frac(2,yn-yp).*(frac(sshn-ssh,yn)-frac(ssh-sshp,-yp));
    d2Zdy2=frac(1,squared((yn-yp)/2)).*(sshn+sshp-2*ssh);
    
    del2hat=d2Zdx2+d2Zdy2-frac(tand(lat),radearth).*dZdy;
    %del2hat=d2Zdx2+d2Zdy2;%-frac(tand(lat),radearth).*dZdy;
        
    if i==1
        index=find(abs(latg)<89.5);
        bool=allall(log10(abs(del2hat(index)-del2(index))./abs(del2(index)))<-5);
        reporttest('POLYMAP Laplacian of linear function of y',bool)
    elseif i==2
        %figure,plot(lat,log10(abs(del2))),hold on,plot(lat,log10(abs(del2-del2hat)))
        index=find(abs(latg)<88.5);
        bool=allall(log10(abs(del2hat(index)-del2(index)))<-6);
        reporttest('POLYMAP Laplacian of quadratic function of y',bool)
    elseif i==3
        %figure, plot(lon,log10(abs(del2-del2hat)./abs(del2))')
        index=find(abs(long)<178.5&isfinite(del2hat));
        %length(find(~isfinite(del2)))
        %length(find(~isfinite(del2hat)))
        bool=allall(log10(abs(del2hat(index)-del2(index))./abs(del2(index)))<-5);
        reporttest('POLYMAP Laplacian of quadratic function of x',bool)
    elseif i==4
        xx=log10(abs(fx-dZdx)./abs(fx));
        index=find(isfinite(xx));
        bool=allall(xx(index)<-5);
        reporttest('POLYMAP longitude gradient of GOLDSNAPSHOT',bool)
        yy=log10(abs(fy-dZdy)./abs(fy));
        index=find(isfinite(yy));
        bool=allall(yy(index)<-5);
        reporttest('POLYMAP latitude gradient of GOLDSNAPSHOT',bool)
        %figure
        %subplot(1,3,1),jpcolor(log10(abs(del2))),caxis([-3 0])
        %subplot(1,3,2),jpcolor(log10(abs(del2hat))),caxis([-3 0])
        %subplot(1,3,3),jpcolor(log10(abs(del2-del2hat))),caxis([-3 0])
    end
end


function[]=polymap_test_peaks_cartesian

%Use peaks for testing... random sampling
[x,y,z]=peaks;
rng(1);
index=randperm(length(z(:)));
index=index(1:600);
[xdata,ydata,zdata]=vindex(x(:),y(:),z(:),index,1);

xo=(-3:.125:3);
yo=(-3:.125:3);

H=0.75; %Slightly different from example in figure
%H=0.3

zhat=zeros(length(yo),length(xo),3);
zhat(:,:,1)=polymap(xdata,ydata,[],[],zdata,xo,yo,{0,H,2,1});
[zhat(:,:,2),beta1]=polymap(xdata,ydata,[],[],zdata,xo,yo,{1,H,2,1});
[zhat(:,:,3),beta2]=polymap(xdata,ydata,[],[],zdata,xo,yo,{2,H,2,1});

[ds,xs,ys,zs]=pm_sort(xdata,ydata,zdata,xo,yo,H);

zhat2=zeros(length(yo),length(xo),3);
for i=1:3
    [zs,amat,xmat,wmat,H,C]=pm_weight(ds,xs,ys,zs,[],[],{i-1,H,2,1});
    zhat2(:,:,i)=pm_apply(zs,amat,xmat);
    %[zhat2(:,:,i),~,~,res]=pm_apply(zs,amat,xmat,wmat,H,C);
end

reporttest('POLYMAP internal and external loop consistency',aresame(zhat,zhat2,1e-10))

zvar=mean(mean(squared(z)));
err0=mean(mean(squared(zhat(:,:,1)-z)))./zvar;
err1=mean(mean(squared(zhat(:,:,2)-z)))./zvar;
err2=mean(mean(squared(zhat(:,:,3)-z)))./zvar;

bool=(err0<0.06)&&(err1<0.055)&&(err2<0.01);

reporttest('POLYMAP Cartesian peaks mapping accuracy',bool)

%next check gradients
dzdx=vdiff(z,2)./(xo(2)-xo(1));
dzdy=vdiff(z,1)./(yo(2)-yo(1));

errx1=mean(mean(squared(beta1(:,:,2)-dzdx)))./mean(mean(squared(dzdx)));
erry1=mean(mean(squared(beta1(:,:,3)-dzdy)))./mean(mean(squared(dzdy)));
errx2=mean(mean(squared(beta2(:,:,2)-dzdx)))./mean(mean(squared(dzdx)));
erry2=mean(mean(squared(beta2(:,:,3)-dzdy)))./mean(mean(squared(dzdy)));

bool=(errx1<0.09)&&(erry1<0.09)&&(errx2<0.07)&&(erry2<0.07);

reporttest('POLYMAP Cartesian peaks first derivative accuracy',bool)


function[]=polymap_test_peaks_spherical

%Use peaks for testing... random sampling
[x,y,z]=peaks;
rng(1);
index=randperm(length(z(:)));
index=index(1:600);
[xdata,ydata,zdata]=vindex(x(:),y(:),z(:),index,1);

londata=xdata*180/3;
latdata=ydata*180/6*85/90;
lono=(-3:.125:3)*180/3;
lato=(-3:.125:3)*180/6*85/90;

H=3000;

zhat=zeros(length(lato),length(lono),3);
zhat(:,:,1)=polymap(londata,latdata,[],[],zdata,lono,lato,{0,H,2,1},'sphere');
[zhat(:,:,2),beta1]=polymap(londata,latdata,[],[],zdata,lono,lato,{1,H,2,1},'sphere');
[zhat(:,:,3),beta2]=polymap(londata,latdata,[],[],zdata,lono,lato,{2,H,2,1},'sphere');

zvar=mean(mean(squared(z)));
err0=mean(mean(squared(zhat(:,:,1)-z)))./zvar;
err1=mean(mean(squared(zhat(:,:,2)-z)))./zvar;
err2=mean(mean(squared(zhat(:,:,3)-z)))./zvar;

bool=(err0<0.09)&&(err1<0.75)&&(err2<0.01);

reporttest('POLYMAP spherical peaks mapping accuracy',bool)

%next check gradients
[dzdx,dzdy]=spheregrad(lato,lono,z);

errx1=mean(mean(squared(beta1(:,:,2)/1000-dzdx)))./mean(mean(squared(dzdx)));
erry1=mean(mean(squared(beta1(:,:,3)/1000-dzdy)))./mean(mean(squared(dzdy)));
errx2=mean(mean(squared(beta2(:,:,2)/1000-dzdx)))./mean(mean(squared(dzdx)));
erry2=mean(mean(squared(beta2(:,:,3)/1000-dzdy)))./mean(mean(squared(dzdy)));

bool=(errx1<0.07)&&(erry1<0.15)&&(errx2<0.06)&&(erry2<0.15);

reporttest('POLYMAP spherical peaks first derivative accuracy',bool)


function[]=polymap_test_peaks_cartesian_robust

%Use peaks for testing... random sampling
[x,y,z]=peaks;
rng(1);
index=randperm(length(z(:)));
index=index(1:600);
[xdata,ydata,zdata]=vindex(x(:),y(:),z(:),index,1);

xo=(-3:.125:3);
yo=(-3:.125:3);

H=0.75; %Slightly different from example in figure

zhat0=polymap(xdata,ydata,[],[],zdata,xo,yo,{0,H,2,1});
zhat1=polymap(xdata,ydata,[],[],zdata,xo,yo,{1,H,2,1});

zdata(end/2)=100;

zhat0r=polymap(xdata,ydata,[],[],zdata,xo,yo,{0,H,2,1},'robust',3);
err0=squeeze(mean(mean(squared(zhat0r-zhat0)))./mean(mean(squared(zhat0))));

zhat1r=polymap(xdata,ydata,[],[],zdata,xo,yo,{1,H,2,1},'robust',3);
err1=squeeze(mean(mean(squared(zhat1r-zhat1)))./mean(mean(squared(zhat1))));

bool0=all(err0(1:3)<0.003);
bool1=all(err1(1:3)<0.002);

reporttest('POLYMAP robustification',bool0&&bool1)


