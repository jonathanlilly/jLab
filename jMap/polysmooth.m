function[varargout]=polysmooth(varargin)
%POLYSMOOTH  Mapping using local polynomial fitting, also known as loess.  
%
%   POLYSMOOTH generates a map from scattered data in two dimensions using
%   a locally weighted least squares fit to a polynomial.
%
%   This method is variously known as local polynomial fitting, local
%   polynomial smoothing, multivariate locally weighted least squares 
%   regression, lowess (originally for LOcally WEighted Scatterplot
%   Smoothing), and loess.  All of these are essentially synonyms.
%
%   POLYSMOOTH has the support for all of the following options:
%
%       --- Cartesian or spherical geometry
%       --- a constant, linear, or quadratic fit in space
%       --- an additional linear or quadartic fit in time
%       --- fixed-bandwidth or fixed population algorithms
%       --- prescribed spatially-varying bandwidth or population
%       --- multiple choices of weighting function or kernel
%       --- additional datapoint weighting factors, e.g. by confidence
%       --- the median-based robustification of Cleveland (1979)
%
%   POLYSMOOTH is implemented using a numerically efficient algorithm that
%   avoids using explicit loops. The data are pre-sorted so that different
%   mapping parameters can be tried out at little computational expense.
%
%   For algorithm details, see Lilly and Lagerloef (2018).
%   __________________________________________________________________
%
%   Local polynomial fitting on the plane
%
%   Let's say we have an array Z of data is at locations X,Y.  X,Y, and Z
%   can be arrays of any size provided they are all the same size.
%
%   The problem is to obtain a mapped field ZHAT on some regular grid 
%   specified by the vectors XO and YO.
%
%   Calling POLYSMOOTH is a two-step process: 
%
%       [DS,XS,YS,ZS]=TWODSORT(X,Y,Z,XO,YO,CUTOFF);
%       ZHAT=POLYSMOOTH(DS,XS,YS,[],ZS,[],RHO,B);
%
%   The empty arrays mark locations of optional arguments described later.
%
%   In the first step, one calls TWODSORT which returns ZS, a 3D array of 
%   data values at each grid point, sorted by increasing distance DS, and 
%   the corresponding positions XS and YS.  See TWODSORT for more details.
%  
%   CUTOFF determines the maximum distance included in the sorting and 
%   should be chosen to be greater than B.  
%
%   In the second step, POLYSMOOTH fits a RHOth order spatial polynomial at 
%   each gridpoint within a neighborhood specified by the "bandwidth" B.
%
%   The fit is found by minimizing the weighted mean squared error between 
%   the fitted surface and the observations.  The bandwidth B sets the 
%   decay of this weighting function, described in more detail shortly.
%
%   The fit order RHO may be chosen as RHO=0 (fit to a constant), RHO=1
%   (fit to a plane), or else RHO=2 (fit to a parabolic surface). 
%
%   The data locations (X,Y) and grid point locations (X0,Y0) shoud have
%   the same units as the bandwidth B (e.g., kilometers).
%
%   Note that B may either be a scalar, or a matrix of size M x N to
%   specify an imposed spatially-varying bandwidth. 
%
%   The dimensions of XO and YO are M x N x J, where J is the maximum
%   number of data points within bandwidth cutoff at any grid point 
%   location.  Then ZHAT is matrix of dimension M x N.
%   __________________________________________________________________
%
%   Choice of weighting function or kernel
%
%   POLYSMOOTH(DS,XS,YS,[],ZS,[],RHO,{B,KERN}) weights the data points in
%   the vicinity of each grid point by some decaying function of distance 
%   called the kernel, specified by KERN. 
%
%        KERN = 'uni' uses the uniform kernel, K=1
%        KERN = 'epa' uses the Epanechnikov (parabolic) kernel K=1-(DS/B)^2
%        KERN = 'bis' uses the bisquare kernel K=(1-(DS/B)^2)^2
%        KERN = 'tri' uses the tricubic kernel K=(1-(DS/B)^3)^3
%        KERN = 'gau' used the Gaussian kernel, K=EXP(-1/2*(3*DS/B)^2)
%
%   Note that all choices of weighting function are set to vanish for DS>B.
%
%   The default behavior is STR='gau' for the Gaussian kernel, which is 
%   specified to have a standard deviation of B/3.
%
%   KERN may also be an integer, in which the standard deviation of the 
%   Gaussian is set to B/KERN, with the default corresponding to KERN = 3.
%
%   KERN may also be a vector (of length greater than one) defining a 
%   custom kernel, with KERN(1) corresponding to DS/B=0 and KERN(end) to 
%   DS/B=1.  Kernel values at any distance are then linearly interpolated.
%   __________________________________________________________________
%
%   Inclusion of temporal variability
%
%   It may be that the data contributing to the map is taken at different
%   times, and that in constructing the map it is important to take into
%   account temporal variability of the underlying field.
%
%   In this case, data points with values Z are taken at locations X,Y and
%   also at times T. POLYSMOOTH is then called as follows:
%
%       [DS,XS,YS,TS,ZS]=TWODSORT(X,Y,T,Z,XO,YO,CUTOFF);
%       ZHAT=POLYSMOOTH(DS,XS,YS,TS,ZS,[],[RHO MU],B,TAU);
%
%   which will fit the data to the sum of spatial polynomial of bandwidth B
%   and order RHO, and a temporal polynomial of bandwidth TAU and order MU.
%
%   TAU, like B, may either be a scalar or an M x N matrix.  Its units
%   should be the same as those of the times T.  
%
%   Note that the times T should be given relative to the center of the 
%   time window, that is, time T=0 should correspond to the time at which
%   you wish to construct the map. 
%   
%   By default the Gaussian kernal is used in time.  One can also employ
%   a cell array, as in POLYSMOOTH(..,[RHO MU],B,{TAU,KERN}), to specify 
%   other behaviors for the time kernel, as described above.
%   __________________________________________________________________
%
%   Additional output arguments
%
%   [ZHAT,BETA,AUX]=POLYSMOOTH(...) returns two additional arguments.
%
%   BETA contains the estimated field, the same as ZHAT, together with 
%   estimates of all spatial derivatives of the field up to RHOth order:
% 
%        BETA(:,:,1) = ZHAT     --- The field z
%        BETA(:,:,2) = ZXHAT    --- The first derivative dz/dx
%        BETA(:,:,3) = ZYHAT    --- The first derivative dz/dy
%        BETA(:,:,4) = ZXXHAT   --- The second derivative d^2z/dx^2
%        BETA(:,:,5) = ZXYHAT   --- The second derivative d^2z/dxdy 
%        BETA(:,:,6) = ZYYHAT   --- The second derivative d^2z/dy^2
%
%   The length of the third dimension of BETA is set by the total number of
%   derivatives of order RHO or less.  This number, called Q, is equal to
%   Q = 1, 3, and 6 for RHO = 0, 1, and 2 respectively. 
%
%   After these, in the case that TS is input, the estimated time 
%   derivatives are returned up to order MU:
%
%       BETA(:,:,Q+1) = ZT   --- The first time derivative dz/dt
%       BETA(:,:,Q+2) = ZTT  --- The second time derivative d^2z/dt^2
%   
%   AUX is an M x N x 5 array of auxiliary fields associated with the fit.
%
%        AUX(:,:,1) = P  --- The total number of datapoints employed
%        AUX(:,:,2) = B  --- The bandwidth used at each point
%        AUX(:,:,3) = E  --- The rms error between the data and the fit
%        AUX(:,:,4) = W  --- The total weight used 
%        AUX(:,:,5) = R  --- The mean distance to the data points
%        AUX(:,:,6) = V  --- The weighted standard deviation of data values
%        AUX(:,:,7) = C  --- The matrix condition number
%
%   P, called the population, is the total number of data points within one
%   bandwidth distance B from each of the (M,N) grid points. 
%
%   The root-mean-squared error E and standard deviation V are both 
%   computed using the same weighted kernal applied to the data.
%
%   The condition number C arises because a matrix must be inverted at each
%   (M,N) location for the RHO=1 or RHO=2 fits. C is equal to 1 for RHO=0. 
%
%   C is computed by COND.  At (M,N) points where C is large, the least 
%   squares solution is unstable, and one should consider using a lower-
%   order fit RHO or a larger value of the bandwidth B.
%   __________________________________________________________________
%
%   Fixed population
%
%   POLYSMOOTH(DS,XS,YS,[],ZS,[],RHO,P,'population') varies the spatial 
%   bandwidth B to be just large enough at each grid point to encompass P 
%   points. This is referred to here as the "fixed population" algorithm.
%
%   Note that the argument P relaces the bandwidth B.  So, if one chooses 
%   to specify a kernel, one use POLYSMOOTH(...,RHO,{P,KERN},'population'). 
%
%   When employed with the option for including a temporal fit, the fixed
%   population algorithm only applies to the spatial kernel.  The temporal 
%   kernel remains specified in terms of a temporal bandwidth. 
%
%   The fixed population algorithm can give good results when the data
%   spacing is uneven, particularly when used with a higher-order fit.
%
%   When using this method, the length of the third dimension of the fields
%   output by TWODSORT or SPHERESORT must be at least P. If it is greater 
%   than P, it may be truncated to exactly P, thus reducing the size of 
%   those fields and speeding up calculations.  This is done internally,
%   and also can be done externally by calling POLYSMOOTH_PRESORT.
%   __________________________________________________________________
%
%   Weighted data points
%
%   POLYSMOOTH can incorporate an additional weighting factor on the data
%   points. Let W be an array of positive values the same size as the data 
%   array Z.  One may form a map incorporating these weights as follows:
%
%       [DS,XS,YS,ZS,WS]=TWODSORT(X,Y,Z,W,XO,YO,CUTOFF);           
%       ZHAT=POLYSMOOTH(DS,XS,YS,[],ZS,WS,RHO,B);
%
%   The weights W could represent the confidence in the measurement values,
%   or an aggregation of invididual measurements into clusters.  The latter 
%   approach may be used to condense large datasets to a managable size.
%   __________________________________________________________________
%
%   Smoothing on the sphere
%
%   POLYSMOOTH supports a local polynomial fit on the sphere, as described
%   in Lilly and Lagerloef (2018).  As before this is a two-step process:
%
%       [DS,XS,YS,ZS]=SPHERESORT(LAT,LON,Z,LATO,LONO,CUTOFF);
%       ZHAT=POLYSMOOTH(DS,XS,YS,[],ZS,[],RHO,B);
%
%   The only different is that one firstly one calls SPHERESORT, the
%   analogue of TWODSORT for the sphere.  See SPHERESORT for more details.
%
%   The bandwidth in this case should have units of kilometers. 
%
%   Note that SPHERESORT and POLYSMOOTH both assume the sphere to be the 
%   radius of the earth, as specified by the function RADEARTH.
%
%   The derivatives appearing in BETA are now the derivatives given in a
%   local tangent plane.  These can be converted into derivatives in terms
%   of latitude and longitude following Lilly and Lagerloef (2018).
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
%
%   There is a simple way to handle this situation without needing to 
%   resort the grid.  First one calls TWODSORT or SPHERESORT as follows:
%
%     [DS,XS,YS,INDEX]=TWODSORT(X,Y,[1:LENGTH(X(:))],XO,YO,CUTOFF);
%     --- or ---
%     [DS,XS,YS,INDEX]=
%            SPHERESORT(LAT,LON,1:LENGTH(LAT(:)),LATO,LONO,CUTOFF);
%
%   INDEX is now an index into the sorted datapoint locations, such that
%
%      ZS=POLYSMOOTH_INDEX(SIZE(DS),INDEX,K);
%
%   returns sorted values of Z that can be passed to POLYSMOOTH.  
%  
%   The virtue of this approach is that one only has to call TWODSORT or 
%   SPHERESORT once, no matter how many variable are to be mapped.
%
%   One tip for using this method: (X,Y) or (LAT,LON) values for which the 
%   corresponding data is always undefined (e.g., altimeter tracks over 
%   land), may be set to NaNs in the calls to TWODSORT or SPHERESORT, such
%   that they will be omitted from DS, XY, YS, and INDEX.
%   __________________________________________________________________
%
%   Robustification 
%
%   ZHAT=POLYSMOOTH(...,'robust',NI) implements the median-based iterative 
%   robust algorithm of Cleveland (1979), p 830--831, using NI iterations.  
%
%   This can be useful when outliers are present in the data.   Typically
%   a single iteration is sufficient to remove most outliers. 
%
%   ZHAT will in this case have NI+1 entries along its third dimension.
%   The iterative estimates are stored in reserve order, with the last 
%   iteration in ZHAT(:,:,1) and the original estimate in ZHAT(:,:,NI+1).
%   __________________________________________________________________
%
%   Parallelization
%
%   POLYSMOOTH(...,'parallel') parallelizes the computation using a PARFOR
%   PARFOR loop, by calling the default speed-optimized algorithm on each 
%   latitude (or matrix row) separately.  
%
%   Whether this is faster or slower than the default algorithm depends on
%   the memory needed by the problem relative to that of your machine.  If 
%   the memory is too large, the parallel method will become very slow.
%
%   Alternatively, it may sometimes be preferable to use the default
%   algorithm but to parallelize an external call to POLYSMOOTH.
%
%   This requires that Matlab's Parallel Computing toolbox be installed.
%   __________________________________________________________________
%
%   'polysmooth --t' runs some tests.
%   'polysmooth --f' generates some sample figures.
%
%   Usage:  [ds,xs,ys,zs]=twodsort(x,y,z,xo,yo,cutoff);  
%           zhat=polysmooth(ds,xs,ys,[],zs,[],rho,B);
%           [zhat,beta,aux]=polysmooth(ds,xs,ys,[],zs,[],rho,B);
%   --or--
%           [ds,xs,ys,zs]=spheresort(lat,lon,z,w,lato,lono,cutoff); 
%           [zhat,beta,aux]=polysmooth(ds,xs,ys,[],zs,[],rho,B);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2018 J.M. Lilly --- type 'help jlab_license' for details
 
%   'polysmooth --f2' with jData installed generates the figure shown
%       above, which may require a relatively powerful computer.

if strcmpi(varargin{1}, '--t')
    polysmooth_test,return
elseif strcmpi(varargin{1}, '--f')
    type makefigs_polysmooth
    makefigs_polysmooth;
    return
elseif strcmpi(varargin{1}, '--f2')
    type makefigs_polysmooth2
    makefigs_polysmooth2;
    return
end

mu=0;
targ=[];
tau=[];
str='speed';
skern='gaussian';
tkern='gaussian';
geostr='cartesian';
varstr='bandwidth';
robstr='non';
iters=0;

%First parse the string arguments
for i=1:4
    if ischar(varargin{end})
        tempstr=varargin{end};
        if strcmpi(tempstr(1:3),'car')||strcmpi(tempstr(1:3),'sph')
            geostr=tempstr;
        elseif strcmpi(tempstr(1:3),'spe')||strcmpi(tempstr(1:3),'loo')||strcmpi(tempstr(1:3),'par')
            str=tempstr;
        elseif strcmpi(tempstr(1:3),'ban')||strcmpi(tempstr(1:3),'pop')
            varstr=tempstr;
        end
        varargin=varargin(1:end-1);
    elseif ~ischar(varargin{end})&&ischar(varargin{end-1})
        tempstr=varargin{end-1};
        if strcmpi(tempstr(1:3),'rob')||strcmpi(tempstr(1:3),'non')
            robstr=tempstr;
            iters=varargin{end};
        end
        varargin=varargin(1:end-2);
    end
end

%These arguments are different between the loop and fast algorithm
if strcmpi(str(1:3),'loo')
    xo=varargin{6};
    yo=varargin{7};
    varargin=varargin([1:5 8:end]);
else
    d=varargin{1};
    varargin=varargin(2:end);
end

x=varargin{1};
y=varargin{2};
t=varargin{3};
z=varargin{4};
w=varargin{5};
rho=varargin{6};
sarg=varargin{7};
if length(varargin)==8
    targ=varargin{8};
end
 
if iscell(sarg)
    B=sarg{1};
    skern=sarg{2};
else
    B=sarg;
end

if iscell(targ)
    tau=targ{1};
    tkern=targ{2};
else
    tau=targ;
end

%this might be of length 2
if length(rho)==2
    mu=rho(2);
    rho=rho(1);
end

%Issue a warning if both T and TAU are empty or non-empty
if ~allall([isempty(t),isempty(tau)])&&...
    ~all([~isempty(t),~isempty(tau)])
    disp('POLYSMOOTH needs TS and TAU to be nonempty for a temporal fit.')
    tau=[];
end

% if isempty(tau)
%     disp(['POLYSMOOTH performing an order ' int2str(rho) ' spatial fit using a ' skern ' kernel.'])
% else
%     disp(['POLYSMOOTH performing an order ' int2str(rho) ' spatial fit using a ' skern ' kernel,'])
%     disp(['plus an order ' int2str(mu) ' temporal fit using a ' tkern ' kernel.'])
% end

%vsize(x,y,t,z,w,xo,yo)

%Set distance to nan for missing data; also swap infs for nans
d(~isfinite(z))=nan;
%--------------------------------------------------------------------------
%This can speed things up if you have a lot of missing data. If you don't
%sort, then you end up doing a lot of extra operations, because you can't
%truncate the matrix.
%vsize(d,x,y,t,z,w,B,tau)
%B,varstr
[d,x,y,t,z,w]=polysmooth_presort(d,x,y,t,z,w,B,tau,varstr);
%vsize(d,x,y,t,z,w)

%--------------------------------------------------------------------------
%vsize(x,y,t,z,w,xo,yo)
%size(B)
%tau,rho,mu,skern,tkern,varstr,str
%vsize(d,x,y,t,z,w)

zhat=nan*zeros(size(d,1),size(d,2),iters+1);
wo=w;
while iters+1>0
    if strcmpi(str(1:3),'loo')
        [beta,aux]=polysmooth_loop(x,y,t,z,w,xo,yo,B,tau,rho,mu,skern,tkern,varstr,geostr);
    elseif strcmpi(str(1:3),'spe')
        [beta,aux,res]=polysmooth_fast(d,x,y,t,z,w,B,tau,rho,mu,skern,tkern,varstr);
    elseif strcmpi(str(1:3),'par')
        [beta,aux]=polysmooth_parallel(d,x,y,t,z,w,B,tau,rho,mu,skern,tkern,varstr);
    else
        error(['Algorithm type ' str ' is not supported.'])
    end
    %vsize(d,beta,zhat)
    zhat(:,:,iters+1)=beta(:,:,1);
    %----------------------------------------------------------------------
    %See Cleveland (1979), p 830--831
    %Note to self, this doesn't work with my loop testing nor with parallel
    if iters>=1
        disp(['Robustification iteration #' int2str(size(zhat,3)-iters) '.'])
        s=vrep(vmedian(abs(res),3),size(res,3),3);
        delta=squared(1-squared(frac(res,6*s)));
        delta(frac(res,6*s)>1)=0;
        %length(find(delta==0))
        %This gives new robustness weights for next iteration
        %vsize(d,x,y,s,z,wo,delta,res)
        if isempty(wo)
            w=delta;
        else
            w=wo.*delta;  
        end
    end
    %----------------------------------------------------------------------
    iters=iters-1;
end
%Adjustment for complex-valued data
if ~allall(isreal(beta(:,:,1)))
    beta(isnan(real(beta(:))))=nan+sqrt(-1)*nan;
end
varargout{1}=zhat;
varargout{2}=beta;
varargout{3}=aux;

function[beta,aux]=polysmooth_parallel(ds,xs,ys,ts,zs,ws,B,tau,rho,mu,skern,tkern,varstr)

M=size(xs,1);
N=size(xs,2);
Q=sum(0:rho+1)+mu;

beta=nan*zeros(M,N,Q);
aux=nan*zeros(M,N,6);

%This way is actually slower than the fast method

parfor i=1:M
   [dsi,xsi,ysi,tsi,zsi,wsi]=vindex(ds,xs,ys,ts,zs,ws,i,1);
   if size(B,1)>1
       Bi=B(i,:);
   else
       Bi=B;
   end
   [beta(i,:,:),aux(i,:,:)]=...
       polysmooth_fast(dsi,xsi,ysi,tsi,zsi,wsi,Bi,tau,rho,mu,skern,tkern,varstr);
end

function[beta,aux,res]=polysmooth_fast(ds,xs,ys,ts,zs,ws,B,tau,rho,mu,skern,tkern,varstr)
 

%maxmax(ts(isfinite(zs))),minmin(ts(isfinite(zs)))
M=size(xs,1);
N=size(xs,2);
Q=sum(0:rho+1)+mu;
%rho,mu, Q

%vsize(ds,xs,ys,ts,zs,ws,B,tau,rho,mu,skern,tkern,varstr)
if strcmpi(varstr(1:3),'pop')
    B=polysmooth_bandwidth(B,ds,ws);  %Input B was actually P
end
%figure,jpcolor(B)

beta=nan*zeros(M,N,Q);
aux=nan*zeros(M,N,6);
 
if size(ds,3)~=0 
    %the spatial weighting kernel
    W=polysmooth_kernel(ds,ws,B,skern);
    %figure,jpcolor(W(:,:,1))    
    %multiply by the temporal weighting kernel, if requested
    if ~isempty(tau)
        %vsize(ts,ws,tau,tkern)
        W=W.*polysmooth_kernel(ts,ws,tau,tkern);
    end
    %figure,jpcolor(log10(W(:,:,1)))    

    vswap(ds,xs,ys,ts,zs,ws,W,nan,0);
    Wsum=sum(W,3);
    %length(ds)
   
    X=polysmooth_xmat(xs,ys,ts,rho,mu);
 
    XtimesW=X.*vrep(W,Q,4);
    mat=zeros(M,N,Q,Q);
 
    %vsize(X,W,XtimesW,mat,zs)
    for i=1:size(X,4)
        mat(:,:,:,i)=sum(XtimesW.*vrep(X(:,:,:,i),Q,4),3);
    end
    %vsize(XtimesW,vrep(zs,Q,4))
    vect=sum(XtimesW.*vrep(zs,Q,4),3);
    vect=permute(vect,[1 2 4 3]);
    
    C=nan*ones(M,N);
    for i=1:N
        for j=1:M
            C(j,i)=cond(squeeze(mat(j,i,:,:)));
        end
    end
    C(isinf(C))=nan;
    
    if Q==1
        vswap(mat,0,nan);
        beta=vect./mat;
    else  
        %vsize(ds,xs,ys,zs,ws,W,vect,mat,XtimesW)
        beta=matmult(matinv(mat),vect,3);
    end
    temp=vrep(permute(beta,[1 2 4 3]),size(zs,3),3);
      
    res=zs-sum(X.*temp,4);
    err=sqrt(sum(W.*squared(res),3)./sum(W,3));
   
    %Form the intercell weighted variance
    zbar=sum(W.*zs,3)./Wsum;
    V=sqrt(sum(W.*squared(zs-vrep(zbar,size(zs,3),3)),3)./Wsum);
    
    aux(:,:,1)=sum((W>0),3);       %population P
    aux(:,:,2)=B(:,:,1);           %bandwidth B
    aux(:,:,3)=err;                %rms error E
    aux(:,:,4)=sum(W,3);           %total weight
    %aux(:,:,5)=sum(W.*ds,3)./Wsum;%weighted mean distance R
    aux(:,:,5)=sum(ds,3);          %mean distance R
    aux(:,:,6)=V;                  %intercell standard deviation V
    aux(:,:,7)=C;                  %condition number C
end
 
function[beta,aux]=polysmooth_loop(x,y,t,z,w,xo,yo,B0,tau,rho,mu,skern,tkern,varstr,geostr)
%Explicit loop without pre-sorting, for testing purposes only
M=length(yo);
N=length(xo);

Q=sum(0:rho+1)+mu;

beta=nan*zeros(M,N,Q);
aux=[];

for j=1:M
    %Note that the slow method loops over the other direction for sphere
    %disp(['Polysmooth computing map for row ' int2str(j) '.'])
    for i=1:N
        if strcmpi(geostr(1:3),'sph')
            [xs,ys,ds]=latlon2xy(x,y,xo(i),yo(j));
            ds(ds>radearth*pi/2)=nan;  %Just in case
        else
            xs=x-xo(i);
            ys=y-yo(j);
            ds=sqrt(xs.^2+ys.^2);
        end
        zs=z;
        ws=w;
        ts=t;
        
        vindex(ds,xs,ys,ts,zs,ws,isfinite(ds)&isfinite(zs),1);
        
        [~,sorter]=sort(ds);
        vindex(ds,xs,ys,ts,zs,ws,sorter,1);
        %%There are some confusing 0 distances coming from +/- 90 lat
        %if sortds(2)==0
        %    figure,plot(x,y,'.')
        %    hold on, plot(xo(i),yo(j),'mo')
        %    hold on,plot(x(ds==0),y(ds==0),'g+')
        %end
        %sortds(1:20)'
        if strcmpi(varstr(1:3),'pop')
            B(j,i)=polysmooth_bandwidth(B0,ds,ws);  %Input B0 was actually P
            %B(j,i)=polysmooth_bandwidth(B0,permute(ds,[3 2 1]),permute(ws,[3 2 1]));  %Input B0 was actually P
        else
            B(j,i)=B0;
        end        
        %Correct for weird effect or repeated distances giving more than B0
        %points, which appears to be happening near the poles
%         if strcmpi(varstr(1:3),'pop')
%              if length(ds)>=B0
%                 sorter=sorter(1:B0);
%              else
%                 %If not enough data points, set all to nan
%                 ws=nan.*ws;
%              end
%         end
                
        %the spatial weighting kernel
        W=polysmooth_kernel(ds,ws,B(j,i),skern);
 
        %multiply by the temporal weighting kernel, if requested
        if ~isempty(tau)
            if ~length(tau)==1
                tauij=tau(j,i);
            else
                tauij=tau;
            end
            W=W.*polysmooth_kernel(ts,ws,tauij,tkern);
        end
        
        X=squeeze(polysmooth_xmat(xs,ys,ts,rho,mu));

        XW=conj(X')*diag(W);
        mat=XW*X;
        
        beta(j,i,:)=mat\(XW*zs); %Following Checkcode's suggestion
    end
end
aux=[];
if strcmpi(geostr(1:3),'sph')
    beta=permute(beta,[2 1 3]);
end

function[X]=polysmooth_xmat(x,y,t,rho,mu)
Q1=sum(0:rho+1);
Q=sum(0:rho+1)+mu;
X=ones(size(x,1),size(x,2),size(x,3),Q);
if rho>=1
    X(:,:,:,2)=x;
    X(:,:,:,3)=y;
end
if rho==2
    X(:,:,:,4)=frac(1,2)*x.^2;  
    X(:,:,:,5)=x.*y;
    X(:,:,:,6)=frac(1,2)*y.^2;   
end
if mu>=1
    X(:,:,:,Q1+1)=t;
end
if mu==2
    X(:,:,:,Q1+2)=frac(1,2)*t.^2;
end
function[]=polysmooth_test

tstart=tic;polysmooth_test_cartesian;toc(tstart)
tstart=tic;polysmooth_test_sphere;toc(tstart)
polysmooth_test_tangentplane;

function[]=polysmooth_test_tangentplane

%Testing tangent plane equations from Lilly and Lagerloef
load goldsnapshot
use goldsnapshot

phip=vshift(lat,-1,1);
phin=vshift(lat,1,1);
thetap=vshift(lon,-1,1);
thetan=vshift(lon,1,1);

[long,latg]=meshgrid(lon,lat);
[thetapg,phipg]=meshgrid(thetap,phip);
[thetang,phing]=meshgrid(thetan,phin);

xp=radearth*cosd(lat).*sind(thetapg-long);
xn=radearth*cosd(lat).*sind(thetang-long);
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
        reporttest('POLYSMOOTH Laplacian of linear function of y',bool)
    elseif i==2
        %figure,plot(lat,log10(abs(del2))),hold on,plot(lat,log10(abs(del2-del2hat)))
        index=find(abs(latg)<88.5);
        bool=allall(log10(abs(del2hat(index)-del2(index)))<-6);
        reporttest('POLYSMOOTH Laplacian of quadratic function of y',bool)
    elseif i==3
        %figure, plot(lon,log10(abs(del2-del2hat)./abs(del2))')
        index=find(abs(long)<178.5&isfinite(del2hat));
        %length(find(~isfinite(del2)))
        %length(find(~isfinite(del2hat)))
        bool=allall(log10(abs(del2hat(index)-del2(index))./abs(del2(index)))<-5);
        reporttest('POLYSMOOTH Laplacian of quadratic function of x',bool)
    elseif i==4
        xx=log10(abs(fx-dZdx)./abs(fx));
        index=find(isfinite(xx));
        bool=allall(xx(index)<-5);
        reporttest('POLYSMOOTH longitude gradient of GOLDSNAPSHOT',bool)
        yy=log10(abs(fy-dZdy)./abs(fy));
        index=find(isfinite(yy));
        bool=allall(yy(index)<-5);
        reporttest('POLYSMOOTH latitude gradient of GOLDSNAPSHOT',bool)
        %figure
        %subplot(1,3,1),jpcolor(log10(abs(del2))),caxis([-3 0])
        %subplot(1,3,2),jpcolor(log10(abs(del2hat))),caxis([-3 0])
        %subplot(1,3,3),jpcolor(log10(abs(del2-del2hat))),caxis([-3 0])
    end
end

function[]=polysmooth_test_cartesian

%Use peaks for testing... random assortment
[x,y,z]=peaks;
index=randperm(length(z(:)));
index=index(1:200);
[xdata,ydata,zdata]=vindex(x(:),y(:),z(:),index,1);

xo=(-3:.5:3);
yo=(-3:.6:3);

B=2;

for i=0:2
    [ds,xs,ys,zs]=twodsort(xdata,ydata,zdata,xo,yo,B);
    tic;[z1,beta1]=polysmooth(xdata,ydata,[],zdata,[],xo,yo,i,{B,'epan'},'loop');etime1=toc;
    tic;[z2,beta2]=polysmooth(ds,xs,ys,[],zs,[],i,{B,'epan'},'speed');etime2=toc;
    
    disp(['POLYSMOOTH was ' num2str(etime1./etime2) ' times faster than loop, zeroth-order fit.'])
    
    tol=1e-8;
    b1=aresame(z1,z2,tol)&&aresame(beta1,beta2,tol);
    reporttest(['POLYSMOOTH speed and loop methods are identical for order ' int2str(i) ' fit'],b1)
end

function[]=polysmooth_test_sphere
 
%Use peaks for testing... random assortment
rng(0);
[x,y,z]=peaks;
index=randperm(length(z(:)));
index=index(1:1000);

%Convert to Lat and Longitude
lon=x.*60;
lat=y.*30;

[latdata,londata,zdata]=vindex(lat(:),lon(:),z(:),index,1);
wdata=10*abs(rand(size(zdata)));  %For heavy data point test

lono=(-180:20:180);
lato=(-90:20:90);

B=2000;

for i=0:2
    [ds,xs,ys,zs,ws]=spheresort(latdata,londata,zdata,wdata,lato,lono,B);
    tic;[z1,beta1]=polysmooth(latdata,londata,[],zdata,[],lato,lono,i,{B,'epan'},'loop','sphere');etime1=toc;
    tic;[z2,beta2]=polysmooth(ds,xs,ys,[],zs,[],i,{B,'epan'});etime2=toc;
    
    disp(['POLYSMOOTH was ' num2str(etime1./etime2) ' times faster than loop, zeroth-order fit.'])
    
    tol=1e-8;
    b1=aresame(z1,z2,tol)&&aresame(beta1,beta2,tol);
    reporttest(['POLYSMOOTH speed and loop methods are identical for order ' int2str(i) ' fit'],b1)
end


for i=0:2
    [ds,xs,ys,zs,ws]=spheresort(latdata,londata,zdata,wdata,lato,lono,B);
    tic;[z1,beta1]=polysmooth(latdata,londata,[],zdata,[],lato,lono,i,{B,'epan'},'loop','sphere');etime1=toc;
    tic;[z2,beta2]=polysmooth(ds,xs,ys,[],zs,[],i,{B,'epan'});etime2=toc;
    
    disp(['POLYSMOOTH was ' num2str(etime1./etime2) ' times faster than loop, zeroth-order fit.'])
    
    tol=1e-8;
    b1=aresame(z1,z2,tol)&&aresame(beta1,beta2,tol);
    reporttest(['POLYSMOOTH speed and loop methods are identical for order ' int2str(i) ' fit with weighted data points'],b1)
end

ok, this test shows large outliers are removed in a single iteration
%znoisy=zdata+randn(size(zdata));
% znoisy=zdata;
% znoisy(500)=200;
% [ds,xs,ys,zs,znoisys]=spheresort(latdata,londata,zdata,znoisy,lato,lono,B);
% tic;[z2,beta2]=polysmooth(ds,xs,ys,[],zs,[],2,{B,'epan'});etime2=toc;
% tic;zhat=polysmooth(ds,xs,ys,[],znoisys,[],0,{B,'epan'},'robust',5);etime2=toc;

% %Zeroth-order fit at 200~km radius
% [zhat0,rbar,beta,b]=polysmooth(ds,xs,ys,[],zs,[],200,0,'sphere');
% 
% figure
% contourf(lono,lato,zhat0,[0:1/2:50]),nocontours,caxis([0 45]),
% latratio(30),[h,h2]=secondaxes;topoplot continents,
% title('Standard Deviation of SSH from TPJAOS.MAT, mapped using 200~km smoothing')
% hc=colorbar('EastOutside');axes(hc);ylabel('SSH Standard Deviation (cm)')
% set(h,'position',get(h2,'position'))
% 
% currentdir=pwd;
% cd([whichdir('jlab_license') '/figures'])
% print -dpng alongtrack_std_constant
% crop alongtrack_std.png
% cd(currentdir)


% function[beta,rbar,B,numpoints,C,Wsum]=polysmooth_fast_experimental(ds,xs,ys,zs,ws,B,R,varstr,kernstr)
% %Additional modifications here to prevent needless multiplications, but the
% %effect appears minimal.  Not worth the trouble.  But I'm leaving this
% %here as a reminder to myself that I tried it. 
% 
% M=size(xs,1);
% N=size(xs,2);
% 
% B=polysmooth_bandwidth(B,M,N,ds,ws,varstr); 
% %This is now an external function
% 
% B=vrep(B,size(xs,3),3);
% 
% P=sum(0:R+1);
% 
% beta=nan*zeros(M,N,P);
% [rbar,C,Wsum]=vzeros(M,N,nan);
% numpoints=zeros(M,N);
% 
% %Initial search to exclude right away points too far away
% ds(ds>B)=nan;
% numgood=squeeze(sum(sum(isfinite(ds),2),1));
% index=find(numgood~=0,1,'last');
% 
% if ~isempty(index)
%     vindex(ds,xs,ys,zs,ws,B,1:index,3);
%     if anyany(ws~=1)
%         W=polysmooth_kernel(dist,B,kernstr).*ws;
%         %This is now an external function
%     else
%         W=polysmooth_kernel(dist,B,kernstr);
%         %This is now an external function
%     end
%     
%     Wsum=sum(W,3,'omitnan');
%     numpoints=sum((W>0)&~isnan(W),3);
%     rbar=sum(W.*ds,3,'omitnan')./Wsum;
% 
%     vswap(ds,xs,ys,zs,W,0,nan);
%     
%     Wmat=vrep(W,P,4);
%     bool=~isnan(Wmat);
%     X=polysmooth_xmat(xs,ys,ts,Q1,Q2);    
%     XtimesW=X(bool).*Wmat(bool);  
%     %Only need to keep that part of XtimesW for which there are good values
%     
%     %vsize(X,W,XtimesW)
%     
%     mat=zeros(M,N,P,P);
%     vect=zeros(M,N,P);
% 
%     for i=1:size(X,4)
%         Xmat=vrep(X(:,:,:,i),P,4);
%         temp=nan*zeros(size(Xmat));
%         temp(bool)=XtimesW.*Xmat(bool);
%         mat(:,:,:,i)=sum(temp,3,'omitnan');
%         %mat(:,:,:,i)=sum(XtimesW.*vrep(X(:,:,:,i),P,4),3);
%    end
%     
%     zsmat=vrep(zs,P,4);
%     temp=nan*zeros(size(zsmat));
%     temp(bool)=XtimesW.*zsmat(bool);
%     vect=sum(temp,3,'omitnan');
%     %vect=sum(XtimesW.*vrep(zs,P,4),3);
%     
%     vect=permute(vect,[1 2 4 3]);
%     for i=1:N
%         for j=1:M
%             C(j,i)=cond(squeeze(mat(j,i,:,:)));
%         end
%     end
%     C(isinf(C))=nan;
%     
%     if P==1
%         vswap(mat,0,nan);
%         beta=vect./mat;
%     else 
%         beta=matmult(matinv(mat),vect,3);
%     end
% end
% B=B(:,:,1);
