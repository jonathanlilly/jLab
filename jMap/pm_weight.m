function[index,A,X,W,H,C]=pm_weight(varargin)
%PM_WEIGHT  Weighting function for local polynomial fitting. 
%
%   PM_WEIGHT computes the weighting function and various related 
%   quantities used in local polynomial fitting.  
%
%   PM_WEIGHT is called internally by POLYMAP.  However, for large problems
%   it may be preferable to call it externally, as documented in POLYMAP.
%   _________________________________________________________________
%
%   [INDEX,AMAT,XMAT,WMAT]=...
%                        PM_WEIGHT(DS,XS,YS,INDEXS,[],[],{P,H,ALPHA,BETA}), 
%
%   where DS, XS, YS and INDEXS are output by PM_SORT, returns quantities
%   needed for a Pth order fit using the kernel specified by ALPHA and BETA.
%
%   The kernel is a beta kernel of the form KAPPA*(1-(DS/H)^ALPHA)^BETA, 
%   where KAPPA is a normalization coefficient depending on ALPHA and BETA.
%
%   DS, XS, YS, and INDEXS describe the sorted distribution of data points 
%   relative to an L x M grid.  These are length L cell arrays, with each 
%   cell containing M columns and a variable number of rows, denoted NO_l,
%   that depends on the L index.  
%
%   Note that H, ALPHA, and BETA may all either be scalars or L x M
%   matrices.  The fit order P must be a scalar. 
%
%   The solution to the least squares problem in local polynomial fitting
%   is of the form
% 
%         BETAHAT = (X^T W X)^(-1) X^T W Z  
%
%   where X and W are matrices defined shortly and Z is an N x 1 array 
%   of observations.  BETAHAT is a vector of the estimated field together
%   with its first P derivatives.  BETAHAT is Q x 1, where Q = 1, 3, and 6
%   for P = 0, 1, and 2, respectively. 
%
%   Here X is an N x Q matrix of powers of coordinate deviations arising 
%   from a Pth order Taylor expansion and W is an N x N matrix of weights. 
%   For convenience we define Q X N matrix A as A = (X^T W X)^(-1) * X^T W, 
%   such that can write simply write BETAHAT = A Z.
%
%   The output of PM_WEIGHT are each length L cell arrays of numerical
%   arrays, the lth element of which has the size
%
%       INDEX --- N_l x M
%        AMAT ---  Q  x N_l x M
%        WMAT ---  1  x N_l x M
%        XMAT --- N_l x  Q  x M
%
%   where N_l is less that or equal to the original size NO_L.  This 
%   dimension is reduced from the original such that any points further
%   away than the bandwidth H are excluded.
%
%   INDEX is a truncated version of INDEXS, while AMAT, WMAT, and XMAT 
%   contain respectively the A, W, and X matrices for each cell. 
%   __________________________________________________________________
%
%   Weighted data points
%
%   PM_WEIGHT(DS,XS,YS,[],WS,{P,H,ALPHA,BETA}) also applies the additional
%   nonnegative weights WS to the data points.  
%   
%   This is useful in at least three circumstances: (i) weighting data 
%   points according to their estimated significance, (ii) weighting as a 
%   part of an iterative robustification scheme, or (iii) weighting
%   aggregated data points according to the number of observations 
%   represented, as a way of downscaling large numbers of observations.  
%   __________________________________________________________________
%
%   Temporal fit
%
%   PM_WEIGHT(DS,XS,YS,TS,[],{P,H,ALPHA,BETA},{MU,TAU,ALPHAT,BETAT}) and
%   PM_WEIGHT(DS,XS,YS,TS,WS,{P,H,ALPHA,BETA},{MU,TAU,ALPHAT,BETAT}) also
%   incorporates an order MU fit in the time domain.  MU may be 1 or 2.
% 
%   Here TS represents the times of the observations, in the same units as
%   the temporal bandwidth TAU, with the center of the weighting function 
%   corresponding to time zero.  
%
%   ALPHAT and BETAT describe the properties of the temporal kernel. 
%
%   In output array sizes, Q is replaced with QTILDE=Q+MU.
%   __________________________________________________________________
%
%   Additional output arguments
%
%   [AMAT,XMAT,WMAT,H,C]=PM_WEIGHT(...) also outputs two L x M matrices.
%
%   H is the matrix of bandwidths used to construct the weighting function.
%   This may differ from the H that is input when the minimum population
%   is specified, or when the 'population' flag is input; see below.
%
%   C contains the condition number of the inverted Q x Q matrix [X^T W X] 
%   at each grid point, a potentially useful diagnostic.
%   __________________________________________________________________
%
%   Other options
%
%   PM_WEIGHT supports a number of other options, as documented more 
%   completely in POLYMAP.
%
%   PM_WEIGHT(...,'minpop',NMIN) uses NMIN for the minimum population.
%
%   PM_WEIGHT(...,{P,HO,ALPHA,BETA},...,'population',N) specifies a fixed
%   population fit with population N and maximum bandwidth HO.  N can be a
%   scalar or an L x M matrix.
% 
%   PM_WEIGHT(...,'verbose') displays the row it is currently working on.
%
%   PM_WEIGHT(...,'parallel') parallelizes the computation using PARFOR.
%
%   PM_WEIGHT(...'parallel',Nworkers) parallelizes the computation with a
%   specified number of workers.
%   __________________________________________________________________
%
%   See also POLYMAP, PM_SORT, and PM_APPLY.
%
%  'pm_weight --t' runs a test.
%
%   Usage: [index,amat,xmat,wmat,H,C]=...
%                               pm_weight(ds,xs,ys,indexs,[],[],{P,H,2,1});
%          [index,amat,xmat,wmat,H,C]=...
%                               pm_weight(ds,xs,ys,indexs,ts,ws,{P,H,2,1});
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2022 J.M. Lilly --- type 'help jlab_license' for details
 

%   Note that one may also call PM_WEIGHT with ZS in the place of INDEXS:
%
%   [Z,AMAT,XMAT,WMAT]=PM_WEIGHT(DS,XS,YS,ZS,[],[],{P,H,ALPHA,BETA})
%
%   where Z is now a truncated version of ZS.  This version is used inside 
%   POLYMAP, but when calling PM_WEIGHT externally the useful thing to do
%   is to call it with 

if strcmp(varargin{1}, '--t')
    pm_bandwidth_test,return
end

str='series';
varstr='bandwidth';
verbstr='quiet';
Nworkers=[];
Pmin=[];
Npop=nan;

%First parse the string arguments
for i=1:6
    if ischar(varargin{end})
        tempstr=varargin{end};
        if strcmpi(tempstr(1:3),'ver')||strcmpi(tempstr(1:3),'qui')
            verbstr=tempstr;
        elseif strcmpi(tempstr(1:3),'ser')||strcmpi(tempstr(1:3),'par')
            str=tempstr;
        end
        varargin=varargin(1:end-1);
    elseif ~ischar(varargin{end})&&ischar(varargin{end-1})
        tempstr=varargin{end-1};
        if strcmpi(tempstr(1:3),'par')
            str=tempstr;
            Nworkers=varargin{end};
        elseif strcmpi(tempstr(1:3),'pop')
            varstr=tempstr;
            Npop=varargin{end};
        elseif strcmpi(tempstr(1:3),'min')
            Pmin=varargin{end};
        end
        varargin=varargin(1:end-2);
    end
end

pool = gcp('nocreate');
if strcmpi(str(1:3),'par')
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

ds=varargin{1};
xs=varargin{2};
ys=varargin{3};
indexs=varargin{4};
ts=varargin{5};
ws=varargin{6};
sarg=varargin{7};
if length(varargin)==8
    targ=varargin{8};
else
    targ=[];
end

L=length(ds);
M=size(ds{1},2);

%spatial fitting arguments
rho=sarg{1};
H=sarg{2};
salpha=sarg{3};
sbeta=sarg{4};
if numel(H)==1,      H=H+zeros(L,M);           end
if numel(salpha)==1, salpha=salpha+zeros(L,M); end
if numel(sbeta)==1,  sbeta=sbeta+zeros(L,M);   end
if numel(Npop)==1,   Npop=Npop+zeros(L,M);     end

%temporal fitting arguments
if isempty(targ)
    mu=[];
    tau=zeros(L,M);
    talpha=cell(L,1);
    tbeta=cell(L,1);    
else
    mu=targ{1};
    tau=targ{2};
    talpha=targ{3};
    tbeta=targ{4};
    if numel(tau)==1,    tau=tau+zeros(L,M);       end
    if numel(talpha)==1, talpha=talpha+zeros(L,M); end
    if numel(tbeta)==1,  tbeta=tbeta+zeros(L,M);   end
end

%convert empties to cell array of empties
if isempty(ts)
    ts=cell(size(ds)); 
end
if isempty(ws)
    ws=cell(size(ds)); 
end

A=cell(length(ds),1);
W=cell(length(ds),1);
X=cell(length(ds),1);

%initialize empty sizes for future reference, in case no data
for i=1:length(A)
    if isempty(mu)
        Q=sum(0:rho+1);
    else
        Q=sum(0:rho+1)+mu;
    end
    index{i}=zeros(0,M);
    A{i}=zeros(Q,0,M);
    W{i}=zeros(1,0,M);
    X{i}=zeros(0,Q,M);
end
C=zeros(size(H));
if ~isempty(ds)
    if  strcmpi(str(1:3),'ser')
        for i=1:L
            if strcmpi(verbstr(1:3),'ver')
                disp(['PM_WEIGHT computing equivalent kernel for row ' int2str(i) ' of ' int2str(L) '.'])
            end
            if ~isempty(ds{i})
                [index{i},A{i},W{i},X{i},H(i,:),C(i,:)]=...
                    pm_weight_one(ds{i},xs{i},ys{i},indexs{i},ts{i},ws{i},H(i,:),tau(i,:),rho,mu,salpha(i,:),sbeta(i,:),talpha(i,:),sbeta(i,:),Npop(i,:),Pmin);
            end
        end
    elseif strcmpi(str(1:3),'par')
        parfor i=1:L
            if strcmpi(verbstr(1:3),'ver')
                disp(['PM_WEIGHT computing equivalent kernel for row ' int2str(i) ' of ' int2str(L) '.'])
            end
            if ~isempty(ds{i})
                [index{i},A{i},W{i},X{i},H(i,:),C(i,:)]=...
                    pm_weight_one(ds{i},xs{i},ys{i},indexs{i},ts{i},ws{i},H(i,:),tau(i,:),rho,mu,salpha(i,:),sbeta(i,:),talpha(i,:),sbeta(i,:),Npop(i,:),Pmin);
            end
        end
    end
end

function[indexs,A,W,X,H,C]=pm_weight_one(ds,xs,ys,indexs,ts,ws,H,tau,rho,mu,salpha,sbeta,talpha,tbeta,Npop,Pmin)

%don't set distance to 0; we need to know the distance
[xs,ys]=vswap(xs,ys,nan,0);
%[xs,ys]=vswap(xs,ys,inf,0);
%length(find(isnan(vcolon(zs{i}))))
%ds{i}=vswap(ds{i},nan,inf);
%vsize(ds,xs,ys,ts,ws,H,tau,rho,mu,salpha,sbeta)

if ~all(~isfinite(Npop)) 
    HO=maxmax(H);
    H=pm_bandwidth(Npop,ds,ws);
    H(H>HO|isnan(H))=HO;
end
%adjust bandwidth to guarantee Pmin points, if input
if ~isempty(Pmin)
    if size(ds,1)>Pmin
        dsPmin=ds(Pmin,:)+1e-9;%set bandwidth to just larger than this point
        %note: can't use Pmin+1 because sometimes the distances are repeated
        H(H<dsPmin)=dsPmin(H<dsPmin);
    end
end

%now form kernel and strip unneeded rows from ds etc.
[W,ds,xs,ys,indexs,ts]=pm_kernel(ds,xs,ys,indexs,ts,ws,H,salpha,sbeta,tau,talpha,tbeta);
%W is N x M

%vsize(zs{i},W),return
if rho == 0 && isempty(mu) %simplification and speedup for constant case
    A=W./sum(W,1);
    %A is N x M
    A(~isfinite(A))=0;

    C=ones(1,size(W,2));
    A=permute(A,[3 1 2]);
    W=permute(W,[3 1 2]);
    X=ones(size(ds,1),1,size(ds,2));
    
    %length(find(isnan(W)))
    %length(find(isinf(W)))
    %vsize(A,W,X)
else
    %truncate xs, ys, and ts according to the number of rows of ds
    XT=pm_weight_xmatT(xs,ys,ts,rho,mu);
    %XT is Q x N x M
    
    W=permute(W,[3  1 2]);
    %XT is Q x N x M
    %W  is 1 x N x M

    XTW=XT.*W;%Since these are compatible, I can multiply elementwise
    %XTW is Q x N x M

    XTWX=pagemtimes(XTW,'none',XT,'transpose');
    %XTWX is Q x Q x M

    A=pagemtimes(matinv(XTWX),XTW);
    %A is Q x N x M
   % vsize(ds,A)
    
    A(~isfinite(A))=0;

    %set all bad data to infs, then masked-out data to NaNs
%    A(~isfinite(A))=inf;
%    A(:,:,isnan(ds(1,:)))=nan;
    
    C=zeros(1,size(XTWX,3));
    for i=1:size(XTWX,3)
        C(:,i)=cond(squeeze(XTWX(:,:,i)));
    end
    X=permute(XT,[2 1 3 4]);
end


function[X]=pm_weight_xmatT(x,y,t,rho,mu)

Q1=sum(0:rho+1);
Q=Q1;
if ~isempty(mu)
    Q=Q+mu;
end
 
X=ones(Q,size(x,1),size(x,2));
if rho>=1
    X(2,:,:)=x;
    X(3,:,:)=y;
end
if rho>=2
    X(4,:,:)=frac(1,2)*x.^2;  
    X(5,:,:)=x.*y;
    X(6,:,:)=frac(1,2)*y.^2;   
end
if mu>=1
    X(Q1+1,:,:)=t;
end
if mu==2
     X(Q1+2,:,:)=frac(1,2)*t.^2;
end

function[W,ds,xs,ys,indexs,ts,ws]=pm_kernel(ds,xs,ys,indexs,ts,ws,H,al,be,tau,tal,tbe)
%PM_KERNEL  Returns the weighting kernel employed by POLYMAP.
%
%   PM_KERNEL is a low-level function called by POLYMAP. 
%
%   [W,DS,XS,YS,TS]=PM_KERNEL(DS,XS,YS,TS,WS,H,AL,BE,TAU,TAL,TBE) where DS
%   is an array of distances and H is the bandwidth, returns the value
%   of the weighting kernel W.
% 
%   W is the generalized beta kernel proportional to (1-(DS/H)^AL)^BE if TS
%   is empty, or to (1-(DS/H)^AL)^BE * (1-(TS/TAU)^TAL)^TBE otherwise.
%   This is multiplied by the additional weights WS if WS is nonempty.
%
%   All input arguments after DS are either scalars, or numeric arrays of
%   length SIZE(DS,2), i.e. the same number of entries as DS has columns.
%
%   In the output, rows of DS, XS, YS, and TS are truncated, as these are 
%   needed subsequently within PM_WEIGHT.  The last row corresponds to the 
%   last row of DS for which any element of DS is less than or equal to H.
%
%   See also POLYMAP.
%
%   Usage: [w,ds,xs,ys,ts]=pm_kernel(ds,xs,ys,ts,ws,H,al,be,tau,tal,tbe);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2022 J.M. Lilly --- type 'help jlab_license' for details

ds=ds./(H(:)');  %convert H to a row vector

%find last row for which all ds<1
minds=min(ds,[],2);
ii=find(minds<1,1,'last');

%truncate ds and ws 
ds=ds(1:ii,:);
xs=xs(1:ii,:);
ys=ys(1:ii,:);
indexs=indexs(1:ii,:);
if ~isempty(ws)
    ws=ws(1:ii,:);
end
if ~isempty(ts)
    ts=ts(1:ii,:);
end

W=zeros(size(ds));

al=vrep(al(:)',size(ds,1),1); %convert to a row vector and replicate
be=vrep(be(:)',size(ds,1),1); %convert to a row vector and replicate

bool=(ds<=1);

kappa=frac(1,2*pi).*frac(al(bool),beta(2./al(bool),be(bool)+1));
W(bool)=kappa.*(1-ds(bool).^al(bool)).^be(bool);

%     older code for deprecated custom kernel interpolation
%     K=kern(:);
%     W(bool)=interp1((0:length(K)-1)'./(length(K)-1),K,ds);

%multiply by the temporal weighting kernel, if requested
if ~isempty(ts)
    ts=ts./(tau(:)');  %convert tau to a row vector
    tal=vrep(tal(:)',size(ts,1),1); %convert to a row vector and replicate
    tbe=vrep(tbe(:)',size(ts,1),1); %convert to a row vector and replicate
    kappa=frac(1,2).*frac(tal(bool),beta(1./tal(bool),tbe(bool)+1));
    %kappa=(al/2)/beta(1/al,be+1) for d=1
    W(bool)=W(bool).*kappa.*(1-ts(bool).^tal(bool)).^tbe(bool);
end

if ~isempty(ws)
    W(bool)=W(bool).*ws(bool);
end

W(~isfinite(W))=0;
%W=vswap(W,nan,0);
%W=vswap(W,inf,0);


function[H]=pm_bandwidth(P,ds,ws)
%PM_BANDWIDTH  Determine bandwidth given population for POLYMAP.
%
%   PM_BANDWIDTH is a low-level function called internally by POLYMAP.  
%   It's unlikely that you would need to call it directly. 
% 
%   H=PM_BANDWIDTH(NPOP,DS) returns an L x M matrix of spatially-varying
%   bandwidths H implied by the population NPOP. DS is the distance output 
%   by PM_SORT and has L cells each containing a matrix with M columns.
%   
%   H is the bandwidth that you would need to have population NPOP at each
%   grid point. 
%
%   H=PM_BANDWIDTH(NPOP,DS,WS) is similar, except that each measurement
%   contributes an amount to the effective population given by WS.
%
%   NPOP may be a scalar or matrix of size L x M, the size of the mapping
%   grid input to PM_SORT.  The output H will always be L x M.
%
%   See also POLYMAP.
%
%   Usage: H=pm_bandwidth(npop,ds,ws);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2018--2022 J.M. Lilly --- type 'help jlab_license' for details

 
if strcmp(P, '--t')
    pm_bandwidth_test 
    return
end

%Convert number of points input to equivalent bandwidth
P=ceil(P);%What I've actually input is the number of points we want to keep

P(P<1)=1;% set minimum population to one

%if minmin(P)<=1
%    error('With fixed population algorithm, population must be greater than one.')
%end

if isempty(ws)
    ws=ones(size(ds));
    ws(isnan(ds))=nan;
    %ws(isinf(ds))=inf;
end
cumpop=cumsum(ws);
ii=cumsum(ones(size(ws)));

%vsize(cumpop,P)

ii(cumpop<=P)=inf;
ii=min(ii,[],1);
H=nan*zeros(size(ii));
colindex=find(isfinite(ii));
index=sub2ind(size(ds),ii(colindex),colindex);
H(colindex)=ds(index);

%Note you may have more than the expected number of points when there are
%exactly repeated distances

function[]=pm_bandwidth_test

rng(1);
[x,y,z]=peaks;
xdata=6*rand(200,1)-3;
ydata=6*rand(200,1)-3;
zdata=interp2(x,y,z,xdata,ydata);

xo=(-3:.125:3);
yo=(-3:.125:3);
[xg,yg]=meshgrid(xo,yo);

[ds,xs,ys,zs]=pm_sort(xdata,ydata,zdata,xo,yo,20);
[zs,amat,xmat,wmat,H,C]=pm_weight(ds,xs,ys,zs,[],[],{0,50,2,1},'population',10);
[zhat,beta,aux]=pm_apply(zs,amat,xmat,wmat,H,C);

reporttest(['PM_WEIGHT population flag without weights'],allall(aux(:,:,1)==10))

ws=ds;
for i=1:length(ds)
    ws{i}=1/2+zeros(size(ws{i}));
end

[ds,xs,ys,zs]=pm_sort(xdata,ydata,zdata,xo,yo,20);
[zs,amat,xmat,wmat,H,C]=pm_weight(ds,xs,ys,zs,[],ws,{0,50,2,1},'population',10);
[zhat,beta,aux]=pm_apply(zs,amat,xmat,wmat,H,C);

reporttest(['PM_WEIGHT population flag with weights'],allall(aux(:,:,1)==20))

