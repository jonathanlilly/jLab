function[zhat,beta,aux,res]=pm_apply(varargin)
%PM_APPLY  Apply weights to data to generate a local polynomial fit.
%
%   PM_APPLY is the final stage in local polynomial fitting as implemented
%   by POLYMAP.  It applies the weighting function computed by PM_WEIGHT 
%   to observational values in order to generate a fit.
%
%   PM_APPLY is called internally by POLYMAP.  However, for large problems
%   it may be preferable to call it externally, as documented in POLYMAP.
%   _________________________________________________________________
%
%   ZHAT=PM_APPLY(ZS,AMAT,XMAT), where ZS contains sorted data values as
%   output by PM_SORT or PM_INDEX, and AMAT and XMAT as output by 
%   PM_WEIGHT, returns the local polynomial fit ZHAT.
%
%   See PM_SORT for a description of PM_WEIGHT, and PM_APPLY for a
%   description of AMAT and XMAT.
%
%   ZHAT will be an L x M matrix, where L and M are arrays of the y- and x-  
%   coordinates of the mapping grid originally input into PM_SORT.
%
%   [ZHAT,BETAHAT]=PM_APPLY(...) also returns the estimated derivatives, 
%   up to the order P specified in the input to POLYMAP or PM_WEIGHT.  
%   See POLYMAP for more details on BETAHAT.  
%
%   [ZHAT,BETAHAT,AUX]=PM_APPLY(ZS,AMAT,XMAT,WMAT,H,C), where WMAT, H, and 
%   C are output by PM_WEIGHT, also returns an L x M x 5 array AUX with 
%   some additional diagnostic information, as described in POLYMAP.
%
%   [ZHAT,BETAHAT,AUX,RES]=PM_APPLY(ZS,AMAT,XMAT,WMAT,H,C) also returns the
%   residual RES, a cell array of numeric arrays of the same size as ZS,
%   which contains the residual Z-ZHAT at all observation points.
%   __________________________________________________________________
%
%   Other options
% 
%   PM_APPLY(...,'verbose') displays the row it is currently working on.
%
%   PM_APPLY(...,'parallel') parallelizes the computation using PARFOR.
%
%   PM_APPLY(...'parallel',Nworkers) parallelizes the computation with a
%   specified number of workers.
%   __________________________________________________________________
%
%   See also POLYMAP, PM_SORT, and PM_WEIGHT.
%
%   Usage: zhat=pm_apply(zs,amat,xmat);
%          [zhat,beta]=pm_apply(zs,amat,xmat);
%          [zhat,beta,aux]=pm_apply(zs,amat,xmat,wmat,H,C);
%          [zhat,beta,aux,res]=pm_apply(zs,amat,xmat,wmat,H,C);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2022 J.M. Lilly --- type 'help jlab_license' for details
 
% if strcmp(varargin{1}, '--t')
%     pm_apply_test,return
% end

str='series';
verbstr='quiet';
Nworkers=[];
Pmin=[];

%First parse the string arguments
for i=1:4
    if ischar(varargin{end})
        tempstr=varargin{end};
        if strcmpi(tempstr(1:3),'ver')||strcmpi(tempstr(1:3),'qui')
            verbstr=tempstr;
        elseif strcmpi(tempstr(1:3),'ser')||strcmpi(tempstr(1:3),'par')||strcmpi(tempstr(1:3),'nol')
            str=tempstr;
        end
        varargin=varargin(1:end-1);
    elseif ~ischar(varargin{end})&&ischar(varargin{end-1})
        tempstr=varargin{end-1};
        if strcmpi(tempstr(1:3),'par')
            str=tempstr;
            Nworkers=varargin{end};
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

%ds=varargin{1};
zs=varargin{1};
amat=varargin{2};
xmat=varargin{3};

L=length(zs);
M=size(zs{1},2);
N=size(zs{1},1);
Q=size(xmat{1},2);

if length(varargin)>3
    wmat=varargin{4};
    H=varargin{5};
    C=varargin{6};
else
    wmat=cell(length(amat),1);
    H=zeros(L,M);
    C=zeros(L,M);
end

if  strcmpi(str(1:3),'nol')
    z=zeros(N,1,L,M);
    A=zeros(Q,N,L,M);
    X=zeros(N,Q,L,M);
%    W=zeros(N,1,L,M);
    for i=1:length(zs)
        z(:,:,i,:)=zs{i};
        X(:,:,i,:)=xmat{i};
        A(:,:,i,:)=amat{i};
    end

tic;    [zhat,beta,aux,res]=pm_apply_one_noloop(z,A,X,[],H,C,nargout);toc

else
 
    zhat=nan*zeros(L,M,size(zs{1},3));
    beta=nan*zeros(L,M,Q,size(zs{1},3));
    aux=nan*zeros(L,M,5,size(zs{1},3));
    res=zs;
    Nout=nargout;

    if ~isempty(zs)
        if  strcmpi(str(1:3),'ser')
            for i=1:L
                if strcmpi(verbstr(1:3),'ver')
                    disp(['PM_APPLY computing fit for row ' int2str(i) ' of ' int2str(L) '.'])
                end
                if ~isempty(zs{i})
%                    [zhat(i,:,:),beta(i,:,:,:),aux(i,:,:,:),res{i}]=pm_apply_one(zs{i},amat{i},xmat{i},wmat{i},H(i,:),C(i,:),Nout);
                    zhat(i,:,:)=pm_apply_one(zs{i},amat{i},xmat{i},wmat{i},H(i,:),C(i,:),Nout);
                end
            end
        elseif strcmpi(str(1:3),'par')
            parfor i=1:L
                if strcmpi(verbstr(1:3),'ver')
                    disp(['PM_APPLY computing fit for row ' int2str(i) ' of ' int2str(L) '.'])
                end
                if ~isempty(zs{i})
                    [zhat(i,:,:),beta(i,:,:,:),aux(i,:,:,:),res{i}]=pm_apply_one(zs{i},amat{i},xmat{i},wmat{i},H(i,:),C(i,:),Nout);
                end
            end
        end
    end
end
%vsize(zs{1},res{1})

function[zhat,beta,aux,res]=pm_apply_one(zs,amat,xmat,wmat,H,C,Nout)

%amat is Q x N x M
%wmat is 1 x N x M
%xmat is N x Q x M

%initialize arrays in case we don't output them
zhat=nan*zeros(1,size(xmat,3),size(zs,3));
beta=nan*zeros(size(amat,1),1,size(amat,3),size(zs,3));
aux=nan*zeros(1,size(amat,3),5);
res=nan*zs;
%vsize(zhat,beta,aux,res)

if size(xmat,1)>0
    %truncate rows of ds and zs to exclude those for which the kernel vanishes
    %vsize(amat,wmat,xmat,zs)
    %ds=ds(1:size(xmat,1),:);
  %  vsize(zs,xmat)

%    zs=zs(1:size(xmat,1),:);
    
    %asum=squeeze(sum(amat(1,:,:),3));
    %ii=find(asum==0,1,'first')
    %size(amat)
    
   % vsize(zs)

    %don't multiply land columns
    bool=~isnan(zs(1,:,1));
%    anyany(~isfinite(zs(:,bool)))

    if Nout == 1
        %we only need to compute z, not all of beta
        %zhat=permute(pagemtimes(amat(1,:,:),permute(zs,[1 3 2])),[2 3 1]);
%        zhat(:,bool,:)=permute(pagemtimes(amat(1,:,bool),permute(zs(:,bool,:),[1 3 2])),[2 3 1]);
        zhat(:,bool,:)=permute(pagemtimes(amat(1,:,bool),permute(zs(:,bool,:),[1 3 2])),[3 2 1]);
        %vsize((amat(1,:,:)),permute(zs,[1 3 2]),pagemtimes(amat(1,:,:),permute(zs,[1 3 2])))
    else
        %beta is Q x 1 x M
        beta(:,:,bool)=pagemtimes(amat(:,:,bool),permute(zs(:,bool,:),[1 3 2]));
        %vsize(beta,zs,amat)
        zhat=beta(1,:,:);
    end

    %vsize(amat,wmat,xmat,beta,zs,res)
    %vsize(zhat,beta,aux)
    %vsize(zs,aux)
    if Nout>2
        res=zs-squeeze(pagemtimes(xmat,beta));

        wmat=permute(wmat,[2 3 1]);
        sumW=sum(wmat,1);
        %wmat is now N x M
        aux(:,:,1)=squeeze(sum(wmat>0,1));                    %population N
        aux(:,:,2)=H;                                         %bandwidth H
        aux(:,:,3)=sum(wmat,1);                               %total weight
        aux(:,:,4)=C;                                         %condition number
        aux(:,:,5)=sqrt(sum(wmat.*squared(res),1,'omitnan')); %rms error E
        %distance to the last nonzero entry
        %ds(isnan(wmat)|wmat==0)=nan;
        %aux(:,:,2)=max(ds,[],1,'omitnan');
        %aux(:,:,5)=sum(wmat.*ds,1,'omitnan')./sumW; %weighted mean distance R
        %not sure I need this last one?
        %zbar=sum(wmat.*zs,1,'omitnan')./sumW;      %this is just the zeroth-order fit
        %zbar=vrep(zbar,size(zs{i},1),1);
        %aux(:,:,6)=sqrt(sum(W.*squared(zs{i}-zbar),1)./sumW); %intercell standard deviation V

        %put back in missing data from zs
        %res is N x M
        res(wasnan)=nan;
        %res(wasinf)=inf;
        %aux(~isfinite(aux))=inf;
        aux(:,~bool,:)=nan;
    end

    %vsize(zs,zhat)
end

%size(beta)
%beta=permute(beta,[2 3 1]);
%vsize(zhat,beta,aux)

%Adjustment for complex-valued data
if ~allall(isreal(beta(:,:,1)))
    beta(isnan(real(beta(:))))=nan+sqrt(-1)*nan;
    % beta(isinf(real(beta(:))))=inf+sqrt(-1)*inf;
end



function[zhat,beta,aux,res]=pm_apply_one_noloop(z,A,X,W,H,C,Nout)


%amat is Q x N x M
%wmat is 1 x N x M
%xmat is N x Q x M

% %initialize arrays in case we don't output them
% zhat=nan*zeros(1,size(xmat,3));
beta=nan*zeros(size(A,1),1,size(A,3));
aux=nan*zeros(1,size(A,3),5);
res=nan*z;
% %vsize(zhat,beta,aux,res)

%if size(xmat,1)>0
    %truncate rows of ds and zs to exclude those for which the kernel vanishes
    %vsize(amat,wmat,xmat,zs)
    %ds=ds(1:size(xmat,1),:);
  %  vsize(zs,xmat)

  %
 %   zs=zs(1:size(xmat,1),:);

   % vsize(zs)

    %swap nans for 0's in zs
    
    %?
    boolnan=squeeze(isnan(z(1,1,:,:)));%figure,plot(boolnan)
    wasnan=isnan(z);
    z(wasnan)=0;

 %   anyany(~isfinite(A))

    % boolinf=isinf(zs(1,:));
    % boolzero=(squeeze(amat(1,1,:))==0)';
    % wasinf=isinf(zs);
    % zs(wasnan)=0;
    % zs(wasinf)=0;

    if Nout == 1
        %we only need to compute z, not all of beta
        zhat=squeeze(pagemtimes(A(1,:,:,:),z));
    else
        %beta is Q x 1 x M
        beta=pagemtimes(A,z);
        zhat=squeeze(beta(1,:,:,:));
    end

    %vsize(amat,wmat,xmat,beta,zs,res)
    %vsize(zhat,beta,aux)
    %vsize(zs,aux)
%     if Nout>2
%         res=z-xmat,beta;
% 
%         wmat=permute(wmat,[2 3 1]);
%         sumW=sum(wmat,1);
%         %wmat is now N x M
%         aux(:,:,1)=squeeze(sum(wmat>0,1));                    %population N
%         aux(:,:,2)=H;                                         %bandwidth H
%         aux(:,:,3)=sum(wmat,1);                               %total weight
%         aux(:,:,4)=C;                                         %condition number
%         aux(:,:,5)=sqrt(sum(wmat.*squared(res),1,'omitnan')); %rms error E
%         %distance to the last nonzero entry
%         %ds(isnan(wmat)|wmat==0)=nan;
%         %aux(:,:,2)=max(ds,[],1,'omitnan');
%         %aux(:,:,5)=sum(wmat.*ds,1,'omitnan')./sumW; %weighted mean distance R
%         %not sure I need this last one?
%         %zbar=sum(wmat.*zs,1,'omitnan')./sumW;      %this is just the zeroth-order fit
%         %zbar=vrep(zbar,size(zs{i},1),1);
%         %aux(:,:,6)=sqrt(sum(W.*squared(zs{i}-zbar),1)./sumW); %intercell standard deviation V
% 
%         %put back in missing data from zs
%         %res is N x M
%         res(wasnan)=nan;
%         %res(wasinf)=inf;
%         %aux(~isfinite(aux))=inf;
%         aux(:,boolnan,:)=nan;
%     end

    %vsize(zs,zhat)
%    zhat(:,boolnan)=nan;
%    beta(:,:,boolnan)=nan;
%end

%put back in missing data from zs
%vsize(zhat,beta,boolinf,boolzero)
%zhat(:,boolinf|boolzero)=inf;
%beta(:,:,boolinf|boolzero)=inf;
%vsize(zhat,beta)

% beta=permute(beta,[2 3 1]);
% %vsize(zhat,beta,aux)
% 
% %Adjustment for complex-valued data
% if ~allall(isreal(beta(:,:,1)))
%     beta(isnan(real(beta(:))))=nan+sqrt(-1)*nan;
%     % beta(isinf(real(beta(:))))=inf+sqrt(-1)*inf;
% end

