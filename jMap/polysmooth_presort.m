function[d,x,y,t,z,w]=polysmooth_presort(varargin)
%POLYSMOOTH_PRESORT  Sort arguments to POLYSMOOTH in case of missing data.
%
%   POLYSMOOTH_PRESORT is a low-level function called by POLYSMOOTH. 
%
%   [DS,XS,YS,TS,ZS,WS]=POLYSMOOTH_PRESORT(DS,XS,YS,TS,ZS,WS,P,TAU,STR),
%   where all the input fields are in the form output by TWODSORT or 
%   SPHERESORT, performs several prelimary processing steps. 
%
%   First, if the temporal bandwidth TAU is nonempty, any datapoints
%   outside of temporal window ABS(TS./TAU)>1 are set to NaNs.
%
%   After this, the data field, ZS, is checked to see if it contains any 
%   NaNs. If it does, then all fields are sorted by distance DS along 
%   third dimension, thus moving any interior NaN value to the bottom.
%
%   In the case that STR='pop', indicating a fixed population algorithm,
%   then finally all input fields are additionally truncated to length P 
%   along their third dimension.
%
%   See also POLYSMOOTH.
%
%   Usage: [ds,xs,ys,zs,ws]=polysmooth_presort(ds,xs,ys,zs,ws,P,tau,str);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2018 J.M. Lilly --- type 'help jlab_license' for details

t=[];
z=[];
w=[];
B=[];
tau=[];
varstr='bandwidth';

d=varargin{1};
x=varargin{2};
y=varargin{3};
if nargin>=4
    t=varargin{4};
end
if nargin>=5
    z=varargin{5};
end
if nargin>=6
    w=varargin{6};
end
if nargin>=7
    B=varargin{7};
end
if nargin>=8
    tau=varargin{8};
end
if nargin>=9
    varstr=varargin{9};
end
%figure,jpcolor(squeeze(z))

if ~isempty(tau)
    if length(tau(:))~=1
        tau=vrep(tau,size(z,3),3);
    end
    z(abs(t./tau)>1)=nan;
    %d(abs(t./tau)>1)=nan;
end

%figure,jpcolor(squeeze(z))
%d(~isfinite(z))=nan;
%x(~isfinite(z))=nan;
%y(~isfinite(z))=nan;

if ~isempty(d)
    if ~strcmpi(varstr(1:3),'pop')
        if ~isempty(B)
            if length(B(:))~=1
                B=vrep(B,size(z,3),3);
            end
            z(abs(d./B)>1)=nan;
        end
    end
    %figure,jpcolor(squeeze(z))

    if ~isempty(z)
        z(~isfinite(d))=nan;  %make sure pattern of nans matches
        d(~isfinite(z))=nan;  %make sure pattern of nans matches
    end
    
    %bool will be true if current point is finite, but previous is not
    bool1=isnan(d);
    bool2=isnan(vshift(d,-1,3));bool2(:,:,1)=false;
    bool=(sum(~bool1&bool2,3)>0)&~(sum(bool1,3)==size(bool1,3));
    
    if anyany(bool)
        disp('POLYSMOOTH_PRESORT detecting interior NaNs; sorting to compress these.')
        [d,kk]=sort(d,3);
        ii=vrep([1:size(d,1)]',[size(d,2) size(d,3)],[2 3]);
        jj=vrep([1:size(d,2)],[size(d,1) size(d,3)],[1 3]);
        index=sub2ind(size(d),ii,jj,kk);
        %Can't use vindex here because i want to preserve the size
        x=x(index);
        y=y(index);
        if ~isempty(z),z=z(index);end
        if ~isempty(t),t=t(index);end
        if ~isempty(w),w=w(index);end
        %vindex(x,y,t,z,w,index,3);
        disp('Sorting complete.')
    end
    %1,vsize(d,x,y,t,z,w)
    
    if strcmpi(varstr(1:3),'pop')
        if ~isempty(B)
            vindex(d,x,y,t,z,w,1:min(maxmax(B),size(d,3)),3);
        end
    end
    %2,vsize(d,x,y,t,z,w)
    %figure,jpcolor(squeeze(d))
    numgood=squeeze(sum(sum(isfinite(d),2),1));
    %figure,plot(numgood)
    index=find(numgood~=0,1,'last');
    %index
    vindex(d,x,y,t,z,w,1:index,3);   
    %3,vsize(d,x,y,t,z,w)
end