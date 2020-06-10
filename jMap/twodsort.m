function[varargout]=twodsort(varargin)
%TWODSORT  Distances from data points to nearby grid points.
%
%   [DS,XS,YS]=TWODSORT(X,Y,XO,YO,CUTOFF) returns sorted distances D 
%   between data points at locations X,Y and grid points at XO,YO.
%
%   X and Y are arrays of the same size into data point locations. XO and
%   YO are arrays of length M and N, say, specifying the bin center
%   locations of an M x N matrix of grid points, i.e.
%
%      XO= [XO_1 XO_2 ... XO_N]      XO =  [YO_1;    
%                                           YO_2; 
%                                            ...
%                                           YO_M]
%
%   CUTOFF is the maximum distance to be included in the output arrays.
%
%   The output arrays are M numerical arrays arranged as a length M cell
%   array.  That is, there is one cell per element of Y0. Each numerical 
%   array has N columns, i.e., the number of elements of X0, with the 
%   number of rows varying between arrays.  
%
%   DS gives the distances SQRT((X-XO)^2+(Y-YO)^2) of all data points less
%   than the CUTOFF distance from the (m,n)th grid point, sorted in order
%   of increasing distance.  Entries farther than CUTOFF in all output
%   fields are filled with NaNs.
%
%   XS and YS are corresponding deviations X-XO and Y-YO from the grid 
%   point location to each data point.
%
%   The choice to put rows into cell arrays is made because for consistency
%   with SPHERESORT and for convenience in parallelizing POLYSMOOTH.
%   _________________________________________________________________
% 
%   Limiting output dimension
%
%   [DS,XS,YS]=TWODSORT(X,Y,XO,YO,[CUTOFF JMAX]), where the fifth input 
%   argument is a 2-vector, additionally specifies that number of rows of
%   in each cell of the output will be no larger than JMAX.  This option is
%   useful for the 'fixed population' algorithm in POLYSMOOTH.
%   _________________________________________________________________
% 
%   Additional input parameters
%
%   Let's say some additional variables Z1, Z2,...,ZK are given at the data
%   locations X,Y.  Then 
%   
%   [DS,XS,YS,Z1S,Z2S,...,ZKS]=
%
%                TWODSORT(X,Y,Z1,Z2,...,ZK,XO,YO,CUTOFF);
%
%   also returns the values of these variables.
%
%   Z1S, Z2S,...,ZKS are the same size as the other output arguments, and 
%   give the values of Z1, Z2,...,ZK sorted according to distance.
%
%   When there are multiple fields to be mapped, one may instead wish to 
%   use the approach described under "One grid, many fields" in POLYSMOOTH.
%   _________________________________________________________________
% 
%   See also SPHERESORT, POLYSMOOTH.
%
%   'twodsort --t' runs a test.
%
%   Usage: [ds,xs,ys,indexs]=twodsort(x,y,xo,yo,cutoff);
%          [ds,xs,ys,zs,indexs]=twodsort(x,y,zs,xo,yo,cutoff);
%          [ds,xs,ys,z1s,z2s,...,zNs]=twodsort(x,y,z1,z2,...,zN,xo,yo,cutoff);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2020 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    twodsort_test,return
end

xdata=varargin{1};
ydata=varargin{2};
xo=varargin{end-2};
yo=varargin{end-1};
cutoff=varargin{end};
varargin=varargin(3:end-3);

%In case max # points is not input
Ncutoff=inf;
if length(cutoff)==2
    Ncutoff=cutoff(2);
    cutoff=cutoff(1);
end

xo=xo(:);
yo=yo(:);

K=nargout;
if length(varargin)<K-3
    error('Not enough input arguments.')
end

%K,length(varargin)

if ~aresame(size(xdata),size(ydata))
    error('XDATA and YDATA must be the same size.')
end

for k=1:K
    varargout{k}=cell(length(yo),1);
end
indexo=find(isfinite(xdata)&isfinite(ydata));

if ~isempty(indexo)
    vcolon(xdata,ydata);    
    vindex(xdata,ydata,indexo,1);
else
    disp(['No finite data values.']), return
end

%indexall=1:length(xdata);

N=length(xo);
for i=1:length(yo)
    xp=vrep(xdata,N,2)-vrep(xo',length(xdata),1);
    yp=vrep(ydata-yo(i),N,2);
    d=sqrt(xp.^2+yp.^2);
    d(d>cutoff)=nan;
    xp(isnan(d))=nan;
    yp(isnan(d))=nan;
    
    %size(d)
    if ~allall(~isfinite(d))
        [dsort,sorter]=sort(d,'ascend');
        jj=vrep(1:N,length(xdata),1);    
        index=sub2ind(size(d),sorter,jj);
        
        xp=xp(index);
        yp=yp(index);
        
        L=min(find(sum(isfinite(dsort),2),1,'last'),Ncutoff);
        varargout{1}{i}=dsort(1:L,:);
        varargout{2}{i}=xp(1:L,:);
        varargout{3}{i}=yp(1:L,:);
        %vsize(d,xp,yp,dsort,sorter)
        for k=4:K
            temp=vrep(varargin{k-3}(indexo),N,2);
            temp=temp(index);
            temp(isnan(dsort))=nan;
            varargout{k}{i}=temp(1:L,:);
        end
    end
end
         
function[]=twodsort_test

[x,y,z]=peaks;
index=randperm(length(z(:)));
index=index(1:200);

[xdata,ydata,zdata]=vindex(x(:),y(:),z(:),index,1);

%Insert some NANs
xdata(1:7:end)=nan;

xo=(-3:.125:3);
yo=(-3:.125:3);
[xg,yg]=meshgrid(xo,yo);

%[ds,xs,ys,xs2,ys2]=twodsort(xdata,ydata,xdata,ydata,xo,yo,[1 20]);
%vsize(ds,xs,ys,xs2,ys2)

[ds,xs,ys,xs2,ys2]=twodsort(xdata,ydata,xdata,ydata,xo,yo,1);
for i=1:length(xs)
    bool(i)=aresame(xs{i}+vrep(xo,size(xs{i},1),1),xs2{i})&aresame(ys{i}+yo(i),ys2{i});
end
reporttest('TWODSORT consistency',allall(bool));
