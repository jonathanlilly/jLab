function[varargout]=twodsort(varargin)
%TWODSORT  Distances from data points to nearby grid points.
%
%   [D,XD,YD]=TWODSORT(X,Y,XO,YO,CUTOFF) returns distances D between data
%   points at locations X,Y and nearby grid points at XO,YO, computed
%   efficiently and organized in a convenient manner.
%
%   XD and YD are deviations from each data point location to the grid 
%   point location. CUTOFF is the maximum distance for which we wish the
%   output arrays to be returned.
% 
%   X and Y are arrays of the same size into data point locations. XO and
%   YO  are arrays of length M and N, say, specifying the bin center
%   locations of an M x N matrix of grid points, i.e.
%
%       XO =  [XO_1;    YO= [YO_1 YO_2 ... YO_N]. 
%              XO_2; 
%               ...
%              XO_M]
%
%   The output arrays are then each M x N x P arrays of column vectors, 
%   where P is the maximum number of points in the CUTOFF neighborhood at
%   any grid point.  Points farther away are filled with NANs.
%
%   D gives the distances of those data points less than CUTOFF from the
%   (m,n)th grid point, sorted in order of increasing distance.  XD and YD 
%   are deviations of each data point from the (m,n)th grid point location.
%   _________________________________________________________________
% 
%   Additional input parameters
%
%   Let's say some additional variables Z1, Z2,...,ZK are given at the data
%   locations X,Y.  Then 
%   
%   [D,XD,YD,Z1D,Z2D,...,ZKD]=
%
%                TWODSORT(X,Y,Z1,Z2,...,ZK,XO,YO,CUTOFF);
%
%   also returns the values of these variables.
%
%   Z1D, Z2D,...,ZKD are the same size as the other output arguments.  Z1D
%   then gives the value of Z1 at data points no more than CUTOFF 
%   kilometers from the (m,n)th grid point, etc.
%
%   To output an index into the data point locations, use
%  
%   [D,X,Y,INDEX]= TWODSORT(X,Y,[1:LENGTH(X(:))],XO,YO,CUTOFF).
%   _________________________________________________________________
% 
%   See also SPHERESORT, POLYSMOOTH.
%
%   'twodsort --t' runs a test.
%
%   Usage: [ds,xs,ys]=twodsort(x,y,xo,yo,cutoff);
%          [ds,xs,ys,zs]=twodsort(x,y,zs,xo,yo,cutoff);
%          [xs,ys,z1s,z2s,...,zNs]=twodsort(x,y,z1,z2,...,zN,xo,yo,cutoff);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    twodsort_test,return
end

xdata=varargin{1};
ydata=varargin{2};
xo=varargin{end-2};
yo=varargin{end-1};
cutoff=varargin{end};

[xo,yo]=meshgrid(xo,yo);

K=nargout;
if length(varargin)<K+2
    error('Not enough input arguments.')
end

if ~aresame(size(xdata),size(ydata))
    error('XDATA and YDATA must be the same size.')
end

for k=1:K
    varargout{k}=cell(size(xo,1),size(xo,2));
end
indexo=find(isfinite(xdata)&isfinite(ydata));

if ~isempty(indexo)
    vcolon(xdata,ydata);    
    vindex(xdata,ydata,indexo,1);
else
    disp(['No finite data values.']), return
end

for j=1:length(xo(:))
     xp=xdata-xo(j);
     yp=ydata-yo(j);
     d=sqrt(xp.^2+yp.^2);
     index=find(d<cutoff);
     if ~isempty(index);
         [dsort,sorter]=sort(d(index),'ascend');
         index=index(sorter);
         %aresame(d(index),dsort)
         varargout{1}{j}=dsort;
         varargout{2}{j}=xp(index);
         varargout{3}{j}=yp(index);
         for k=4:K
             varargout{k}{j}=varargin{k-1}(indexo(index));
         end
     end
end

cells=varargout;
xs=cells{1};
L=zeros(size(xs));
for i=1:size(xs,1)
    for j=1:size(xs,2)
        L(i,j)=length(xs{i,j});
    end
end
maxL=maxmax(L);

for k=1:K
    varargout{k}=vzeros(size(xs,1),size(xs,2),maxL,'nan');
end

for i=1:size(xs,1)
    for j=1:size(xs,2)
        if L(i,j)>=1
            for k=1:K
                varargout{k}(i,j,1:L(i,j))=cells{k}{i,j};
            end
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

[ds,xs,ys,xs2,ys2]=twodsort(xdata,ydata,xdata,ydata,xo,yo,1);
reporttest('TWODSORT consistency',aresame(xs+vrep(xg,size(xs,3),3),xs2)&&aresame(ys+vrep(yg,size(xs,3),3),ys2))

ds2=nan*zeros(size(ds));
for i=1:size(ds,1)
    for j=1:size(ds,2)
        for k=1:size(ds,3)
            ds2(i,j,k)=sqrt((xs(i,j,k)).^2+(ys(i,j,k)).^2);
        end
    end
end
reporttest('TWODSORT distance',aresame(ds,ds2,1e-10))
 
