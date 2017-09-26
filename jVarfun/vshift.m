function[varargout]=vshift(varargin)
% VSHIFT  Cycles the elements of an array along a specified dimension.
%
%   Y=VSHIFT(X,N,DIM) cycles the elements of X N places along dimension DIM.
%  
%   Example: x=[1 2 3 4 5];
%            vshift(x,+1,2)=[2 3 4 5 1]           
%            vshift(x,-1,2)=[5 1 2 3 4]           
%
%   Note shifting by N and then by -N recovers the original array. 
%
%   [Y1,Y2,...YN]=VSHIFT(X1,X2,...XN,N,DIM) also works.
%
%   VSHIFT(X1,X2,...XN,N,DIM); with no arguments overwrite the original 
%   input variables.
%
%   Note that VSHIFT is similar to Matlab's CIRCSHIFT, but has the opposite
%   convection for positive and negative shifts. 
%
%   ------------------------------------------------------------------
%   Y=VSHIFT(X,N,DIM,INDEX,DIM2) applies this shift selectively, only to
%   that subset of X obtained by indexing X with INDEX along DIM2, i.e.
%
%		    1 2      DIM2     DIMS(X)
%		    | |        |         |
%		  X(:,:, ... INDEX, ..., :)	
%
%   is cycled N places along dimension DIM, but the remainder of X is not. 
%   DIM and DIM2 cannot be the same.  The above extensions to multiple 
%   output varibles work in this case as well.  
%   ------------------------------------------------------------------
%
%   See also: VINDEX, CIRCSHIFT.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2016 J.M. Lilly --- type 'help jlab_license' for details    
  

if strcmpi(varargin{1}, '--t')
  vshift_test,return
end

%/********************************************************
%Sort out input arguments
nax=2;
if nargin>4
  if  length(varargin{end-3}(:))==1
     nax=4;
     dim2=(varargin{end});
     jj=varargin{end-1}(:);
     dim=varargin{end-2};
     n=varargin{end-3};
     if dim==dim2
       error('DIM and DIM2 cannot be the same.')
     end
     %n, jj, dim,dim2
  end
end

if nax==2
  dim=varargin{end};
  n=varargin{end-1};
end
%\********************************************************

for i=1:length(varargin)-nax
  if nax==2
    varargout{i}=vshift1(varargin{i},n,dim);
  else
    varargout{i}=vshift1(varargin{i},n,dim,jj,dim2);
  end  
end
eval(to_overwrite(nargin-nax))


function[y]=vshift1(x,n,ndim,jj,ndim2)


N=size(x,ndim);
if n>0
    ii=[(n+1:N) (1:n)];
elseif n<0
    n=-n;
    ii=[(N-(n-1):N) (1:N-n)];
elseif n==0
    ii=(1:N);
end

if nargin==3
%   Same as using circshift... these are the same speed
    %array=zeros(max(ndim,numel(size(x))),1);
    %array(ndim)=-n;
    %ndim
 %   y=circshift(x,-n,ndim);
    y=vindex(x,ii,ndim);
else
    y1=vindex(x,jj,ndim2);
    y1=vindex(y1,ii,ndim);
    %vsize(x,y1,jj,ndim2);
    %x,y1,ii,jj,ndim2
    y=vindexinto(x,y1,jj,ndim2);
end

function[]=vshift_test
x=(1:10);
ans1=[(2:10) 1];
reporttest('VSHIFT col case', aresame(vshift(x,1,2),ans1))

x=(1:10)';
ans1=[10 (1:9)]';
reporttest('VSHIFT row case', aresame(vshift(x,-1,1),ans1))

clear x ans1
x(:,:,1)=[1 2; 3 4];
x(:,:,2)=2*[1 2; 3 4];
ans1(:,:,2)=x(:,:,1);
ans1(:,:,1)=x(:,:,2);
reporttest('VSHIFT mat case', aresame(vshift(x,1,3),ans1))

clear x ans1
x(:,:,1)=[1 2; 3 4];
x(:,:,2)=2*[1 2; 3 4];
ans1(:,:,1)=[3 2;1 4];
ans1(:,:,2)=2*[3 2;1 4];
reporttest('VSHIFT mat selective case one', aresame(vshift(x,1,1,1,2),ans1))

clear x ans1
x(:,:,1)=[1 2; 3 4];
x(:,:,2)=2*[1 2; 3 4];
ans1(:,:,1)=[3 4;1 2];
ans1(:,:,2)=2*[3 4;1 2];
reporttest('VSHIFT mat selective case two', aresame(vshift(x,1,1,1:2,2),ans1))

