function[varargout]=vrep(varargin)
%VREP  Replicates an array along a specified dimension.
%
%   Y=VREP(X,N,DIM) replicates the array by N times along dimension
%   dimension DIM.  For instance:   
%                                                                         
%        VREP([1:4]',3,2)=[ [1:4]' [1:4]' [1:4]' ]                            
%                                                                         
%   This is often useful in array algebra.            
%
%   IF N and DIM are arrays of length M, then X is replicated along each of
%   the M different dimensions: N(1) times along dimensions DIM(1), etc.
%
%   [Y1,Y2,...,YP]=VREP(X1,X2,...,XP,N,DIM) also works.
%                                                                         
%   See also VINDEX, DIM.      
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2014 J.M. Lilly --- type 'help jlab_license' for details    

if strcmpi(varargin{1}, '--t')
  vrep_test,return
end

n=varargin{end-1};
ndim=varargin{end};

for i=1:length(varargin)-2
   varargout{i}=varargin{i};
end

for i=1:length(varargin)-2
    for j=1:length(n)
        varargout{i}=vrep1(varargout{i},n(j),ndim(j));
    end
end


eval(to_overwrite(nargin-2))
 
%You would think Matlab would provide a simpler way to do this.
function[y]=vrep1(x,n,dim)
  
str='y=repmat(x,[';
ndx=ndims(x);
for i=1:max(ndx,dim)
    if i~=dim
        str=[str '1,'];
    else
	str=[str 'n,'];
    end
end
str=[str(1:end-1) ']);'];
eval(str);


function[]=vrep_test

ans1=vrep((1:4)',3,2);
ans2=[ (1:4)' (1:4)' (1:4)' ];
reporttest('VREP', aresame(ans1,ans2))

x1=(1:4)';x2=(1:4)';
vrep(x1,x2,3,2);
reporttest('VREP output redirect', aresame(x1,ans1) && aresame(x2,ans2))
