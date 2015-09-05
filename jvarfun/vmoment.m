function[varargout] = vmoment(varargin)
%VMOMENT Central moment over finite elements along a specfied dimension.
%
%   Y=VMOMENT(X,N,DIM) finds the Nth central moment of all finite elements 
%   of X along dimension DIM. 
%                                                                         
%   [Y,NUM]=VMOMENT(X,N,DIM) also outputs the number of good data points NUM, 
%   which has the same dimension as X.                                
%
%   [Y1,Y2,...YN]=VMOMENT(X1,X2,...XN,N,DIM) also works.
%
%   VMOMENT(X1,X2,...XN,N,DIM);  with no output arguments overwrites the 
%   original input variables.      
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2012 J.M. Lilly --- type 'help jlab_license' for details    
    
if strcmpi(varargin{1}, '--t')
  vmoment_test,return
end

ndim=varargin{end};
n=varargin{end-1};

for i=1:length(varargin)-2
  x=varargin{i};
  m=vmean(x,ndim);
  m=vrep(m,size(x,ndim),ndim);
  
  [varargout{i},numi{i}]=vmean(abs(x-m).^n,ndim);
end

for i=length(varargin)-1:nargout
  varargout{i}=numi{i-length(varargin)+2};
end

eval(to_overwrite(nargin-2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=vmoment_test
x1=[1 2 3 nan];
x2=x1;
ans1=2/3;

vmoment(x1,x2,2,2);
reporttest('VMOMENT output overwrite', aresame(x1,ans1) && aresame(x2,ans1))

x1=[1 2 3 nan];
ans2=3;
[y1,y2]=vmoment(x1,2,2);
reporttest('VMOMENT moment & num', aresame(y1,ans1) && aresame(y2,ans2))



