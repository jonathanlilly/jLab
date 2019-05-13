function[n]=vsize(varargin)
%VSIZE  Returns the sizes of multiple arguments.
%
%   VSIZE(X1,X2,... XN) returns the sizes of the input variables.
%   These are organized in a matrix of size M x N where M is the number
%   of input variables and M is the maximum number of dimensions.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2018 J.M. Lilly --- type 'help jlab_license' for details    
  
if strcmpi(varargin,'--t')
   vsize_test,return
end
   
for i=1:length(varargin)
  x=size(varargin{i});
  xall{i}=x(:)';
  L(i)=length(x);
end

n=ones(length(varargin),max(L));

for i=1:length(varargin)
   index=1:length(xall{i});
   n(i,index)=xall{i};
end

function[]=vsize_test
x=ones(1,1);
y=ones(2,2);
z=ones(2,2,3);

x=vsize(x,y,z);
bool=allall(x==[1 1 1 ; 2 2 1 ;2 2 3]);
reporttest('VSIZE', bool)
