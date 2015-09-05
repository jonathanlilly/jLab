function[varargout]=lnsd(varargin)
% LNSD  Last non-singleton dimension of an array.
%  
%   LNSD(X) returns the number of the last non-singleton dimension of X.
%
%   This provides a useful definition of the dimensionality of X.  Unlike 
%   Matlab's NDIMS, which thinks that a column vector and a scalar both 
%   have dimension 2, LNSD defines the dimension of a scalar to be zero and
%   that of a column vector to be one, while a row vector has an LNSD of 2.
%  
%   Note that LNSD(X) changes if X is permuted. 
%  
%   [N1,N2,...NM]=LNSD(X1,X2,...XM) returns the dimensions of multiple 
%   input arguments.  If zero or one output arguments are given, a single 
%   row array [N1 N2 ...NM] is output.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2015 J.M. Lilly --- type 'help jlab_license' for details        
  
if strcmpi(varargin{1}, '--t')
  lnsd_test,return
end


for i=1:nargin
  x=varargin{i};
  sx=size(x);
  if isempty(x)
    nd1=nan;
  elseif isscalar(x)
    nd1=0;
  else
   nd1=find(sx~=1,1,'last');
  end

  nd(i)=nd1;
  varargout{i}=nd(i);
end

if nargout==0 || nargout==1
   varargout{1}=nd;
end

function[]=lnsd_test

x=lnsd([],1,(1:10),[ (1:10)' (1:10)']);
bool(1)=aresame(x,[nan 0 2 2]);
%disp('Should be NAN 0 2 2')

x=(1:10)';
z(:,:,3)=x;
q(:,:,:,1)=z;
q(:,:,:,2)=z;
c(1,1,:)=x;  

x=lnsd(x,z,q,c,permute(c,[3 2 1]));
bool(2)=aresame(x,[1 3 4 3 1]);
%disp('Should be 1 3 4 3 1')
reporttest('LNSD',all(bool))

