function[varargout]=vrms(varargin)
%VRMS  Root-mean-square of non-NaN elements along a specified dimension.
%
%   Y=VRMS(X,DIM) takes the root-mean-square of all non-NaN elements of X 
%   along dimension DIM. 
%
%   That is, it is equivalent to Y=SQRT(VMEAN(SQUARED(X),DIM));
%
%   [Y,NUM]=VRMS(X,DIM) also outputs the number of non-NaN data points 
%   NUM, which has the same dimension as X.              
%
%   [Y1,Y2,...YN]=VRMS(X1,X2,...XN,DIM) or also works, where all the XN
%   are the same size.  
%
%   VRMS(X1,X2,...XN,DIM);  with no output arguments overwrites the 
%   original input variables.
%   __________________________________________________________________
%
%   Weighted means
%
%   Y=VRMS(X,DIM,W) forms the weighted root-mean-square of X along 
%   dimension DIM using weights W, an array  of the same size as X. 
%
%   This is equivalent to Y=SQRT(VMEAN(SQUARED(X),DIM,W));
%
%   In this case [Y,NUM]=VRMS(X,DIM,W) returns total weight used in 
%   forming the mean at each point, rather than the number of data points.
%
%   [Y1,Y2,...YN]=VRMS(X1,X2,...XN,DIM,W) also works. 
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2018 J.M. Lilly --- type 'help jlab_license' for details
 
if ~isscalar(varargin{end})
    w=varargin{end};
    varargin=varargin(1:end-1);
else
    w=[];
end
dim=varargin{end};

for i=1:length(varargin)-1
  if isempty(w)
      [varargout{i},numi{i}]=vsum(squared(varargin{i}),dim);
      numi{i}(numi{i}==0)=nan;
      varargout{i}=sqrt(varargout{i}./numi{i});
  else
      varargout{i}=vsum(squared(varargin{i}).*w,dim);
      numi{i}=vsum(w,dim);
      varargout{i}=sqrt(varargout{i}./numi{i});
  end
end

for i=length(varargin):nargout
  varargout{i}=numi{i-length(varargin)+1};
end

eval(to_overwrite(nargin-1))

