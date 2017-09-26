function[varargout]=vtranspose(varargin)
%VTRANSPOSE  Transpose multiple input arguments simultaneously.
%
%   [Y1,Y2, ... YN]=VTRANSPOSE(X1,X2, ... XN) is equivalent to
%
%       Y1=CONJ(X1'); Y2=CONJ(X2'); ... YN=CONJ(XN');
%
%   for matrix arguments X1,X2, ... XN.  This is the matrix transpose.
%   Recall X1' is the conjugate transpose of X1 in Matlab's notation. 
%
%   VTRANSPOSE(X1,X2,...XN);  with no output arguments overwrites the 
%   original input variables.    
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014 J.M. Lilly --- type 'help jlab_license' for details

for i=1:length(varargin)
  %varargout{i}=permute(varargin{i},[2 1]);
  varargout{i}=conj(varargin{i}');
end

eval(to_overwrite(nargin))


