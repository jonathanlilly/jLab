function[varargout]=cellmean(varargin)
%CELLMEAN  Mean value of each element a cell array or set of cell arrays.
%
%   M=CELLMEAN(X) where X is a cell array of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   returns the N x 1 array of mean values M with
%  
%      M(1)=MEAN(X1(:)),  M(2)=MEAN(X2(:)),...,  M(N)=MEAN(XN(:)).
%
%   In taking the mean, non-finite values are ignored, as in VMEAN.
%   M is a column vector of the same length as X.  
%
%   [M1,M2,...,MP]=CELLMEAN(X1,X2,...,XP) also works for P different 
%   input arguments. 
%
%   CELLMEAN(X1,X2,...XP);  with no output arguments overwrites the 
%   original input variables.
%
%   See also JCELL, VTOOLS.
%
%   Usage: m=cellmean(x);
%          [m1,m2,m3]=cellmean(x1,x2,x3);
%          cellmean(x1,x2,x3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014 J.M. Lilly --- type 'help jlab_license' for details

if ~iscell(varargin{1})
    error('X must be a cell array.')
end
 
for i=1:nargin
    varargout{i}=cellmean_one(varargin{i});
end

eval(to_overwrite(nargin))


function[y]=cellmean_one(x)
y=vmean(col2mat(cell2col(x)),1);

if length(y)~=length(x)
   disp('Fast version of CELLMEAN failed, perhaps because X contained NANs.')
   disp('Be patient, implementing slow version.')
   y=nan*zeros(length(x),1);
   for i=1:length(x)
      y(i)=vmean(vcolon(x{i}),1);
   end
end

