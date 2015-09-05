function[x]=celllog10(x)
%CELLLOG10 Base ten logarithm of each element in a cell array.
%
%   XA=CELLLOG10(X) where X is a cell array of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   returns the length N cell array of XA with 
%  
%      XA(1)=LOG10(X1), XA(2)=LOG10(X2),..., XA(N)=LOG10(XN).
%
%   See also JCELL.
% 
%   Usage: xa=celllog10(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014 J.M. Lilly --- type 'help jlab_license' for details
 
if ~iscell(x)
    error('X must be a cell array.')
end
 
for i=1:length(x)
    x{i}=log10(x{i});
end
 