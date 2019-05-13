function[x]=celllog10(x)
%CELLLOG10 Base ten logarithm of each element in a cell array.
%
%   X=CELLLOG10(X) where X is a cell array of N arrays, is equivalent to
%  
%      M(1,1)=LOG10(X{1}(:)), ...,  M(N,1)=LOG10(X{N}(:))
%
%   thus returning an N x 1 array containing the based 10 logarithm of all
%   values in each element of the cell array.
%
%   See also JCELL.
% 
%   Usage: xa=celllog10(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2019 J.M. Lilly --- type 'help jlab_license' for details
 
if ~iscell(x)
    error('X must be a cell array.')
end
 
for i=1:length(x)
    x{i}=log10(x{i});
end
 