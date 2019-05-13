function[x]=cellimag(x)
%CELLIMAG  Imaginary part of each element in a cell array.
%
%   XR=CELLIMAG(X) where X is a cell array of N arrays, is equivalent to 
%  
%      XR{1}=IMAG(X{1}), XR{2}=IMAG(X{2}), ..., XR{N}=IMAG(X{N})
%
%   thus returning the imaginary part of each element in the cell array.
%
%   See also JCELL.
% 
%   Usage: xi=cellimag(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2019 J.M. Lilly --- type 'help jlab_license' for details
 
if ~iscell(x)
    error('X must be a cell array.')
end

for i=1:length(x)
    x{i}=imag(x{i});
end
