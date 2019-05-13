function[x]=cellreal(x)
%CELLREAL  Real part of each element in a cell array.
%
%   XR=CELLREAL(X) where X is a cell array of N arrays, is equivalent to 
%  
%      XR{1}=REAL(X{1}), XR{2}=REAL(X{2}), ..., XR{N}=REAL(X{N})
%
%   thus returning the real part of each element in the cell array.
%
%   See also JCELL.
%
%   Usage: xr=cellreal(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2019 J.M. Lilly --- type 'help jlab_license' for details
 
if ~iscell(x)
    error('X must be a cell array.')
end
for i=1:length(x)
    x{i}=real(x{i});
end