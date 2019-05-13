function[x]=cellabs(x)
%CELLABS  Absolute value of each element in a cell array.
%
%   XA=CELLABS(X) where X is a cell array of N arrays, is equivalent to
%  
%      XA{1}=ABS(X{1}), XA{2}=ABS(X{2}), ..., XA{N}=ABS(X{N})
%
%   thus returning the  absolute value of each element in the cell array.
%
%   See also JCELL.
% 
%   Usage: xa=cellabs(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2019 J.M. Lilly --- type 'help jlab_license' for details
 
if ~iscell(x)
    error('X must be a cell array.')
end
 
for i=1:length(x)
    x{i}=abs(x{i});
end
 