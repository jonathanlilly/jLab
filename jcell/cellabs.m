function[x]=cellabs(x)
%CELLABS  Absolute value of each element in a cell array.
%
%   XA=CELLABS(X) where X is a cell array of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   returns the length N cell array of XA with 
%  
%      XA(1)=ABS(X1), XA(2)=ABS(X2),..., XA(N)=ABS(XN).
%
%   See also JCELL.
% 
%   Usage: xa=cellabs(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009 J.M. Lilly --- type 'help jlab_license' for details
 
if ~iscell(x)
    error('X must be a cell array.')
end
 
for i=1:length(x)
    x{i}=abs(x{i});
end
 