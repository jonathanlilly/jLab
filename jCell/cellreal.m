function[x]=cellreal(x)
%CELLREAL  Real part of each element in a cell array.
%
%   XR=CELLREAL(X) where X is a cell array of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   returns the length N cell array XR with 
%  
%      XR(1)=REAL(X1), XR(2)=REAL(X2),..., XR(N)=REAL(XN).
%
%   See also JCELL.
%
%   Usage: xr=cellreal(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009 J.M. Lilly --- type 'help jlab_license' for details
 
if ~iscell(x)
    error('X must be a cell array.')
end
for i=1:length(x)
    x{i}=real(x{i});
end