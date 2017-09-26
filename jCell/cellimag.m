function[x]=cellimag(x)
%CELLIMAG  Imaginary part of each element in a cell array.
%
%   XI=CELLIMAG(X) where X is a cell array of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   returns the length N cell array XI with 
%  
%      XI(1)=IMAG(X1), XI(2)=IMAG(X2),..., XI(N)=IMAG(XN).
%
%   See also JCELL.
% 
%   Usage: xi=cellimag(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009 J.M. Lilly --- type 'help jlab_license' for details
 
if ~iscell(x)
    error('X must be a cell array.')
end

for i=1:length(x)
    x{i}=imag(x{i});
end
