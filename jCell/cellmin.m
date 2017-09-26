function[minx]=cellmin(x)
%CELLMIN  Minimum of each element in a cell array.
%
%   M=CELLMIN(X) where X is a cell array of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   returns the N x 1 array of minimum values M with
%  
%      M(1)=MIN(X1(:)),  M(2)=MIN(X2(:)),...,  M(N)=MIN(XN(:)).
%
%   M is the same size as X.
%
%   CELLMIN requires that X have four or fewer dimensions.
%
%   See also JCELL.
%
%   Usage: m=cellmin(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009 J.M. Lilly --- type 'help jlab_license' for details

if ~iscell(x)
    error('X must be a cell array.')
end
 
minx=nan*zeros(size(x));
for i=1:size(x,1)
    for j=1:size(x,2)
        for k=1:size(x,3)
            for l=1:size(x,4)
                m=minmin(x{i,j,k,l});
                if ~isempty(m)
                       minx(i,j,k,l)=m;
                end
            end
        end
    end
end