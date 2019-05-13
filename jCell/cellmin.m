function[minx]=cellmin(x)
%CELLMIN  Minimum of each element in a cell array.
%
%   M=CELLMIN(X) where X is a cell array of N arrays, is equivalent to
%  
%      M(1,1)=MIN(X{1}(:)),  M(2,1)=MIN(X{2}(:)), ...,  M(N,1)=MIN(X{N}(:))
%
%   thus returning an N x 1 array containing the minimum values of each
%   element in the cell array.
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