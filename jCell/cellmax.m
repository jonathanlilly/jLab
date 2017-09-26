function[maxx]=cellmax(x)
%CELLMAX  Maximum of each element in a cell array.
%
%   M=CELLMAX(X) where X is a cell array of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   returns the N x 1 array of maximum values M with
%  
%      M(1)=MAX(X1(:)),  M(2)=MAX(X2(:)),...,  M(N)=MAX(XN(:)).
%
%   M is the same size as X.
%
%   CELLMAX requires that X have four or fewer dimensions.
%
%   See also JCELL.
%
%   Usage: m=cellmax(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009 J.M. Lilly --- type 'help jlab_license' for details

if ~iscell(x)
    error('X must be a cell array.')
end
 
maxx=nan*zeros(size(x));
for i=1:size(x,1)
    for j=1:size(x,2)
        for k=1:size(x,3)
            for l=1:size(x,4)
                m=maxmax(x{i,j,k,l});
                if ~isempty(m)
                       maxx(i,j,k,l)=m;
                end
            end
        end
    end
end