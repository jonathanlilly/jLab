function[s]=cellength(x)
%CELLENGTH  Length of each element in a cell array.
%
%   L=CELLENGTH(X) where X is a cell array of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   returns the N x 1 array of lengths L with 
%  
%      L(1)=LENGTH(X1), L(2)=LENGTH(X2),..., L(N)=LENGTH(XN).
%
%   L is the same size as X. 
%
%   CELLENGTH requires the cell array X to have four or fewer dimensions.
%
%   See also CELLSIZE.
%
%   Usage: len=cellength(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(x, '--t')
     cellength_test,return
end

if ~iscell(x)
    error('X must be a cell array.')
end

s=zeros(size(x));
for i=1:size(x,1)
    for j=1:size(x,2)
        for k=1:size(x,3)
            for l=1:size(x,4)
                s(i,j,k,l)=length(x{i,j,k,l});
            end
        end
    end
end

function[]=cellength_test
 
x{1}=[1 2; 2 2];
x{2}=[];
x{3}=[2 2];

reporttest('CELLENGTH',aresame(cellsize(x,2),[2 0 2]))
