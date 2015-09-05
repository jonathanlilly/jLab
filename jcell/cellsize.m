function[s]=cellsize(x,n)
%CELLSIZE  Size of each element in a cell array along specified dimension.
%
%   S=CELLSIZE(X,DIM) where X is a cell array of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   returns the N x 1 array of lengths S with 
%  
%      S(1)=SIZE(X1,DIM), S(2)=SIZE(X2,DIM),..., S(N)=SIZE(XN,DIM).
%
%   S is the same size as X. 
%
%   CELLSIZE requires the cell array X to have four or fewer dimensions.
%
%   See also CELLENGTH.
%
%   Usage: s=cellsize(x,dim);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(x, '--t')
     cellsize_test,return
end
 
s=zeros(size(x));
for i=1:size(x,1)
    for j=1:size(x,2)
        for k=1:size(x,3)
            for l=1:size(x,4)
                s(i,j,k,l)=size(x{i,j,k,l},n);
            end
        end
    end
end

function[]=cellsize_test
 
x{1}=[1 2; 2 2];
x{2}=[];
x{3}=[2 2];

reporttest('CELLSIZE',aresame(cellsize(x,1),[2 0 1])&&aresame(cellsize(x,2),[2 0 2]))
