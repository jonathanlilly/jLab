function[varargout]=cellindex(varargin)
%CELLINDEX  Applies a cell array of indices to a cell array of column vectors.
%
%   Y=CELLINDEX(X,INDEX) where X is a cell array of N column vectors, and
%   INDEX is a cell array of indices into those vectors
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
%       INDEX{1}=I1, INDEX{2}=I2,..., INDEX{N}=IN
% 
%   applies INDEX to each element of X, returning Y defined by
%
%      Y{1}=X1(I1), Y{2}=X2(I2),..., Y{N}=XN(IN).
%
%   The output arrays YN are the same size as the indices IN. 
%   __________________________________________________________________
%
%   Multiple input arguments
%
%   [Y1,Y2,...,YM]=CELLINDEX(X1,X2,...,XM,INDEX) also works, with M 
%   different cell arrays XM of identical size input, leading to M cell
%   arrays YM all having the same size as INDEX.
%   __________________________________________________________________
%
%   See also JCELL, VINDEX.
%
%   Usage: y=cellindex(x,index);
%          [y1,y2,y3]=cellindex(x1,x2,x3,index);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2014 J.M. Lilly --- type 'help jlab_license' for details
 

index=varargin{end};

for j=1:length(index)  %Loop over number of cells
    for i=1:length(varargin)-1  %Loop over variables
        varargout{i}{j,1}=varargin{i}{j}(index{j}); 
    end
end

