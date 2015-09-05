function[y]=cellmult(x,y)
%CELLMULT  Multiplication acting on each element in a cell array.
%
%   Z=CELLMULT(X,Y) where X and Y are both cell arrays of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
%  
%   with the Nth element in X and in Y having the same size, returns the
%   cell array Z containing their products, 
%
%       Z{1}=X{1}.*Y{1}, Z{2}=X{2}.*Y{2},..., Z{N}=X{N}.*Y{N}.
%
%   X may also be a scalar or an array of N=LENGTH(Y) numbers. 
%
%   Usage: z=cellmult(x,y);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2015 J.M. Lilly --- type 'help jlab_license' for details

if ~iscell(x)
    if length(x)==1
        x=x.*ones(length(y),1);
    end
    if length(x)==length(y)
        for i=1:length(x)
            y{i}=x(i).*y{i};
        end
    else 
        error('LENGTH(X) and LENGTH(Y) must be the same.')
    end
else
    for i=1:length(x)
        y{i}=x{i}.*y{i};
    end
end

    
function[]=cellmult_test
 
