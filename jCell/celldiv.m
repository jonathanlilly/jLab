function[z]=celldiv(x,y)
%CELLDIV  Division acting on each element in a cell array.
%
%   Z=CELLDIV(X,Y) where X and Y are both cell arrays of N arrays, with the
%   corresponding elements in X and in Y having the same size, returns the
%   cell array Z containing their ratios, 
%
%       Z{1}=Y{1}./X{1}, Z{2}=Y{2}./X{2},..., Z{N}=Y{N}./X{N}.
%
%   One of X or Y may also be a scalar or a numeric array of the same 
%   length as the other input argument.
%
%   Usage: z=celldiv(x,y);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015--2019 J.M. Lilly --- type 'help jlab_license' for details

if ~iscell(x)&&iscell(y)
    z=y;
    if length(x)==1
        x=x.*ones(length(y),1);
    end
    if length(x)==length(y)
        for i=1:length(x)
            z{i}=y{i}./x(i);
        end
    else 
        error('LENGTH(X) and LENGTH(Y) must be the same.')
    end
elseif iscell(x)&&~iscell(y)
    z=x;
    if length(y)==1
        y=y.*ones(length(x),1);
    end
    if length(x)==length(y)
        for i=1:length(x)
            z{i}=y(i)./x{i};
        end
    else 
        error('LENGTH(X) and LENGTH(Y) must be the same.')
    end
elseif iscell(x)&&iscell(y)
    z=y;
    for i=1:length(x)
        z{i}=y{i}./x{i};
    end
end

 