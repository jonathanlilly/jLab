function[z]=celladd(x,y)
%CELLADD  Addition acting on each element in a cell array.
%
%   Z=CELLADD(X,Y) where X and Y are both cell arrays of N arrays, with 
%   corresponding elements in X and in Y having the same size, returns the
%   cell array Z containing their sums, 
%
%       Z{1}=X{1}+Y{1}, Z{2}=X{2}+Y{2},..., Z{N}=X{N}+Y{N}.
%
%   One of X or Y may also be a scalar or a numeric array of the same 
%   length as the other input argument.
%
%   Usage: z=celladd(x,y);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2019 J.M. Lilly --- type 'help jlab_license' for details
 
if ~iscell(x)&&iscell(y)
    if length(x)==1
        x=x.*ones(length(y),1);
    else
        x=x(:);
    end
elseif iscell(x)&&~iscell(y)
    if length(y)==1
        y=y.*ones(length(x),1);
    else
        y=y(:);
    end
end

if iscell(x)
    z=x;
else
    z=y;
end

for i=1:length(y)
     if iscell(x)&&iscell(y)
         z{i}=x{i}+y{i};
     elseif ~iscell(x)&&iscell(y)
         z{i}=x(i)+y{i};
     elseif iscell(x)&&~iscell(y)
         z{i}=x{i}+y(i);
     end
end
    
