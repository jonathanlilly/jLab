function[z]=cellpair(x,y)
%CELLPAIR  Complex pairing for elements in a cell array.
%
%   Z=CELLPAIR(X,Y) where X and Y are cell arrays of N arrays, is
%   equivalent to 
%  
%      Z{1}=X{1}+1i*Y{1}, Z{1}=X{2}+1i*Y{2}, ..., Z{N}=X{N}+1i*Y{N}
%
%   so that is X and Y are both real-valued, the real part of each element
%   of Z is taken from X, and the imaginary part is taken from Y>
%
%   See also JCELL.
%
%   Usage: z=cellpair(x,y);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2019 J.M. Lilly --- type 'help jlab_license' for details
 
if ~iscell(x)||~iscell(y)
    error('X and Y must be cell arrays.')
end
for i=1:length(x)
    z{i,1}=x{i}+1i*y{i};
end