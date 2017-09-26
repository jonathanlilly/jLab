function[z]=oprod(x,y)
%OPROD  Outer product:  OPROD(X,Y)=X*CONJ(Y') 
%
%   OPROD(X,Y) calculates the outer product of two column vectors X and Y. 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2003--2015 J.M. Lilly --- type 'help jlab_license' for details        
  
if ~(iscolumn(x) || isscalar(x)) || ~(iscolumn(y) || isscalar(y))
  error('X and Y must both be column vectors or scalars.')
else
  z=x*conj(y');  
end

