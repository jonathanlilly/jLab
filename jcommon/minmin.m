function[b]=minmin(x)
%MINMIN  MINMIN(X)=MIN(X(ISFINITE(X))) 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2012 J.M. Lilly --- type 'help jlab_license' for details        
b=min(x(isfinite(x(:))));
