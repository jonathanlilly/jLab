function[b]=maxmax(x)
%MAXMAX  MAXMAX(X)=MAX(X(~ISNAN(X(:))))
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2012 J.M. Lilly --- type 'help jlab_license' for details        
b=max(x(~isnan(x(:))));