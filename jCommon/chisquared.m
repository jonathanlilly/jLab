function[fx]=chisquared(x,k)
%CHISQUARED  The chi-squared distribution.
%
%   FX=CHISQUARED(X,K) returns the chi-squared probability distribution
%   with K degrees of freedom at values X, which may be an array.
%
%   Usage: fx=chisquared(x,k);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2020 J.M. Lilly --- type 'help jlab_license' for details
 
fx=frac(x.^(k/2-1).*exp(-x/2),2.^(k/2).*gamma(k/2));
fx(x<=0)=0;

