function[c]=materncfun(alpha)
%MATERNCFUN Returns the normalization function C_ALPHA for a Matern process.
%
%   MATERNCFUN is a low-level function used in the JMATERN package.
%
%   MATERNCFUN(ALPHA) returns 
%
%         (1/2/PI)*GAMMA(1/2)*GAMMA(ALPHA-1/2)/GAMMA(ALPHA) 
%
%   where GAMMA is the usual gamma function.
%
%   Usage: c=materncfun(alpha);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
 
%c=frac(1,2*pi)*frac(gamma(1/2)*gamma(alpha-1/2),gamma(alpha));
c=frac(1,2*sqrt(pi))*frac(gamma(alpha-1/2),gamma(alpha));



