function[b]=choose(n,k)
%CHOOSE  Binomial coefficient: CHOOSE(N,K) = N!K!/(N-K)!
%
%   CHOOSE(N,M) returns the binomial coefficient "N choose K", written
%
%         ( N )         N!    
%         (   )   =  --------
%         ( K )      K!(N-K)!
%
%   with a little imagination for the left-hand side notation.  
%
%   N and K can each either be an array or a scalar, and must of course
%   be the same size if they are both arrays. 
%
%   Usage: b=choose(n,k);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2015 J.M. Lilly --- type 'help jlab_license' for details
 
b=(factorial(n))./(factorial(k).*factorial(n-k));
