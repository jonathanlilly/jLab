function[b]=choose(n,k)
%CHOOSE  Binomial coefficient: CHOOSE(N,K) = N!K!/(N-K)!
%
%   CHOOSE(N,K) returns the binomial coefficient "N choose K", written
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
 

if length(k)==1&&length(n)~=1
    k=k+0*n;
elseif length(k)~=1&&length(n)==1
    n=n+0*k;
end
bool=k>n;
k(k>n)=0;

% b=zeros(size(n));
% n,k
% for i=1:length(n)
%     b(i)=(factorial(n(i)))./(factorial(k(i)).*factorial(n(i)-k(i)));
% end

b=(factorial(n))./(factorial(k).*factorial(n-k));
b(bool)=nan;