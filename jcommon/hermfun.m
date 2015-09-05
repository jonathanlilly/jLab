function h=hermfun(t,j)
%HERMFUN  Orthonormal Hermite functions. [with F. Rekibi]
%
%   H=HERMFUN(T,N) generates the fisrt N+1 orthonormal Hermite functions
%   [H0,...HN] on a time axis specfied by the column vector T.
%
%   HERMFUN uses the expression of Simons et al. 2003.
%  
%   Note that H(:,1) is the zeroth-order Hermite function, which is equal
%   to a Gaussian.  H(:,2) is the first-order function, and so forth.  
%  
%   See also HERMPOLY.
%
%   'hermfun --f' generates a sample figure; compare with the Hermite 
%   function figure at
%
%     http://en.wikipedia.org/wiki/Hermite_polynomials#Definition
%
%   Usage:  h=hermfun(t,n);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2004--2015 F. Rekibi and J. M. Lilly
%                         --- type 'help jlab_license' for details

% 05.08.07  JML fixed bug to include N+1 columns

%   'hermfun --f' generates a sample figure; compare with Figure 2 
%   of Simons, van der Hilst, and Zuber (2003), JGR 108 B5.

if strcmpi(t,'--f')
  type makefigs_hermfun
  makefigs_hermfun;
  return
end

if size(t,1)==1
    t=t';
end

H=hermpoly(t,j);

E=exp(-t.^2/2)*ones(1,j+1);
HE=H.*E;

h=zeros(length(t),j);
for k=1:j+1
	h(:,k)=HE(:,k)*frac(1,pi^(1/4))*frac(1,sqrt(2^(k-1)*factorial(k-1)));
end

