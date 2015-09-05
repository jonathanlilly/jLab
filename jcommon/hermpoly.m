function h=hermpoly(t,n)
% HERMPOLY  Hermite polynomials. [with F. Rekibi]
%
%   H=HERMPOLY(T,N) generates the first N+1 orthonormal Hermite polynomials 
%   [H0,...HN] on a time axis specfied by the column vector T.
%
%   In the nomenclature of Wikipedia,
%  
%      http://en.wikipedia.org/wiki/Hermite_polynomials#Definition 
%
%   the polynomials returned by HERMPOLY are the "physicists'" definition.
% 
%   Note that H(:,1) is the zeroth-order Hermite polynomial, which is equal
%   to one.  H(:,2) is the first-order polynomial, and so forth.  
%
%   See also HERMFUN, HERMEIG.
%
%   'hermpoly --t' runs a test.  
%
%   Usage:  h=hermpoly(t,n);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 F. Rekibi and J. M. Lilly 
%                         --- type 'help jlab_license' for details  
  
if strcmpi(t,'--t')
  hermpoly_test;return
end

if (size(t,1)==1)
	t=t(:);	% t vecteur colonne
end

h=zeros(length(t),n+1);

h(:,1)=ones(size(t,1),1);

if n+1>1
	h(:,2)=2.*t;
	if n+1>=2
		for i=3:n+1
			h(:,i)=2.*t.*h(:,i-1)-2.*(i-2).*h(:,i-2);
		end
	end
end
%h=h(:,2:end);
   
function[]=hermpoly_test
t=(-2:.01:2)';
n=4;
h=hermpoly(t,n);
h2(:,1)=1+0*t;
h2(:,2)=2.*t;     %H1
h2(:,3)=4*t.^2-2;     %H2
h2(:,4)=8.*t.^3-12.*t;   %H3
h2(:,5)=16.*t.^4-48.*t.^2+12;   %H4

tol=1e-10;
b=aresame(h,h2,tol);
reporttest('HERMPOLY Hermites 0-4 match analytic expressions',b);

