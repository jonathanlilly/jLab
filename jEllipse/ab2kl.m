function[K,L]=ab2kl(A,B)
%AB2KL  Converts A and B to ellipse parameters Kappa and Lambda.
%
%   [K,L]=AB2KL(A,B)  Convert ellipse semi-major and semi-minor
%   axes A and B to ellipse amplitude K and ellipse parameter L.
%
%   K and L are related to A and B by
%
%          K   = SQRT( (A^2 +B^2) / 2 )
%          L   = SIGN(B) * (A^2 - B^2 ) / (A^2 + B^2)
%
%   as discussed in Lilly and Gascard (2006).  
% 
%   A and B are either arrays of the same size, or one may be an array
%   and one may be a scalar.
%
%   Usage: [kappa,lambda]=ab2kl(a,b)
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2019 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(A, '--t')
    ab2kl_test,return
end

if numel(A)==1&&numel(B)~=1
    A=A+zeros(size(B));
elseif numel(B)==1&&numel(A)~=1
    B=B+zeros(size(A));
end

L=sign(B).*frac(A.^2-B.^2,A.^2+B.^2);
K=sqrt(frac(A.^2+B.^2,2));
L(B==0)=1;

function[]=ab2kl_test

x=rand(100,1)*2-1;

[a,b]=kl2ab(1,x);
[k,l]=ab2kl(a,b);
l2=ecconv(b./a,'ell2lin');

tol=1e-6;
reporttest('AB2KL matches ECCONV',aresame(l,l2,tol))
