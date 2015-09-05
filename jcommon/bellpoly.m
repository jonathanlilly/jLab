function[varargout]=bellpoly(varargin)
%BELLPOLY  Complete Bell polynomials.
%
%   BN=BELLPOLY(K1,K2,...KN) with N arguments returns the Nth order 
%   complete Bell polynomial BN.  
%
%   For details, see the article at Wikipedia:
%  
%        http://en.wikipedia.org/wiki/Bell_polynomials.
%
%   The Bell polynomials are used by Lilly and Olhede (2008b) to
%   create terms called "instantaneous modulation functions" 
%   quantifying signal time variability.  
%
%   BCELL=BELLPOLY(KCELL) also works, where KCELL is a cell array 
%   containing N elements of identical size.  
%
%   BELLPOLY uses an iterative algorithm.
%
%   See also CUM2MOM, INSTMOM.
%
%   'bellpoly --t' runs some tests.
%
%   Usage: [b1,b2,b3]=bellpoly(k1,k2,k3);
%          bcell=bellpoly(kcell);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2009 J.M. Lilly --- type 'help jlab_license' for details

if strcmpi(varargin{1}, '--t')
     bellpoly_test,return
end

if nargin==1&&iscell(varargin{1})
    cum=varargin{1};
else
    cum=varargin;
end

mom{1}=1;
for n=1:length(cum)
   mom{n+1}=zeros(size(mom{1}));  
   for p=0:n-1
       mom{n+1}=mom{n+1}+choose(n-1,p).*cum{n-p}.*mom{p+1};
   end
end

mom=mom(2:end);

if nargin==1&&iscell(varargin{1})
    varargout{1}=mom;
else
    for i=1:nargout
        varargout{i}=mom{i};
    end
end

function[b]=bellpoly_algebraic(n,k)

switch n
    case 1 
        b=k{1};
    case 2
        b=k{1}.^2+k{2};
    case 3
        b=k{1}.^3+3*k{1}.*k{2}+k{3};
    case 4
        b=k{1}.^4+6*(k{1}.^2).*k{2}+4*k{1}.*k{3}+3*k{2}.^2+k{4};
end




function[b]=bellpoly_determinant(n,kcell)
b=zeros(size(kcell{1}));
for k=1:length(b(:));
    b(k)=bellpoly_innerloop(n,k,kcell);
end

function[b]=bellpoly_innerloop(n,k,kcell)
mat=zeros(n,n);

for i=1:n
    for j=i:n
        mat(i,j)=choose(n-i,j-i).*kcell{j-i+1}(k);
    end
    if i<n
        mat(i+1,i)=-1;
    end
end
b=det(mat);
 

function[]=bellpoly_test

N=10;
kcell=cell(4,1);
for i=1:4
   kcell{i}=randn(N,1);
end

tol=1e-6;
bool1=zeros(4,1);
bool2=zeros(4,1);
bool3=zeros(4,1);

[d1,d2,d3,d4]=bellpoly(kcell{1},kcell{2},kcell{3},kcell{4});
for i=1:4
    b=bellpoly_algebraic(i,kcell);
    c=bellpoly_determinant(i,kcell);
    d=bellpoly(kcell);
    bool1(i)=aresame(b,c,tol);
    bool2(i)=aresame(b,d{i},tol);
    bool3(i)=aresame(eval(['d' int2str(i)]),b,tol);
end


reporttest('BELLPOLY algebraic versus determinant algorithms',allall(bool1))
reporttest('BELLPOLY recursive versus algebraic algorithms, cell array form',allall(bool2))
reporttest('BELLPOLY recursive versus algebraic algorithms, output form',allall(bool2))



