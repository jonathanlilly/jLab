function[varargout]=mom2cum(varargin)
%MOM2CUM  Convert moments to cumulants. 
%
%   [K0,K1,...KN]=MOM2CUM(M0,M1,...MN) converts the first N moments 
%   M0,M1,...MN into the first N cumulants K0,K1,...KN.
%
%   The MN and KN are all scalars or arrays of the same size.
%
%   Note for a probability density function, M0=1 and K0=0.
%
%   KCELL=MOM2CUM(MCELL), where MCELL is a cell array whose (N+1)th
%   is the Nth moment, returns a similar cell array of cumulants.
%
%   MOM2CUM is inverted by CUM2MOM.
%
%   See also CUM2MOM, BELLPOLY.
%
%   'mom2cum --t' runs a test.
%
%   Usage: [k0,k1,k2]=mom2cum(m0,m1,m2);
%          kcell=mom2cum(mcell);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    mom2cum_test,return
end

if nargin==1&&iscell(varargin{1})
    mom=varargin{1};
else
    mom=varargin;
end
warning('off','MATLAB:log:logOfZero')
cum{1}=log(mom{1});
warning('on','MATLAB:log:logOfZero')

for n=1:length(mom)-1
   coeff=zeros(size(cum{1}));  
   for k=1:n-1
       coeff=coeff+choose(n-1,k-1).*cum{k+1}.*frac(mom{n-k+1},mom{1});
   end
   cum{n+1}= frac(mom{n+1},mom{1}) - coeff;
end

if nargin==1&&iscell(varargin{1})
    varargout{1}=cum;
else
    varargout=cum;
end

function[]=mom2cum_test
 
m=randn(4,1);
k=0*m;
k(1)=m(1);
k(2)=m(2)-m(1).^2;
k(3)=2*m(1).^3-3*m(1).*m(2)+m(3);
k(4)=-6*m(1).^4+12*m(1).^2.*m(2)-3.*m(2).^2-4*m(1).*m(3)+m(4);

[k0,k1,k2,k3,k4]=mom2cum(1,m(1),m(2),m(3),m(4));
reporttest('MOM2CUM array input',aresame([k1 k2 k3 k4]',k,1e-6))

mcell{1}=1;
mcell{2}=m(1);
mcell{3}=m(2);
mcell{4}=m(3);
mcell{5}=m(4);

kcell{1}=0;
kcell{2}=k(1);
kcell{3}=k(2);
kcell{4}=k(3);
kcell{5}=k(4);

k=mom2cum(mcell);
bool=false(5,1);
for i=1:length(mcell)
    bool(i)=aresame(k{i},kcell{i},1e-6);
end
reporttest('MOM2CUM cell array input',allall(bool))

[m0a,m1a,m2a,m3a,m4a]=cum2mom(k0,k1,k2,k3,k4);
reporttest('MOM2CUM inverted by CUM2MOM, array input',aresame([1;m],[m0a m1a m2a m3a m4a]',1e-6))

mcella=cum2mom(kcell);
bool=false(5,1);
for i=1:length(kcell)
    bool(i)=aresame(mcella{i},mcell{i},1e-6);
end

reporttest('MOM2CUM inverted by CUM2MOM, cell array input',allall(bool))


ga1=(1:.1:11);
be1=(1:.1:10);
[ga,be]=meshgrid(ga1,be1);

[m0,n0]=morsemom(0,ga,be);
[m1,n1]=morsemom(1,ga,be);
[m2,n2]=morsemom(2,ga,be);
[m3,n3]=morsemom(3,ga,be);
[m4,n4]=morsemom(4,ga,be);

[k0,k1,k2,k3,k4]=mom2cum(m0,m1,m2,m3,m4);
[l0,l1,l2,l3,l4]=mom2cum(1,m1./m0,m2./m0,m3./m0,m4./m0);

%skew=k3./(k2.^(3/2));
%kurt=k4./(k2.^2);

bool=false(4,1);
tol=1e-6;
bool(1)=aresame(k1,l1,tol);
bool(2)=aresame(k2,l2,tol);
bool(3)=aresame(k3,l3,tol);
bool(4)=aresame(k4,l4,tol);

reporttest('MOM2CUM for nonzero zeroth moment using Morse wavelets',allall(bool))
