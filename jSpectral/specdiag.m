function[l1,l2,th,nu,z]=specdiag(varargin)
% SPECDIAG  Diagonalize a 2 x 2 spectral matrix.
%  
%   [D1,D2,TH,NU,Z]=SPECDIAG(S) for a 2 x 2 spectral matrix S returns the
%   eigenvalues D1 and D2 and the angles TH and NU and eigenvector matrix Z
%   which diagonalize S according to the decomposition
%
%          Z = JMAT2(TH) *  KMAT(PI/4) * JMAT2(NU) 
%          S = Z [D1 0; 0 D2] Z'
%
%   where JMAT2 is the two by two rotation matrix and KMAT(PI/4) is a
%   diagonal matrix with main diagonal [(1+i) (1-i)]/SQRT(2). 
%
%   If S is a 2 x 2 x M matrix, then D1, D2, NU, and TH are all Mx1 column
%   vectors, and Z is 2 x 2 x M. More generally, if S is 2 x 2 x M x ... N, 
%   then D1, D2, NU, and TH are all M x ... N  and Z is 2 x 2 x M x ... N.  
%   Note that singleton dimensions will be squeezed out.
%  
%   [...]=SPECDIAG(S11,S22,S12) also works. In this case all input and
%   output variables are matrices of the same size.  Z is not output. 
%
%   'SPECDIAG --t' runs some tests.
%  
%   Usage: [d1,d2,th,nu,z]=specdiag(s);
%          [d1,d2,th,nu]=specdiag(s11,s22,s12);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2016 J.M. Lilly --- type 'help jlab_license' for details        

if strcmpi(varargin{1},'--t')
   specdiag_test;return
end

breshape=0;
if nargin==1
  s=varargin{1};
  if ndims(s)>3
     sizes=size(s);
     sizes=sizes(3:end);
     breshape=1;
     s=reshape(s,[2,2,prod(sizes)]);
  end
elseif nargin==3
  if ndims(varargin{1})>1
     sizes=size(varargin{1});
     breshape=1;
     for i=1:nargin
       varargin{i}=varargin{i}(:);
     end
  end 
  s(1,1,:)=varargin{1};
  s(2,2,:)=varargin{2};
  s(1,2,:)=varargin{3};
  s(2,1,:)=conj(varargin{3});
  if nargout>4
    error('Z cannot be output with S11, S22, S12 input format.')
  end
end
  
[e,p,alpha,beta]=polparams(s);

alphamax=sqrt(alpha.^2+real(beta).^2);
twoth=atan(frac(real(beta),alpha));

 
index=find(alpha.*cos(twoth)+real(beta).*sin(twoth)<0);
if ~isempty(index)
  twoth(index)=twoth(index)+pi;
end
th=squeeze(angle(rot(twoth))/2);

if ~isreal(s)
    twonu=atan(frac(imag(beta),alphamax));
    index=find(alphamax.*cos(twonu)+imag(beta).*sin(twonu)<0);
    if ~isempty(index)
        twonu(index)=twonu(index)+pi;
    end
    nu=squeeze(angle(rot(twonu))/2);
else
    nu=0*th;
end


c=real(s(2,1,:));
d=imag(s(2,1,:));
a=s(1,1,:);
b=s(2,2,:);
%figure,plot([squeeze(a) squeeze(b) squeeze(c) squeeze(d)] )
dets=a.*b-c.^2-d.^2;%figure,plot(abs(squeeze(dets)))
trs=a+b;

%ßdets(dets<0)=0;  %Sometimes I get negative values of numerical noise


l1=squeeze(frac(1,2).*(trs+sqrt(trs.^2-4*dets)));
l2=squeeze(frac(1,2).*(trs-sqrt(trs.^2-4*dets)));


if nargout>4
  z1=jmat2(th);
  z2=kmat(pi/4);
  z3=jmat2(nu);
  z=0*z1;
  for i=1:size(z1,3)
    z(:,:,i)=z1(:,:,i)*z2*z3(:,:,i);
  end
  z=squeeze(z);
end


if breshape
  l1=squeeze(reshape(l1,sizes));
  l2=squeeze(reshape(l2,sizes));
  th=squeeze(reshape(th,sizes));
  nu=squeeze(reshape(nu,sizes));
  if nargout>4
     z=squeeze(reshape(z,[2 2 sizes]));
  end
end


function[]=specdiag_test
M=100;
s=randspecmat(M);
[l1,l2,th,rho,z]=specdiag(s);

[l1b,l2b,thb,rhob]=specdiag(s(1,1,:),s(2,2,:),s(1,2,:));
tol=1e-6;
b=aresame([l1,l2,th,rho],[l1b,l2b,thb,rhob],tol);
reporttest('SPECDIAG two input formats match',b)

for i=1:M
  [z2,d]=eigs(s(:,:,i));  %In 2015 eigs is flipping the eigenvalue matrix
  d1=max(d(1),d(4));
  d2=min(d(1),d(4));
  d=[d1 0 ; 0 d2];
  b1(i,1)=aresame(l1(i),d1,tol) && aresame(l2(i),d2,tol);
  b2(i,1)=aresame(d,z(:,:,i)'*s(:,:,i)*z(:,:,i),tol);  
end

%disp(['Testing ' int2str(M) ' random iterations'])
reporttest('SPECDIAG eigenvalues match Matlab', b1)
reporttest('SPECDIAG eigenvector matrix diagonalizes S', b2)

P=frac(l1-l2,l1+l2);
dets=squeeze(s(1,1,:).*s(2,2,:)-s(1,2,:).*s(2,1,:));
detreals=squeeze(s(1,1,:).*s(2,2,:)-real(s(1,2,:)).*real(s(2,1,:)));
trs=squeeze(s(1,1,:)+s(2,2,:));
P2=sqrt(1-4.*frac(dets,trs.^2));
reporttest('SPECDIAG P lambda vs det/tr', aresame(P,P2,tol))

%/********************************************************
disp('Alpha, beta, and gamma')
alpha=squeeze(frac(s(1,1,:)-s(2,2,:),s(1,1,:)+s(2,2,:)));
beta=squeeze(frac(2.*s(1,2,:),s(1,1,:)+s(2,2,:)));
gamma=squeeze(frac(s(1,2,:),sqrt(s(1,1,:).*s(2,2,:))));

b3=aresame(alpha.^2+abs(beta).^2,P.^2,tol);
reporttest('SPECDIAG alpha beta P', b3)

alpha2=sqrt(frac(P.^2-abs(gamma).^2,1-abs(gamma).^2));
gamma2=sqrt(frac(P.^2-alpha.^2,1-alpha.^2));

index=find(abs(1-abs(gamma).^2)>.00001);
b4=aresame(abs(alpha(index)),alpha2(index),tol);
reporttest('SPECDIAG alpha in terms of gamma', b4)

index=find(abs(1-abs(alpha).^2)>.00001);
b5=aresame(1./abs(gamma(index)),1./gamma2(index),tol);
reporttest('SPECDIAG gamma in terms of alpha', b5)

if 0
thm=squeeze(1/2*atan(frac(2.*real(s(1,2,:)),s(1,1,:)-s(2,2,:))));

for i=1:M
  stemp=jmat2(thm(i))*s(:,:,i)*jmat2(thm(i))';
  if stemp(1,1)<stemp(2,2);
      stemp=jmat2(thm(i)+pi/2)*s(:,:,i)*jmat2(thm(i)+pi/2)';
  end
  
  sp(:,:,i)=stemp;
  %This is rotated into normal coordinates
end
alphap=squeeze(frac(sp(1,1,:)-sp(2,2,:),sp(1,1,:)+sp(2,2,:)));
alphap=real(alphap);  %some small numerical noise
alphap2=sqrt(1-4.*frac(detreals,trs.^2));
b6=aresame(alphap,alphap2,tol); 
reporttest('Maximum alpha in terms of det/tr', b6)

alphap2=sqrt(P.^2-4.*frac(imag(squeeze(s(1,2,:))).^2,trs.^2));
b7=aresame(alphap,alphap2,tol); 
reporttest('Maximum alpha in terms of Q/tr', b7)

b8=all(abs(real(sp(1,2,:)))<tol);
reporttest('Imaginary cross-terms', b8)
%\********************************************************
end


function[S,C]=randspecmat(m)
%RANDSPECMAT  Generates random 2x2 spectral matrices for testing.
%
%   S=RANDSPECMAT generates a random 2 x 2 spectral matrix.  S is Hermitian
%   with complex-valued off-diagnonal elements, with a minimum determinant 
%   of zero.
%  
%   S=RANDSPECMAT(M) generates M such matrices, with S being 2x2xM.
%
%   [S,C]=RANDSPECMAT returns both a random spectral matrix and a random 
%   complementary spectral matrix.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details        

a=randn(m,1).^2;  
b=randn(m,1).^2;  
c=(2*rand(m,1)-1+2*sqrt(-1).*rand(m,1)-sqrt(-1)).*sqrt(a.*b./2);

S=zeros(2,2,m);
S(1,1,:)=a;
S(2,2,:)=b;
S(2,1,:)=c;
S(1,2,:)=conj(c);

if nargout ==2
  a=randn(m,1).^2.*rot(rand(m,1)*2*pi);
  b=randn(m,1).^2.*rot(rand(m,1)*2*pi);
  c=randn(m,1)+sqrt(-1).*randn(m,1);

  C=zeros(2,2,m);
  C(1,1,:)=a;
  C(2,2,:)=b;
  C(2,1,:)=c;
  C(1,2,:)=conj(c);
end


function[x]=kmat(theta)
%KMAT 2x2 phase shift matrix through specified angle.
%
%   K=KMAT(PHI) creates a phase shift matrix
%
%       K=[e.^(i*PHI)     0;
%          0     e.^(-i*PHI)]
%
%   such that K*X phase-shifts the first element of the column vector X
%   forward by PHI radians, and phase-shfits the second the second
%   element backwards by PHI radians.
%
%   If LENGTH(PHI)>1, then K will have dimension 2 x 2 x LENGTH(PHI).
%
%   See also JPOLY, VECTMULT, JMAT2, IMAT, and TMAT.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2010 J.M. Lilly --- type 'help jlab_license' for details        

if length(theta)>1;
  theta=theta(:);
end

x=zeros(2,2,length(theta));
x(1,1,:)=rot(theta);
x(2,2,:)=rot(-theta);


if size(theta,1)~=length(theta(:));
    x=reshape(x,[2,2,size(theta)]);
end

