function[varargout]=vectmult(varargin)
%VECTMULT  Matrix multiplication for arrays of vectors. 
%
%   [Y1,Y2,...,YM]=VECTMULT(A,X1,X2,...,XN) is equivalent to the usual
%   multiplication of a matrix by a vector,
%  
%           Y = A * X
%
%   where A is an M x N matrix, X is an N-vector, and Y is an M-vector.
%   X1--XN and Y1--YM are the components of X and Y respectively.
%
%   Often however one has arrays of vectors.  Here we group the N 
%   components of X into N arrays, X1,X2,...,XN, which can be of any size
%   provided they are all the same size.
%
%   VECTMULT then returns the components of the output vector Y as the M
%   arrays Y1,Y2,...,YM.
%
%   The matrix A may be a M x N matrix, or alternatively, it may be an 
%   array of size M x N x SIZE(XN).
%
%   See also MATMULT.
%   
%   Usage: [y1,y2]=vectmult(a,x1,x2);
%          [y1,y2,y3]=vectmult(a,x1,x2,y3);
%
%   'vectmult --t' runs a test.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details


if strcmpi(varargin{1},'--t')
  vectmult_test; nnsd_test;return
end

A=varargin{1};
varargin=varargin(2:end);
N=length(varargin);
M=nargout;

size1=size(varargin{1});

bool=true;
for i=2:N
    bool=bool.*aresame(size1,size(varargin{i}));
end
if ~bool
    error('Input arrays must be all the same size.')
end

nnsd1=nnsd(varargin{1});
if size(A,2)~=N
    error('Number of columns in A must equal number of vector components input.')
end
    
for i=1:M
   varargout{i}=zeros(size1);
   for j=1:N
       
       %All those colons since if I write A(i,j,:) Matlab will collapse
       %remaining dimensions
       
       varargout{i}=varargout{i}+squeeze(A(i,j,:,:,:,:,:,:,:)).*varargin{j};
   end
end


function[varargout]=nnsd(varargin)
%NNSD  Number of nonsingleton dimensions
%  
%   N=NNSD(X) returns the number of nonsingleton dimensions of X.
%   Unlike Matlab's NDIMS, which thinks that an array and a scalar
%   both have dimension 2, NNSD defines the dimension of a scalar to
%   be zero and that of an array to be one.  Singleton dimensions are
%   never counted.
% 
%   [N1,N2,...NM]=NNSD(X1,X2,...XM) returns the dimensions of multiple 
%   input arguments.  If zero or one output arguments are given, a
%   single row array [N1 N2 ...NM] is output.  
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002,2004 J.M. Lilly --- type 'help jlab_license' for details        
  
  
if strcmpi(varargin{1}, '--t')
  nnsd_test,return
end


for i=1:nargin
  x=varargin{i};
  sx=size(x);
  nd(i)=length(sx)-length(find(sx==1));
  varargout{i}=nd(i);
end

if nargout==0 || nargout==1
   varargout{1}=nd;
end

function[]=nnsd_test

x=nnsd(1,(1:10),[ (1:10)' (1:10)']);
bool(1)=aresame(x,[0 1 2]);
%disp('Should be 0 1 2')

x=(1:10)';
z(:,:,3)=x;
q(:,:,:,1)=z;
q(:,:,:,2)=z;

x=nnsd(x,z,q);
bool(2)=aresame(x,[1 2 3]);
%disp('Should be 1 2 3')
reporttest('NNSD',all(bool))


function[]=vectmult_test
tol=1e-10;
M=100;
th=rand(M,1)*2*pi;

xo=randn(M,2);
xo1=squeeze(xo(:,1));
xo2=squeeze(xo(:,2));

[xf1,xf2]=vectmult(jmat2(th(1)),xo1,xo2);
xfo(:,1)=cos(th(1)).*xo1-sin(th(1)).*xo2;
xfo(:,2)=sin(th(1)).*xo1+cos(th(1)).*xo2;

b=aresame(squeeze(xfo(:,1)),xf1,tol) && aresame(squeeze(xfo(:,2)),xf2,tol);
reporttest('VECTMULT 2x2 with single random rotation matrix',b)

[xf1,xf2]=vectmult(jmat2(th),xo1,xo2);
xfo=0*xo;
for i=1:M
  xfo(i,:)=jmat2(th(i))*conj(xo(i,:)');
end
b=aresame(squeeze(xfo(:,1)),xf1,tol) && aresame(squeeze(xfo(:,2)),xf2,tol);
reporttest('VECTMULT 2x2 with array of random rotation matrices',b)


tol=1e-10;
M=100;
th=rand(M,1)*2*pi;

xo=randn(M,3);
xo1=squeeze(xo(:,1));
xo2=squeeze(xo(:,2));
xo3=squeeze(xo(:,3));

rotdim=ceil(rand(M,1)*3);
[xf1,xf2,xf3]=vectmult(jmat3(th,rotdim),xo1,xo2,xo3);

xfo=0*xo;
for i=1:M
  xfo(i,:)=jmat3(th(i),rotdim(i))*conj(xo(i,:)');
end

b=aresame(squeeze(xfo(:,1)),xf1,tol) && aresame(squeeze(xfo(:,2)),xf2,tol)&& aresame(squeeze(xfo(:,3)),xf3,tol);

reporttest('VECTMULT 3x3 with array of random rotation matrices',b)

mat=jmat3(th,rotdim);
mat=mat(:,1:2,:);

[xf1,xf2,xf3]=vectmult(mat,xo1,xo2);

xfo=0*xo;
for i=1:M
  xfo(i,:)=squeeze(mat(:,:,i))*conj(xo(i,1:2)');
end

b=aresame(squeeze(xfo(:,1)),xf1,tol) && aresame(squeeze(xfo(:,2)),xf2,tol);

reporttest('VECTMULT 3x2 with array of random rotation matrices',b)


