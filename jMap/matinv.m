function[invmat]=matinv(varargin)
%MATINV  Fast inversion of arrays of small matrices.
%
%   MATINV is a low-level function called by POLYSMOOTH.
%  
%   Let MAT be an array of K different M x M matrices A1,A2,...,AK.
%   INV=MATINV(MAT) then returns an array of N inverse matrices.
%    
%   If MAT has dimensions M x M x K1 x ... KK, then MATINV returns an 
%   array of the same size containing the inverses of the M x M matrices.
%
%   For example, MAT could be 4 x 4 x 10 x 10, in which case the inverses
%   of one hundred 4 x 4 matrices are found.
%
%   MAT can have any dimensionality so long as the matrices to be inverted
%   occupy the first two dimensions. 
%
%   Note that MATINV only works matrices with M=2 through M=12.
%   ____________________________________________________________
%
%   Algorithms
%
%   MATINV can use either of two different algorithms.  This is specified
%   with INV=MATINV(MAT,STR). 
%    
%   MATINV(MAT,'direct') uses algebraic expressions for 2 x 2 and 3 x 3 
%   matrix inverses together with Boltz's block diagonal recursion formula.
%   For details, see the following links
%
%        http://mathworld.wolfram.com/MatrixInverse.html
%        http://en.wikipedia.org/wiki/Invertible_matrix.
%
%   MATINV(MAT,'loop') uses Matlab's INV function together with a 
%   straightforward loop.
%
%   The direct algorithm, which is the default, can be much faster when
%   MAT is large and the dimension to be inverted is small.  For arrays of 
%   matrices larger than 8 x 8, the straightforward loop may be faster.
%   ____________________________________________________________
%
%   'matinv --t' runs some tests.
%
%   Usage: inv=matinv(mat);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2022 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    matinv_test,return
end

str='direct';
if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
end

mat=varargin{1};
D=2;

K=size(mat,D);

if K~=size(mat,D-1)
    error('Dimensions D and D-1 of MAT must be the same for it to be invertible.')
end

if K>12
    error('MATINV number of dimensions of matrix should be no more than M=12.')
end

sizemat=size(mat);

if length(sizemat)>2
    %[mat,index]=matinv_strip(mat); 
    
    if strcmpi(str(1:3),'dir')
        invmat=matinv_direct(mat,D,K);
    elseif strcmpi(str(1:3),'loo')
        invmat=matinv_loop(mat,D);
    end

else
    invmat=inv(mat);
end


function[mat,index]=matinv_strip(mat)

sizemat=size(mat);
mat=reshape(mat,prod(sizemat(1:end-2)),sizemat(end-1),sizemat(end));
bool=~isnan(sum(sum(mat,3),2));
index=find(bool);
mat=mat(index,:,:);

function[invmat]=matinv_loop(mat,D)

%[lastmsg,lastid]=lastwarn;
warning('off','MATLAB:illConditionedMatrix');
warning('off','MATLAB:nearlySingularMatrix');
%warning('off','MATLAB:singularMatrix');

sizemat=size(mat);

if D==lnsd(mat)
    mat=reshape(mat,prod(sizemat(1:end-2)),sizemat(end-1),sizemat(end));
    invmat=0*mat;

    for i=1:size(mat,1)
        invmat(i,:,:)=inv(squeeze(mat(i,:,:)));
    end
elseif D==2
    mat=reshape(mat,sizemat(1),sizemat(2),prod(sizemat(3:end)));
    invmat=0*mat;

    for i=1:size(mat,3)
        invmat(:,:,i)=inv(squeeze(mat(:,:,i)));
    end
end

invmat=reshape(invmat,sizemat);
warning('on','MATLAB:illConditionedMatrix');
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:singularMatrix');

function[invmat]=matinv_direct(mat,D,K)

%if K==1
%    invmat=1./mat;
if K==2
    invmat=matinv_twoxtwo(mat,D);
elseif K==3
    invmat=matinv_threexthree(mat,D);
else
    invmat=matinv_block(mat,D,K);
end
    
function[invmat]=matinv_twoxtwo(mat,D)
ac=vindex(mat,1,D);
bd=vindex(mat,2,D);

a=vindex(ac,1,D-1);
c=vindex(ac,2,D-1);
b=vindex(bd,1,D-1);
d=vindex(bd,2,D-1);

deta=a.*d-b.*c;
a=a./deta;
b=b./deta;
c=c./deta;
d=d./deta;

ac=0*ac;
bd=0*bd;

ac=vindexinto(ac,d,1,D-1);
ac=vindexinto(ac,-c,2,D-1);

bd=vindexinto(bd,-b,1,D-1);
bd=vindexinto(bd,a,2,D-1);

invmat=zeros(size(mat));
invmat=vindexinto(invmat,ac,1,D);
invmat=vindexinto(invmat,bd,2,D);

function[mat]=matinv_threexthree(mat,D)

c1=vindex(mat,1,D);
c2=vindex(mat,2,D);
c3=vindex(mat,3,D);

a11=vindex(c1,1,D-1);
a21=vindex(c1,2,D-1);
a31=vindex(c1,3,D-1);

a12=vindex(c2,1,D-1);
a22=vindex(c2,2,D-1);
a32=vindex(c2,3,D-1);

a13=vindex(c3,1,D-1);
a23=vindex(c3,2,D-1);
a33=vindex(c3,3,D-1);

b11=a22.*a33-a23.*a32;
b21=a23.*a31-a21.*a33;
b31=a21.*a32-a22.*a31;

b12=a13.*a32-a12.*a33;
b22=a11.*a33-a13.*a31;
b32=a12.*a31-a11.*a32;

b13=a12.*a23-a13.*a22;
b23=a13.*a21-a11.*a23;
b33=a11.*a22-a12.*a21;

deta=a11.*b11+a12.*b21+a13.*b31;  
%This looks funny because b21 is minus 

b11=b11./deta;
b21=b21./deta;
b31=b31./deta;

b12=b12./deta;
b22=b22./deta;
b32=b32./deta;

b13=b13./deta;
b23=b23./deta;
b33=b33./deta;

clear deta

c1=vindexinto(c1,b11,1,D-1);
c1=vindexinto(c1,b21,2,D-1);
c1=vindexinto(c1,b31,3,D-1);

c2=vindexinto(c2,b12,1,D-1);
c2=vindexinto(c2,b22,2,D-1);
c2=vindexinto(c2,b32,3,D-1);

c3=vindexinto(c3,b13,1,D-1);
c3=vindexinto(c3,b23,2,D-1);
c3=vindexinto(c3,b33,3,D-1);

mat=vindexinto(mat,c1,1,D);
mat=vindexinto(mat,c2,2,D);
mat=vindexinto(mat,c3,3,D);

function[invmat]=matinv_block(mat,D,K)

%This keeps the two matrices about the same size, which seems the 
%fastest option in tests

i1=1:floor(K/2);
i2=floor(K/2)+1:K;

%i1,i2

ac=vindex(mat,i1,D);
bd=vindex(mat,i2,D);

a=vindex(ac,i1,D-1);
c=vindex(ac,i2,D-1);
b=vindex(bd,i1,D-1);
d=vindex(bd,i2,D-1);

%vsize(a,b,c,d)

%Recursion
ainv=matinv(a);

%Note with D=2 I can use pagemtimes, which is considerably faster
dcab=matinv(d-pagemtimes(pagemtimes(c,ainv),b),D);
bnew=-pagemtimes(pagemtimes(ainv,b),dcab);
anew=ainv+pagemtimes(pagemtimes(-bnew,c),ainv);
cnew=-pagemtimes(pagemtimes(dcab,c),ainv);
    
dnew=dcab;

ac=vindexinto(ac,anew,i1,D-1);
ac=vindexinto(ac,cnew,i2,D-1);

bd=vindexinto(bd,bnew,i1,D-1);
bd=vindexinto(bd,dnew,i2,D-1);

invmat=zeros(size(mat));
invmat=vindexinto(invmat,ac,i1,D);
invmat=vindexinto(invmat,bd,i2,D);

function[z]=matmult(mata,matb,K)
%MATMULT  Matrix multiplication for arrays of matrices. 
%
%   Often we wish to form the matrix product A*B for two matrices A and B. 
%   MATMULT performs matrix multiplication for arrays of such matrices.
%
%   Let A be an array of K different M x N matrices, and similarly
%   let B be an array of K different N x P matrices.
%
%   C=MATMULT(A,B,DIM) returns an array of the K products CK=AK*BK, where 
%   DIM gives the dimension in A on which the M x N matrix begins.
%
%   Thus A, B, and C have the following dimensions
%
%                               DIM
%                                |
%         A --   K1 x K2 x ... x M x N x ... x KN
%         B --   K1 x K2 x ... x N x P x ... x KN
%         C --   K1 x K2 x ... x M x P x ... x KN
%
%   The usual matrix multiplication is then C=MATMULT(A,B,1).  
%
%   A and B can have any dimensionality so long as the dimensions
%   containing the matrices are adjacent to each other.
%
%   See also VECTMULT.
%
%   'matmult --t' runs some tests.
%
%   Usage: matc=matmult(mata,matb,dim); 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details        

if strcmpi(mata, '--t')
  matmult_test,return
end

sizea=size(mata);
sizeb=size(matb);

if size(mata,K+1)~=size(matb,K)
    error('Sorry, MATA and MATB do not have the right dimensions for matrix multiplication.')
end

M=size(mata,K);
P=size(matb,K+1);

z=zeros([sizea(1:K) sizeb(K+1:end)]);

evalme='z(';
for i=1:length(sizea)
    if i~=K&&i~=K+1
        evalme=[evalme ':,'];
    elseif i==K
        evalme=[evalme 'i,'];
    elseif i==K+1;
        evalme=[evalme 'j,'];
    end
end
evalme=[evalme(1:end-1) ')=zsub;'];
    
for i=1:M
    suba=vindex(mata,i,K);
    suba=permute(suba,[1:K-1 K+1 K K+2:length(sizea)]);

    for j=1:P
        subb=vindex(matb,j,K+1);        
        zsub=sum(suba.*subb,K); 
        eval(evalme)
    end
end

function[]=matmult_test

matm=[1 3; 4 5];
matn=[2 -7; 3 -4];

mat2a=matmult(matm,matn,1);
mat2b=matm*matn;

reporttest('MATMULT with single 2x2 matrix vs',aresame(mat2a,mat2b,1e-6))


matm=[1 3; 4 5; 2 3];
matn=[2 -7 2;  3 3 -4];

mat2a=matmult(matn,matm,1);
mat2b=matn*matm;

reporttest('MATMULT of 2x3 matrix with 3x2 matrix vs',aresame(mat2a,mat2b,1e-6))

matm=[1 3; 4 5; 2 3];
matn=[2 -7 2;  3 3 -4];

mat2a=matmult(matm,matn,1);
mat2b=matm*matn;

reporttest('MATMULT of 3x2 matrix with 2x3 matrix vs',aresame(mat2a,mat2b,1e-6))


matm=[1 ; 4 ; 2 ];
matn=[2 -7 2;  3 3 -4];

mat2a=matmult(matn,matm,1);
mat2b=matn*matm;

reporttest('MATMULT of 2x3 matrix with 3x1 matrix vs',aresame(mat2a,mat2b,1e-6))

matm=[1 ; 4 ; 2 ]';
matn=[2 -7 2;  3 3 -4]';

mat2a=matmult(matm,matn,1);
mat2b=matm*matn;

reporttest('MATMULT of 1x3 matrix with 3x2 matrix vs',aresame(mat2a,mat2b,1e-6))


matm=randn(10000,2,2);
matn=randn(10000,2,2);

tic
mat=0*matm;
for i=1:size(matm,1)
    mat(i,:,:)=squeeze(matm(i,:,:))*squeeze(matn(i,:,:));
end
t1=toc;

tic 
mat2=matmult(matm,matn,2);
t2=toc;

reporttest('MATMULT with 10000 random 2x2 matrices',aresame(mat,mat2,1e-6))
disp(['MATMULT was ' num2str(t1./t2) ' times faster than loop.'])


matm(1:14:end)=nan;
matn(1:17:end)=nan;
tic
mat=0*matm;
for i=1:size(matm,1)
    mat(i,:,:)=squeeze(matm(i,:,:))*squeeze(matn(i,:,:));
end
t1=toc;

tic 
mat2=matmult(matm,matn,2);
t2=toc;

reporttest('MATMULT with 10000 random 2x2 matrices and nans',aresame(mat,mat2,1e-6))
disp(['MATMULT was ' num2str(t1./t2) ' times faster than loop.'])


function[]=matinv_test
disp('Testing that direct and looping algorithms match for non-singular matrices.')

rng(0);


for n=2:12
    disp(['MATINV testing case of ' int2str(n) 'x' int2str(n) ' matrices.'])

    mat=randn(n,n,10000);
   % mat=randn(n,n,10);
    tic
    invmat=matinv(mat,'loop');
    t1=toc;
    
    tic
    invmat2=matinv(mat,'direct');
    t2=toc;
    
    tol=1e-1;
    bool=false(size(mat,3),1);
    for i=1:size(mat,3)
        bool(i)=rcond(squeeze(mat(:,:,i)))>=tol;
    end
    
    vindex(invmat,invmat2,find(bool),3);
    
    reporttest(['MATINV M=' int2str(n) ' and D=2'],aresame(invmat,invmat2,1e-2))
    disp(['MATINV was ' num2str(t1./t2) ' times faster than INV.'])
end
