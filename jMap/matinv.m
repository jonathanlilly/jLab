function[invmat]=matinv(mat,str)
%MATINV  Fast inversion of arrays of small matrices.
%
%   MATINV is a low-level function called by POLYSMOOTH.
%  
%   Let MAT be an array of K different M x M matrices A1,A2,...,AK.
%   INV=MATINV(MAT) then returns an array of N inverse matrices.
%    
%   If MAT has dimensions K1 x K2 x .... M x M, then MATINV returns an 
%   array of the same size containing the inverses of the M x M matrices.
%
%   For example, MAT could be 10 x 10 x 4 x 4, in which case the inverses
%   of one hundred 4 x 4 matrices are found.
%
%   MAT can have any dimensionality so long as the matrices to be inverted
%   occupy the last two dimensions.  The last dimension is interpreted as
%   "columns" and the second to last as "rows."
%
%   Note that MATINV only works matrices with M=2 through M=6.
%   ____________________________________________________________
%
%   Algorithms
%
%   MATINV can use either of two different algorithms.   This is specified
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
%   MAT is large and the dimension to be inverted is small.
%   ____________________________________________________________
%
%   See also MATMULT.
%
%   'matinv --t' runs some tests.
%
%   Usage: inv=matinv(mat);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(mat, '--t')
    matinv_test,return
end

ndims=lnsd(mat);
K=size(mat,ndims);
if K~=size(mat,ndims-1)
    error('The last two dimensions of MAT must be the same for it to be invertible.')
end

if K>6
    error('Sorry, MATINV only supports matrices of size M>1 and M<=6.')
end
if nargin==1   
    str='direct';
end

 %size(mat),ndims
 
sizemat=size(mat);
if length(sizemat)>2
    %[mat,index]=matinv_strip(mat); 
    
    if ~isempty(strfind(str,'dir'))
        invmat=matinv_direct(mat,ndims);
    elseif ~isempty(strfind(str,'loo'))
        invmat=matinv_loop(mat);
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



function[invmat]=matinv_loop(mat)

sizemat=size(mat);
mat=reshape(mat,prod(sizemat(1:end-2)),sizemat(end-1),sizemat(end));
invmat=0*mat;

%[lastmsg,lastid]=lastwarn;
warning('off','MATLAB:illConditionedMatrix');
warning('off','MATLAB:nearlySingularMatrix');
%warning('off','MATLAB:singularMatrix');
for i=1:size(mat,1)
     invmat(i,:,:)=inv(squeeze(mat(i,:,:))); 
end

invmat=reshape(invmat,sizemat);
warning('on','MATLAB:illConditionedMatrix');
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:singularMatrix');

function[invmat]=matinv_direct(mat,ndims)
%ndims
%size(mat)

if size(mat,ndims)==1&&size(mat,ndims-1)==1
        invmat=matinv_onexone(mat,ndims);
elseif size(mat,ndims)==2&&size(mat,ndims-1)==2
    invmat=matinv_twoxtwo(mat,ndims);
elseif size(mat,ndims)==3&&size(mat,ndims-1)==3
    invmat=matinv_threexthree(mat,ndims);
elseif (size(mat,ndims)==4&&size(mat,ndims-1)==4)...
        ||(size(mat,ndims)==5&&size(mat,ndims-1)==5)...
        ||(size(mat,ndims)==6&&size(mat,ndims-1)==6)
    invmat=matinv_block(mat,ndims);
end

%function[invmat]=matinv_onexone(mat,ndims)
%invmat=1./mat;    

function[invmat]=matinv_twoxtwo(mat,ndims)
ac=vindex(mat,1,ndims);
bd=vindex(mat,2,ndims);

a=vindex(ac,1,ndims-1);
c=vindex(ac,2,ndims-1);
b=vindex(bd,1,ndims-1);
d=vindex(bd,2,ndims-1);

deta=a.*d-b.*c;
a=a./deta;
b=b./deta;
c=c./deta;
d=d./deta;

ac=0*ac;
bd=0*bd;

ac=vindexinto(ac,d,1,ndims-1);
ac=vindexinto(ac,-c,2,ndims-1);

bd=vindexinto(bd,-b,1,ndims-1);
bd=vindexinto(bd,a,2,ndims-1);

invmat=zeros(size(mat));
invmat=vindexinto(invmat,ac,1,ndims);
invmat=vindexinto(invmat,bd,2,ndims);

function[mat]=matinv_threexthree(mat,ndims)

c1=vindex(mat,1,ndims);
c2=vindex(mat,2,ndims);
c3=vindex(mat,3,ndims);

a11=vindex(c1,1,ndims-1);
a21=vindex(c1,2,ndims-1);
a31=vindex(c1,3,ndims-1);

a12=vindex(c2,1,ndims-1);
a22=vindex(c2,2,ndims-1);
a32=vindex(c2,3,ndims-1);

a13=vindex(c3,1,ndims-1);
a23=vindex(c3,2,ndims-1);
a33=vindex(c3,3,ndims-1);

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

c1=vindexinto(c1,b11,1,ndims-1);
c1=vindexinto(c1,b21,2,ndims-1);
c1=vindexinto(c1,b31,3,ndims-1);

c2=vindexinto(c2,b12,1,ndims-1);
c2=vindexinto(c2,b22,2,ndims-1);
c2=vindexinto(c2,b32,3,ndims-1);

c3=vindexinto(c3,b13,1,ndims-1);
c3=vindexinto(c3,b23,2,ndims-1);
c3=vindexinto(c3,b33,3,ndims-1);

mat=vindexinto(mat,c1,1,ndims);
mat=vindexinto(mat,c2,2,ndims);
mat=vindexinto(mat,c3,3,ndims);

function[invmat]=matinv_block(mat,ndims)

if size(mat,ndims)==4
    i1=1:2;
    i2=3:4;
elseif size(mat,ndims)==5
    i1=1:3;
    i2=4:5;
elseif size(mat,ndims)==6
    i1=1:3;
    i2=4:6;
end


ac=vindex(mat,i1,ndims);
bd=vindex(mat,i2,ndims);

a=vindex(ac,i1,ndims-1);
c=vindex(ac,i2,ndims-1);
b=vindex(bd,i1,ndims-1);
d=vindex(bd,i2,ndims-1);

%Recursion
ainv=matinv(a);

dcab=matinv(d-matmult(matmult(c,ainv,ndims-1),b,ndims-1));

bnew=-matmult(matmult(ainv,b,ndims-1),dcab,ndims-1);
anew=ainv+matmult(matmult(-bnew,c,ndims-1),ainv,ndims-1);
cnew=-matmult(matmult(dcab,c,ndims-1),ainv,ndims-1);
dnew=dcab;


ac=vindexinto(ac,anew,i1,ndims-1);
ac=vindexinto(ac,cnew,i2,ndims-1);

bd=vindexinto(bd,bnew,i1,ndims-1);
bd=vindexinto(bd,dnew,i2,ndims-1);

invmat=zeros(size(mat));
invmat=vindexinto(invmat,ac,i1,ndims);
invmat=vindexinto(invmat,bd,i2,ndims);

function[]=matinv_test
disp('Testing that direct and looping algorithms match for non-singular matrices.')
matinv_test_twoxtwo
matinv_test_threexthree
matinv_test_fourxfour
matinv_test_fivexfive
matinv_test_sixxsix


function[]=matinv_test_twoxtwo
mat=[1 3; 4 5];
invmat=inv(mat);
invmat2=matinv(mat);

reporttest('MATINV with single 2x2 matrix',aresame(invmat,invmat2,1e-6))

mat=randn(10000,2,2);

tic
invmat=matinv(mat,'loop');
t1=toc;

tic 
invmat2=matinv(mat,'direct');
t2=toc;

tol=1e-4;
bool=false(size(mat,1),1);
for i=1:size(mat)
    bool(i)=rcond(squeeze(mat(i,:,:)))>=tol;
end

vindex(invmat,invmat2,find(bool),1);

reporttest('MATINV with 10000 random 2x2 matrices',aresame(invmat,invmat2,1e-6))
disp(['MATINV was ' int2str(t1./t2) ' times faster than INV.'])

function[]=matinv_test_threexthree
mat=[1 3 7; -3 4 5; 6 2 -9];
invmat=inv(mat);
invmat2=matinv(mat);

reporttest('MATINV with single 3x3 matrix',aresame(invmat,invmat2,1e-6))

%mat=permute(mat,[3 1 2]);
%mat=vrep(mat,10000,1);
mat=randn(10000,3,3);

tic
invmat=matinv(mat,'loop');
t1=toc;

tic 
invmat2=matinv(mat,'direct');
t2=toc;

tol=1e-3;
bool=false(size(mat,1),1);
for i=1:size(mat)
    bool(i)=rcond(squeeze(mat(i,:,:)))>=tol;
end

vindex(invmat,invmat2,find(bool),1);

reporttest('MATINV with 10000 random 3x3 matrices',aresame(invmat,invmat2,1e-3))
disp(['MATINV was ' int2str(t1./t2) ' times faster than INV.'])


function[]=matinv_test_fourxfour
mat=[1 3 7 -3 ; -3 4 5 10 ; 6 2 -9 2; 1 9 -3 7];
invmat=inv(mat);
invmat2=matinv(mat);

reporttest('MATINV with single 4x4 matrix',aresame(invmat,invmat2,1e-3))

mat=randn(10000,4,4);

tic
invmat=matinv(mat,'loop');
t1=toc;

tic 
invmat2=matinv(mat,'direct');
t2=toc;

tol=1e-3;
bool=false(size(mat,1),1);
for i=1:size(mat)
    bool(i)=rcond(squeeze(mat(i,:,:)))>=tol;
end

vindex(invmat,invmat2,find(bool),1);

reporttest('MATINV with 10000 random 4x4 matrices',aresame(invmat,invmat2,1e-3))
disp(['MATINV was ' num2str(t1./t2) ' times faster than INV.'])



function[]=matinv_test_fivexfive
mat=[1 3 7 -3 17; -3 4 5 10 2; 4 6 2 -9 2;-5 1 9 -3 7; -18 5 9 -4 10];

invmat=inv(mat);
invmat2=matinv(mat);

reporttest('MATINV with single 5x5 matrix',aresame(invmat,invmat2,1e-3))
mat=randn(10000,5,5);
tic
invmat=matinv(mat,'loop');
t1=toc;

tic 
invmat2=matinv(mat,'direct');
t2=toc;

tol=1e-1;
bool=false(size(mat,1),1);
for i=1:size(mat)
    bool(i)=rcond(squeeze(mat(i,:,:)))>=tol;
end

vindex(invmat,invmat2,find(bool),1);

reporttest('MATINV with 10000 random 5x5 matrices',aresame(invmat,invmat2,1e-2))
disp(['MATINV was ' num2str(t1./t2) ' times faster than INV.'])



function[]=matinv_test_sixxsix
mat=[1 3 7; -3 4 5; 6 2 -9];
mat=[mat mat' ; mat flipud(mat)];
invmat=inv(mat);
invmat2=matinv(mat);

reporttest('MATINV with single 6x6 matrix',aresame(invmat,invmat2,1e-3))
mat=randn(10000,6,6);
tic
invmat=matinv(mat,'loop');
t1=toc;

tic 
invmat2=matinv(mat,'direct');
t2=toc;

tol=1e-1;
bool=false(size(mat,1),1);
for i=1:size(mat)
    bool(i)=rcond(squeeze(mat(i,:,:)))>=tol;
end

vindex(invmat,invmat2,find(bool),1);

reporttest('MATINV with 10000 random 6x6 matrices',aresame(invmat,invmat2,1e-2))
disp(['MATINV was ' num2str(t1./t2) ' times faster than INV.'])


