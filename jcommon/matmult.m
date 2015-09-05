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



