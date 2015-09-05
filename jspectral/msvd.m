function[d,u1,v1,trS,u2,v2]=msvd(mmat,i2,i3)
%MSVD  Singular value decomposition for polarization analysis.
%
%   MSVD computes the singular value decomposition for multiple-
%   transform polarization analysis at all or selected frequencies.
%
%   [D,U1,V1]=MSVD(W) computes the singular value decomposition of 
%   sub-matrices of the eigentransform W.
%
%   The input matrix W is an J x N x K matrix where 
%
%          J is the number of transforms (frequencies or scales)
%          N is the number of dataset components
%          K is the number of eigentransforms (tapers or wavelets)
%
%   The output is
%
%	       D  --  J x MIN(K,N) matrix of singular values
%	      U1  --  J x N matrix of first left singular vectors
%	      V1  --  J x K matrix of first right singular vectors
%
%   The left singular vectors U1 are the eigenvectors of the spectral
%   matrix, while the right singular vectors V1 are the eigenvectors of 
%   the "structure matrix".  
%
%   W may optionally be of size M x J x N x K, in which case the
%   output arguments are all three-dimensional matrices with sizes
%
%	       D  --  M x J x MIN(K,N) matrix of singular values
%	      U1  --  M x J x N matrix of first left singular vectors
%	      V1  --  M x J x K matrix of first right singular vectors
%
%   [D,U1,V1]=MSVD(W,INDEX) for three-dimensional W optionally 
%   performs the SVD only at the rows of MMAT indicated by INDEX.  The 
%   output arguments then all have LENGTH(INDEX) rows rather than M.
%
%   [D,U1,V1,TR]=MSVD optionally outputs TR, the trace of the spectral 
%   matrix, which is of size J (for 3-D W) or M x J (for 4-D W).  
%
%   MSVD(..., 'quiet') suppresses a progress report.
%	
%   Usage:  [d,u1,v1]=msvd(w);
%           [d,u1,v1]=msvd(w,index);
%           [d,u1,v1,tr]=msvd(w,index);
%
%   'msvd --t' runs a test.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 1993--2014 J.M. Lilly --- type 'help jlab_license' for details        



%   Note that these definitions of the left and 
%   right singular vectors use the convention of Lilly and Olhede (20XX) 
%   and not of Park et al. (1987), who use a transposed version of W. 

if strcmpi(mmat,'--t')
  msvd_test;return
end

index=1:size(mmat,1);
str='loud';

if nargin==3
   str=i3;
   index=i2;
elseif nargin==2
   if ischar(i2)
       str=i2;
   else 
       index=i2;
   end
end

    
if ndims(mmat)==3
   mmat1=mmat(index,:,:);
   mmat=zeros([1 size(mmat)]);
   mmat(1,:,:,:)=mmat1;
else
   mmat=mmat(index,:,:,:);
end

 
[N,J,M,K]=size(mmat);
    
mmat=mmat./sqrt(K*M);

d=zeros(N,J,min(M,K));
u1=zeros(N,J,M);
v1=zeros(N,J,K);

if nargout>3
  u2=zeros(N,J,M);
  v2=zeros(N,J,K);
end

for j=1:J
    if ~strcmpi(str(1:3),'qui')
        disp(['Performing SVD of W at band j=' int2str(j) '.'])
    end
    for n=1:N
       %note--for some reason SVD is much faster than SVDS
       mmattemp=squeeze(mmat(n,j,:,:));
       [utemp tempd vtemp]=svd(mmattemp,0);
       d(n,j,:)=diag(tempd);
       u1(n,j,:)=utemp(:,1);
       v1(n,j,:)=vtemp(:,1);
       if nargout> 3
          u2(n,j,:)=utemp(:,2);
          v2(n,j,:)=vtemp(:,2);
       end
   end
end

d=squeeze(d);
u1=squeeze(u1);
v1=squeeze(v1);
if nargout >3
    u2=squeeze(u2);
    v2=squeeze(v2);
end

if nargout ==4
   trS=vsum(vsum(abs(mmat).^2,4),3);
    if N==1
        trS=permute(trS,[2 1]);
    end
end

function[]=msvd_test

[x,t]=testseries_lillypark_array;
x=x(1:100,:);

N=size(x,1);
M=size(x,2);

%Calculate wavelet matrix
J=50;
K=3;
fs=1./(logspace(log10(20),log10(600),J)');
psi=morsewave(N,K,2,4,fs,'bandpass');


%Compute wavelet transforms
wx=wavetrans(x,psi,'mirror');

[d,u1,v1,tr]=msvd(wx,'quiet');

[N2,J2,K2]=size(d);
[N3,J3,M3]=size(u1);
[N4,J4,K4]=size(v1);
[N5,J5]=size(tr);

bool=aresame([N,J,K],[N2,J2,K2]) && aresame([N,J,K],[N4,J4,K4]) && aresame([N,J,M],[N3,J3,M3]) && aresame([N,J],[N5,J5]);
reporttest('MSVD output matrices have correct sizes, J > 1 case',bool)


[d,u1,v1,tr]=msvd(squeeze(wx(:,1,:,:)),'quiet');

[N2,K2]=size(d);
[N3,M3]=size(u1);
[N4,K4]=size(v1);
N5=size(tr);

bool=aresame([N,K],[N2,K2]) && aresame([N,K],[N4,K4]) && aresame([N,M],[N3,M3]) && aresame([N 1],N5);
reporttest('MSVD output matrices have correct sizes, J = 1 case',bool)


function[]=msvd_test2

[x,t]=testseries_lillypark_array;
N=size(x,1);
M=size(x,2);

%Calculate wavelet matrix
J=50;
K=3;
fs=1./(logspace(log10(2),log10(60),J)');
psi=morsewave(N,K,2,4,fs,'bandpass');

%Compute clean wavelet transforms
wx0=wavetrans(x,psi,'mirror');
[d0,u10,v10]=msvd(wx0);

%Compute noisy wavelet transforms
xn=x+3*randn(size(x));
wx=wavetrans(xn,psi,'mirror');
[d,u1,v1,u2,v2]=msvd(wx);



p1=[1 2 3 3 2.5 1.5]'.*rot(dom*(0:5)'); 
p1=p1./sqrt(p1'*p1);
p1mat=vrep(vrep(reshape(p1,[1 1 6]),J,2),N,1);

u1proj=vsum(u1.*p1mat,3);
u2proj=vsum(u2.*p1mat,3);

trwx0=squeeze(vsum(vsum(abs(wx0).^2,3),4));
trwx=squeeze(vsum(vsum(abs(wx).^2,3),4));

snr0=sqrt(frac(squared(d0(:,:,1)),trwx0-squared(d0(:,:,1))));
snr1=sqrt(frac(squared(d(:,:,1)),trwx-squared(d(:,:,1))));
snr2=sqrt(frac(squared(d(:,:,2)),trwx-squared(d(:,:,1))-squared(d(:,:,2))));

efrac1=frac(squared(d(:,:,1)),trwx);
efrac2=frac(squared(d(:,:,2)),trwx);

figure,
h=wavespecplot(t,x,1./fs,d0(:,:,1),d(:,:,1),d(:,:,2),0.5);
axes(h(2));cax=caxis;axes(h(3));caxis(cax);axes(h(4));caxis(cax);

figure,
h=wavespecplot(t,x(:,1),1./fs,wx0(:,:,1,1),wx0(:,:,1,2),wx0(:,:,1,2),0.5);
axes(h(2));cax=caxis;axes(h(3));caxis(cax);axes(h(4));caxis(cax);

figure,
h=wavespecplot(t,xn(:,1),1./fs,wx(:,:,1,1),wx(:,:,1,2),wx(:,:,1,2),0.5);
axes(h(2));cax=caxis;axes(h(3));caxis(cax);axes(h(4));caxis(cax);

figure,
h=wavespecplot(t,xn,1./fs,abs(u1proj),abs(u2proj),0.5);
axes(h(2));cax=caxis;axes(h(3));caxis(cax);

figure,
h=wavespecplot(t,xn,1./fs,d(:,:,1),snr1,snr2,0.5);
axes(h(3));cax=caxis;axes(h(4));

figure,
h=wavespecplot(t,xn(:,1),1./fs,d(:,:,1),efrac1,efrac2,0.5);
axes(h(3));cax=caxis;axes(h(4));

figure,
h=wavespecplot(t,xn,1./fs,abs(u1proj).*efrac1,abs(u2proj).*efrac2,0.5);
axes(h(2));cax=caxis;axes(h(3));caxis(cax);




figure,
h=wavespecplot(t,xn,1./fs,abs(u1proj),abs(u2proj));
axes(h(2));cax=caxis;axes(h(3));caxis(cax);


function[x,t]=testseries_lillypark_array
[x1,t]=testseries_lillypark;
x1=anatrans(x1);
x1=x1(1:10:end);
t=t(1:10:end);
dom=pi/6;
p1=[1 2 3 3 2.5 1.5]'.*rot(dom*(0:5)');
p1=p1./sqrt(p1'*p1);
x=real(oprod(x1,p1));
x=10*x./maxmax(x);

function[x,t,xo]=testseries_lillypark
t=(1:4096)';
om=zeros(size(t));
x2=0*om;
index=1000:3000;
om(index)=(t(index))*2/100000 - .019999;
x=sin(om.*t);
x2(index)=-sin((t(index)-1600)*2*pi*1/2800);
x=x.*x2;
x(3500:3515)=-1;
x(3516:3530)=1;
t=linspace(-100,100,length(x))';
xo=x2;

