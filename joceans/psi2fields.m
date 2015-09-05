function[cv,zeta,n,s,p]=psi2fields(varargin)
%PSI2FIELDS  Velocity and other fields from the streamfunction. [with P.E. Isachsen]
%
%   CV=PSI2FIELDS(PSI) where PSI is a matrix of streamfunction values, 
%   returns the complex velocity U+iV.
%
%   [CV,ZETA,N,S,P]=PSI2FIELDS(PSI) also returns velocity gradient fields:
%
%          ZETA    - Vorticity               dv/dx-du/dy
%             N    - Normal strain           du/dx-dv/dy
%             S    - Shear strain            dv/dx+du/dy
%             P    - Okubo-Weiss parameter   S^2+N^2-ZETA^2
%
%   PSI should be oriented with the X-dimension along *columns* and the 
%   Y-dimensions along the *rows*.  Both of these must be of even length.
%
%   The output fields are all the same size as PSI.  
%  
%   PSI can have additional dimensions after the first two, for example, 
%   with time in the third dimensions.  Up to four dimensions are possible.
%
%   PSI2FIELDS(DX,PSI) uses the grid spacing DX, in units of kilometers, in
%   computing the gradient fields. By default, DX is set to unity. 
%   
%   PSI is taken to have units of meters squared per second.  The units of
%   CV are given in centimeters per second, N and S have units of inverse
%   seconds, and P has units of inverse seconds squared. 
%   __________________________________________________________________
%
%   Algorithm choice
%
%   PSI2FIELDS can use two different algorithms.
%
%   PSI2FIELDS(PSI,'spectral'), the default hehavior, computes the gradient
%   fields in the spectral domain.  This is particularly appropriate when
%   PSI is doubly periodic.
%
%   PSI2FIELDS(PSI,'spatial') alternately computes the gradient fields
%   in the spatial domain, using first central differences for derivatives.
%
%   PSI2FIELDS(PSI,'arakawa') computes gradients in the spatial domain 
%   with the weighted central difference scheme used by Arakawa (1966).
%
%   PSI2FIELDS(PSI,'spatial',STR) or PSI2FIELDS(PSI,'arakawa',STR) 
%   specifies the boundary conditions for derivatives at the domain edges.  
%
%   The default behavior for spatial derivatives is STR='periodic'.  Other
%   choices are 'mirror' and 'zeros'.  See VDIFF for more details. 
%   __________________________________________________________________
%
%   'psi2fields --t' runs some tests.
%   'psi2fields --f' generates two sample figures.
%
%   Usage: cv=psi2fields(psi);
%          [cv,zeta,N,S,P]=psi2fields(psi);
%          [cv,zeta,N,S,P]=psi2fields(dx,psi);
%          [cv,zeta,N,S,P]=psi2fields(dx,psi,'spatial');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2015 J.M. Lilly and P.E. Isachsen 
%                             --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    psi2fields_test,return
elseif strcmpi(varargin{1}, '--f')
    type makefigs_psi2fields
    makefigs_psi2fields;
    return
end

if length(varargin{1})==1
    dx=varargin{1};
    varargin=varargin(2:end);
else
    dx=1;
end
dx=dx*1000;
psi=varargin{1};

str='spectral';
edgestr='periodic';
algstr='fast';

for i=1:3
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'spe')||strcmpi(varargin{end}(1:3),'spa')||strcmpi(varargin{end}(1:3),'ara')
            str=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'fas')||strcmpi(varargin{end}(1:3),'loo')
            algstr=varargin{end};
        else
            edgestr=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end


if strcmpi(str(1:3),'spe')
    if strcmpi(algstr(1:3),'loo')
        [cv,zeta,n,s,p]=psi2fields_spectral_loop(dx,psi);
    elseif strcmpi(algstr(1:3),'fas')
        %We don't actually need to loop, since fft2 is smart
        %Not really much faster though
        [cv,zeta,n,s,p]=psi2fields_spectral(dx,psi,nargout);
    end
elseif strcmpi(str(1:3),'spa')||strcmpi(str(1:3),'ara')
    [cv,zeta,n,s,p]=psi2fields_spatial(dx,psi,edgestr,str,nargout);
end


function[cv,zeta,n,s,p]=psi2fields_spectral_loop(dx,psi)

for i=1:size(psi,3)
    for j=1:size(psi,4)
        [cv(:,:,i,j),zeta(:,:,i,j),n(:,:,i,j),s(:,:,i,j),p(:,:,i,j)]=psi2fields_spectral(dx,psi(:,:,i,j),5);
    end
end

% Great idea, doesn't actually speed things up 
% function[x]=fft2fast(x)
% 
% M=size(x,1);
% N=size(x,2);
% sizex=size(x);
% L=length(x(:));
%
% x=fft(reshape(x,[M L./M]),[],1);
% x=reshape(x,sizex);
% x=permute(x,[2 1 3]);
% x=fft(reshape(x,[N L./N]),[],1);
% x=reshape(x,sizex);
% x=permute(x,[2 1 3]);


function[cv,zeta,n,s,p]=psi2fields_spectral(dx,psi,na)
%This part was written by Pål Erik Isachsen, with changes by JML

[zeta,n,s,p]=vempty;

if mod(size(psi,1),2)==0
    Lx=size(psi,2);
    Ly=size(psi,1);
    nx=[0:Lx/2-1 -Lx/2:-1]./dx*frac(2*pi,Lx);
    ny=[0:Ly/2-1 -Ly/2:-1]./dx*frac(2*pi,Ly);
    %n=[0:size(psi,1)/2 -size(psi,1)/2+1:-1];
else
    error('Sorry, the size of the PSI must be even in the first two dimensions.')
end

[k,l]=meshgrid(nx,ny);
k=vrep(vrep(k,size(psi,3),3),size(psi,4),4);
l=vrep(vrep(l,size(psi,3),3),size(psi,4),4);
Psi=fft2(psi);

%There is a very small imaginary part coming from handling of the Nyquist
%Psi(ceil(size(Psi,1)/2),:)=0;
%Psi(:,ceil(size(Psi,2)/2))=0;

%vsize(k,l,Psi)

U=-l.*Psi*sqrt(-1);
V=k.*Psi*sqrt(-1);
cv=ifft2(U,'symmetric')+sqrt(-1)*ifft2(V,'symmetric');
cv=100*cv;%Convert to cm/s

if na>1
    Zeta=k.*V*sqrt(-1)-l.*U*sqrt(-1);
    zeta=ifft2(Zeta,'symmetric');
end
if na>2
    N=k.*U*sqrt(-1)-l.*V*sqrt(-1);
    n=ifft2(N,'symmetric');
end
if na>3
    S=k.*V*sqrt(-1)+l.*U*sqrt(-1);
    s=ifft2(S,'symmetric');
end

if na>4
    p=(s.^2+n.^2)-zeta.^2;
end


function[cv,zeta,n,s,p]=psi2fields_spatial(dx,psi,edgestr,schemestr,na)

[zeta,n,s,p]=vempty;
u=-psi2fields_vdiff(dx,psi,1,edgestr,schemestr);
v=psi2fields_vdiff(dx,psi,2,edgestr,schemestr);
cv=100*(u+sqrt(-1)*v);  %Convert to cm/s

if na>1
    zeta=psi2fields_vdiff(dx,v,2,edgestr,schemestr)-psi2fields_vdiff(dx,u,1,edgestr,schemestr);
end
if na>2
    n=psi2fields_vdiff(dx,u,2,edgestr,schemestr)-psi2fields_vdiff(dx,v,1,edgestr,schemestr);
end
if na>3
    s=psi2fields_vdiff(dx,v,2,edgestr,schemestr)+psi2fields_vdiff(dx,u,1,edgestr,schemestr);
end
if na>4
    p=(s.^2+n.^2)-zeta.^2;
end

function[df]=psi2fields_vdiff(dx,f,dim,edgestr,schemestr)
%This is basically just to implement what I think is the Arakawa advection
%scheme.  But does it do the same for vorticity?  

f0ca=f(:,1,:,:);
f0cb=f(:,end,:,:);
f0ra=f(1,:,:,:);
f0rb=f(end,:,:,:);

if strcmpi(schemestr(1:3),'ara')
    if dim==1
        dim2=2;
    elseif dim==2
        dim2=1;
    end
    f=frac(1,4)*(2*f + vshift(f,1,dim2) + vshift(f,-1,dim2));
    if dim==1
        %size(f),dim,size(frac(1,3)*(2*f + vshift(f,1,dim2)))
        %f(:,1,:,:)=0;
        %f(:,end,:,:)=0;
        f(:,1,:,:)=frac(1,3)*(2*f0ca + vshift(f0ca,1,dim2));
        f(:,end,:,:)=frac(1,3)*(2*f0cb + vshift(f0cb,-1,dim2));
        %f(:,end,:,:)=f0(:,end,:,:);
    elseif dim==2
        %f(1,:,:,:)=0;
        %f(end,:,:,:)=0;
        f(1,:,:,:)=frac(1,3)*(2*f0ra + vshift(f0ra,1,dim2));
        f(end,:,:,:)=frac(1,3)*(2*f0rb + vshift(f0rb,-1,dim2));
        %f(end,:,:,:)=f0(end,:,:,:);
    end
end
% dim
% size(f)
% edgestr
df=vdiff(dx,f,dim,edgestr);


function[]=psi2fields_test
psi2fields_test_turbulence
%psi2fields_test_jet

function[]=psi2fields_test_jet
load jetsnapshot
cv=jetsnapshot.cv;
dx=jetsnapshot.x(2)-jetsnapshot.x(1);
[cv2,zeta2]=psi2fields(dx,jetsnapshot.psi,'arakawa','endpoint');
%[cv2,zeta2]=psi2fields(dx,jetsnapshot.psi,'spatial','endpoint');

erru=real(cv(2:end-1,2:end-1)-cv2(2:end-1,2:end-1))./maxmax(abs(cv));
bool(1)=allall(erru<2e-7);

errv=imag(cv(2:end-1,2:end-1)-cv2(2:end-1,2:end-1))./maxmax(abs(cv));
bool(2)=allall(errv<2e-7);

reporttest('PSI2FIELDS Arakawa differentiation gives expected velocity for JETSNAPSHOT to 2 parts in 10^-7',allall(bool))

[cv3,zeta2]=psi2fields(dx,jetsnapshot.psi','arakawa','endpoint');
cv3=-sqrt(-1)*cv3';  
reporttest('PSI2FIELDS Arakawa differentiation transpose for JETSNAPSHOT',aresame(cv2,cv3,1e-10))


% This is leftover code from an earlier bugfix 
% if 0
% use jetsnapshot
% [cv3,zeta2]=psi2fields(dx,jetsnapshot.psi,'spatial','endpoint')
% erru3=real(cv(2:end-1,2:end-1)-cv3(2:end-1,2:end-1))./maxmax(abs(cv));
% errv3=imag(cv(2:end-1,2:end-1)-cv3(2:end-1,2:end-1))./maxmax(abs(cv));
% figure,
% subplot(2,2,1)
% jpcolor(x(2:end-1),y(2:end-1),log10(abs(erru)))
% title('Log10 U error, model vs. Arakawa'),colorbar('south'),caxis([-6 -1])
% subplot(2,2,2)
% jpcolor(x(2:end-1),y(2:end-1),log10(abs(errv)))
% title('Log10 V error, model vs. Arakawa'),colorbar('south'),caxis([-6 -1])
% subplot(2,2,3)
% jpcolor(x(2:end-1),y(2:end-1),log10(abs(erru3)))
% title('Log10 U error, model vs. spatial'),colorbar('south'),caxis([-6 -1])
% subplot(2,2,4)
% jpcolor(x(2:end-1),y(2:end-1),log10(abs(errv3)))
% title('Log10 V error, model vs. spatial'),colorbar('south'),caxis([-6 -1])
% packfig(2,2)
% 
% set(gcf,'paperposition',[1 1 11 4])
% fontsize 12 10 10 10
% %print -dpng uverror.png
% end



function[]=psi2fields_test_turbulence
load qgsnapshot
dx=qgsnapshot.x(2)-qgsnapshot.x(1);

[cv,zeta,N,S,P]=psi2fields(dx,qgsnapshot.psi);
[cv2,zeta2,N2,S2,P2]=psi2fields(dx,qgsnapshot.psi,'spatial');

bool(1)=sum(abs(cv(:)-cv2(:)).^2)<sum(abs(cv(:)).^2)/100;
bool(2)=sum(abs(zeta(:)-zeta2(:)).^2)<sum(abs(zeta(:)).^2)/100;
bool(3)=sum(abs(N(:)-N2(:)).^2)<sum(abs(N(:)).^2)/100;
bool(4)=sum(abs(S(:)-S2(:)).^2)<sum(abs(S(:)).^2)/100;
bool(5)=sum(abs(P(:)-P2(:)).^2)<sum(abs(P(:)).^2)/100;

reporttest('PSI2FIELDS spatial and spectral calculations match to 1% for QGSNAPSHOT',allall(bool))

xx=qgsnapshot.psi;
xx(:,:,2)=-2*qgsnapshot.psi';
tic;[cv,zeta,N,S,P]=psi2fields(dx,xx,'loop');t1=toc;
tic;[cv2,zeta2,N2,S2,P2]=psi2fields(dx,xx,'fast');t2=toc;

bool(1)=sum(abs(cv(:)-cv2(:)).^2)<sum(abs(cv(:)).^2)/1e6;
bool(2)=sum(abs(zeta(:)-zeta2(:)).^2)<sum(abs(zeta(:)).^2)/1e6;
bool(3)=sum(abs(N(:)-N2(:)).^2)<sum(abs(N(:)).^2)/1e6;
bool(4)=sum(abs(S(:)-S2(:)).^2)<sum(abs(S(:)).^2)/1e6;
bool(5)=sum(abs(P(:)-P2(:)).^2)<sum(abs(P(:)).^2)/1e6;

reporttest('PSI2FIELDS loop and loopless calculations match to 1 part in 10^6 for QGSNAPSHOT',allall(bool))

[cv,zeta,N,S,P]=psi2fields(dx,qgsnapshot.psi);
bool=sum(abs(zeta(:)-qgsnapshot.zeta(:)).^2)<sum(abs(zeta(:)).^2)/1e9;
reporttest('PSI2FIELDS matches source to 1 part in 10^9 for QGSNAPSHOT vorticity',allall(bool))

bool=sum(abs(cv(:)-qgsnapshot.cv(:)).^2)<sum(abs(cv(:)).^2)/1e9;
reporttest('PSI2FIELDS matches source to 1 part in 10^9 for QGSNAPSHOT velocity',allall(bool))

bool=sum((N(:)-qgsnapshot.N(:)).^2)<sum(N(:).^2)/1e9;
reporttest('PSI2FIELDS matches source to 1 part in 10^9 for QGSNAPSHOT normal strain',allall(bool))

bool=sum((S(:)-qgsnapshot.S(:)).^2)<sum(S(:).^2)/1e9;
reporttest('PSI2FIELDS matches source to 1 part in 10^9 for QGSNAPSHOT shear strain',allall(bool))

