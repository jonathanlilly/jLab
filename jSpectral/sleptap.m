function[v,lambda]=sleptap(varargin)
%SLEPTAP  Calculate Slepian tapers.
%
%   [PSI,LAMBDA]=SLEPTAP(M,P,K) calculates the K lowest-order Slepian
%   tapers PSI of length M and time-bandwidth product P, together with
%   their eigenvalues LAMBDA. PSI is M x K and LAMBDA is K x 1.
%
%   K is optional and defaults to 2P-1.  
%   P is optional and defaults to 4.
%   
%   For M<=512, SLEPTAP uses the tridiagonal method described in Percival 
%   and Walden (1993).  For M>512, it first computes tapers for M=512 and 
%   then spline-interpolates.
%   
%   M may also be an array of lengths.  In this case PSI is a cell array of 
%   matrices, with PSI{1} being M(1) x K, PSI{2} being M(2) x K, etc., 
%   while LAMBDA is K x LENGTH(M).  See 'Cell array input' under MSPEC.
%   _____________________________________________________________________
%
%   Normalization
%
%   By default, the tapers are set to have unit energy. Alternatively
%   SLEPTAP(...,'bandpass') uses the "bandpass" normalization in which the
%   tapers are rescaled so that the maximum value of the Fourier transform
%   of the first taper is set to two. 
%
%   See WAVETRANS for details on bandpass normalization.  
%   _____________________________________________________________________
%
%   Parallelization
%
%   SLEPTAP(M, ...,'parallel') when M is an array of lengths parallelizes 
%   the taper computation using a PARFOR loop.  This requires that Matlab's
%   Parallel Computing Toolbox be installed.
%   _____________________________________________________________________
%
%   See also MSPEC, MSVD, TWOSPECPLOT.
%   
%   'sleptap --t' runs some tests.  
%
%   Usage:  [psi,lambda]=sleptap(n); 
%           [psi,lambda]=sleptap(n,p,k); 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2017 J.M. Lilly --- type 'help jlab_license' for details        

if strcmpi(varargin{1}, '--t')
  sleptap_test,return
end
n=varargin{1};
str='energy';
cores='serial';

for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'par')||strcmpi(varargin{end}(1:3),'ser')
            cores=varargin{end};
        else
        str=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end
na=length(varargin);

if strcmpi(cores(1:3),'par')
    if exist('parpool')~=2
        disp('Sorry, parallel algorithm requires the Parallel Computing Toolbox.')
        disp('Defaulting to the standard algorithm.')
        cores='serial';
    end
end



%Nmax=256;
Nmax=512;
if na==1
	p=4;
else
    p=varargin{2};
end

if na<=2
	k=2*p-1;
else
    k=varargin{3};
end


%if n>Nmax, interpolate
if anyany(n>=Nmax)
    [v1,d]=sleptap_one(Nmax,p,k);
end
     
if strcmpi(cores(1:3),'ser')||(length(n)==1)  
    for j=1:length(n)
        if length(n)>1
            disp(['SLEPTAP computing tapers for time series #' int2str(j) ' of ' int2str(length(n)) '.'])
        end
        if n(j)==Nmax
            v{j,1}=v1;
        elseif n(j)<=Nmax
            v{j,1}=sleptap_one(n(j),p,k);
        elseif n(j)>Nmax
            disp(['SLEPTAP interpolating to length ',int2str(n(j)),'.'])
            for i=1:k
                vi=sinterp(v1(:,i),n(j));
                vi=vi./sqrt(vi'*vi);
                v{j,1}(:,i)=vi;
            end
        end
    end
    if nargout==2
        lambda=zeros(k,length(n));
        for j=1:length(n)
            if n(j)>=Nmax
                lambda(:,j)=sleptap_lambda_one(Nmax,p,k,v1);
            else
                lambda(:,j)=sleptap_lambda_one(n(j),p,k,v{j});
            end
        end
    end
else
    %Exactly the same as the above but with two changes to parfors
    parfor j=1:length(n)  %Parfor
        if length(n)>1
            disp(['SLEPTAP computing tapers for time series #' int2str(j) ' of ' int2str(length(n)) '.'])
        end
        if n(j)==Nmax
            v{j,1}=v1;
        elseif n(j)<=Nmax
            v{j,1}=sleptap_one(n(j),p,k);
        elseif n(j)>Nmax
            disp(['SLEPTAP interpolating to length ',int2str(n(j)),'.'])
            for i=1:k
                vi=sinterp(v1(:,i),n(j));
                vi=vi./sqrt(vi'*vi);
                if vi(round(end/2))<0
                    vi=-vi;
                end
                v{j,1}(:,i)=vi;
            end
        end
    end
    if nargout==2
        lambda=zeros(k,length(n));
        parfor j=1:length(n)  %Parfor
            if n(j)>=Nmax
                lambda(:,j)=sleptap_lambda_one(Nmax,p,k,v1);
            else
                lambda(:,j)=sleptap_lambda_one(n(j),p,k,v{j});
            end
        end
    end
end






if length(v)==1
    v=v{1};
end

%Bandpass normalization
if strfind(str,'ban')
    Vmax=max(abs(fft(v(:,1))));
    v=v*frac(2,Vmax);
end

function[v,d]=sleptap_one(n,p,k)
disp(['SLEPTAP calculating tapers of length ',int2str(n),'.'])


w=p./n;

%taper calculation using tridiagonal matrix
%tic 
mat=zeros(n,n);
index=(1:n+1:n*n);
mat(index)=((n-1-2*(0:n-1))./2).^2.*cos(2*pi*w);
index2=index(1:length(index)-1)+1;
index3=index(2:length(index))-1;
mat(index2)=(1:n-1).*(n-(1:n-1))./2;
mat(index3)=(1:n-1).*(n-(1:n-1))./2;
%toc

%tic
%t=[0:n-1]';
%tmat=osum(t,-t);
%mat=frac(sin(2*pi*w*tmat),pi*tmat);
%vswap(mat,nan,2*w);
%toc

OPTIONS.disp=0;
OPTIONS.maxit=2000;
%OPTIONS.tol=1e-8;

[v,d]=eigs(mat,k,'LA',OPTIONS);
%[v,d]=eigs(double(mat),max(2*p-1,k),'LM',OPTIONS);
%v=v(:,1:k);d=d(1:k);


for i=1:size(v,2)
    if v(round(end/2),i)<0
        v(:,i)=-v(:,i);
    end
end


function[lambda]=sleptap_lambda_one(n,p,k,v)
w=p./n;
i=(0:n-1)'*ones(1,n);
j=i';
A=pi*(i-j);
index=find(A==0);
A(index)=1;
A=sin(2*pi*w*(i-j))./A;
A(index)=2*w;
%find eigenvalues
for i=1:k
    %note unexplained matlab quirk: dividing
    %two column vectors gives you a column vector,
    %dividing two row vectors gives you a scalar
    lambda(i,1)=(A*v(:,i))'/v(:,i)';
end

%/***************************************************
%here's some garbage
if 0
%see how much spline-interpolated ones vary from others
xx=v(:,1);xx=xx/sqrt(xx'*xx);figure,plot(xx)
xx=diff(xx);xx=xx/sqrt(xx'*xx);hold on,plot(xx,'g')
xx=diff(xx);xx=xx/sqrt(xx'*xx);plot(xx,'r')
xx=diff(xx);xx=xx/sqrt(xx'*xx);plot(xx,'c');

for i=1:4
v(:,i)=v(:,i)/sqrt(v(:,i)'*v(:,i));
end

figure,plot(v)


l1=lambda;
v1=v;


k=4;
n=256;
w=4/n;
mat=zeros(n,n);
index=(1:n+1:n*n);
mat(index)=((n-1-2*(0:n-1))./2).^2.*cos(2*pi*w);
index2=index(1:length(index)-1)+1;
index3=index(2:length(index))-1;
mat(index2)=(1:n-1).*(n-(1:n-1))./2;
mat(index3)=(1:n-1).*(n-(1:n-1))./2;
[v,d]=eigs(mat,k);
%[d,index]=sort(diag(d));
%v=v(:,index);
%d=flipud(d);
%v=fliplr(v);
for i=1:4
v(:,i)=v(:,i)/sqrt(v(:,i)'*v(:,i));
end

A=calcdefmat(n,w);

for i=1:k
	lambda(i,1)=mean((A*v(:,i))\v(:,i));
end
v=v(:,1:k);
 
for i=1:4
	v1i(:,i)=interp1((1:100)/100',v1(:,i),(1:256)'/256,'cubic');
end
for i=1:4
v1i(:,i)=v1i(:,i)/sqrt(v1i(:,i)'*v1i(:,i));
end
end
%end garbage
%\***************************************************



function[a]=calcdefmat(n,w)
i=(0:n-1)'*ones(1,n);
j=i';
a=sin(2*pi*w*(i-j))./(pi*(i-j));
a(isnan(a))=2*w;

function [yi] = sinterp(y,n2)
%SINTERP  Spline-interpolates a column vector to a new length.
%
%   YI=SINTERP(Y,NI) spline-interpolates the column vector Y to the
%   new length NI.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information (C)
%   1993, 2004 J.M. Lilly --- type 'help jlab_license' for details

n1=length(y);
x=zeros(1,n1);
xi=zeros(1,n2);
x=(1:n1)';
xi=(1:(n1-1)/(n2-1):n1)';
yi=interp1(x,y,xi,'spline');  %This is much better than linear

function[]=sleptap_test

[psi,lambda]=sleptap(200);

tol=1e-6;
bool=false(1,size(psi,2));
for j=1:size(psi,2)       
        bool(1,j)=aresame(vsum(psi(:,j).^2,1),1,tol);
end
reporttest('SLEPTAP unit energy',allall(bool))

[psi,lambda]=sleptap([200 512 1024]);

tol=1e-6;
bool=false(length(psi),size(psi{1},2));
for i=1:length(psi)
    for j=1:size(psi{1},2)       
        bool(i,j)=aresame(vsum(psi{i}(:,j).^2,1),1,tol);
    end
end

reporttest('SLEPTAP unit energy with interpolation & cell output',allall(bool))

[psi,lambda]=sleptap(200,4,1,'bandpass');

reporttest('SLEPTAP bandpass normalization',aresame(maxmax(abs(fft(psi))),2,1e-10))
