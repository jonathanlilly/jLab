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
%   then spline-interpolates.  (Tests show spline interpolation is far
%   superior to linear interpolation for this problem.)   
%
%   The tapers are normalized to have unit energy.
%   _____________________________________________________________________
%
%   Computing multiple taper lengths simultaneously
%
%   M may also be an array of lengths.  In this case PSI is a cell array of 
%   matrices, with PSI{1} being M(1) x K, PSI{2} being M(2) x K, etc., 
%   while LAMBDA is again K x 1.
%
%   By default, SLEPTAP will down-interpolate the M=512 tapers to all 
%   shorter lengths, an approximation that is very good for data lengths
%   greater than say M=64, and resulting in a vast computational savings. 
%   Because short data segments are difficult to extract reliable spectral
%   information from anyway, this should be sufficient for most purposes.
%
%   Alternatively, SLEPTAP(...,'exact') will directly compute solutions of 
%   PSI and LAMBDA for entries m for which M(m)<=512.  In this case, LAMBDA
%   will be a K x M matrix.  This algorithm, like the default behavior, 
%   will spline-interpolate PSI from the M=512 solution to larger M values.
%
%   SLEPTAP(...,'exact','parallel') will parallelize this computation using
%   a PARFOR loop.  This requires Matlab's Parallel Computing Toolbox.
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
%   (C) 2000--2021 J.M. Lilly --- type 'help jlab_license' for details        

if strcmpi(varargin{1}, '--t')
  sleptap_test,return
end
n=varargin{1};
cores='serial';
computestr='interpolate';

for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'par')||strcmpi(varargin{end}(1:3),'ser')
            cores=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'exa')||strcmpi(varargin{end}(1:3),'int')
            computestr=varargin{end};
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

if length(n)==1
    computestr='exact';
end


if length(n)==1
    [v,lambda]=sleptap_one(min(n,Nmax),p,k,n);
else
    %In this case we have a cell array input
    %--------------------------------------------------------------------------
    %this code computes all tapers for a multi-component dataset using a
    %one-time interpolation, much faster than looping
    
    [v1,lambda1]=sleptap_one(Nmax,p,k);
    x=[0:size(v1,1)-1]'./(size(v1,1)-1);
    xnew=cell(length(n),1);
    for i=1:length(n)
        xnew{i}=[0:n(i)-1]'./(n(i)-1);
    end
    vnew=cell(size(v1,2),1);
    for j=1:size(v1,2) %these are the K eigenvectors
        vnew{j}=interp1(x,v1(:,j),cell2col(xnew),'spline');
        vnew{j}(isnan(cell2col(xnew)))=nan;
        vnew{j}=col2cell(vnew{j});
        for i=1:length(vnew{j})%normalize to unit energy
            vnew{j}{i}=vnew{j}{i}./sqrt(sum(squared(vnew{j}{i})));
        end
    end
    %need to reorganize to be a cell array of matrices
    vnew_former=vnew;
    vnew=vnew{1};%that's the first eigenvector
    for j=2:size(v1,2)
        for i=1:length(vnew)
            vnew{i}(:,j)=vnew_former{j}{i};
        end
    end
    
    %compute exactly if requested, otherwise use interpolation
    %---------------------------------------------------------------------
    lambda=zeros(k,length(n));
    if strcmpi(cores(1:3),'ser')
        if strcmpi(computestr(1:3),'exa')
            for j=1:length(n)
                disp(['SLEPTAP computing tapers for time series #' int2str(j) ' of ' int2str(length(n)) '.'])
                if n(j)==Nmax
                    v{j,1}=v1;
                    lambda(:,j)=lambda1;
                elseif n(j)<Nmax
                    v{j,1}=sleptap_one(n(j),p,k);
                    lambda(:,j)=sleptap_lambda_one(n(j),p,k,v{j,1});
                elseif n(j)>Nmax
                    v{j,1}=vnew{j};
                    lambda(:,j)=lambda1;
                end
            end
        else
            v=vnew;
            lambda=lambda1;
        end
    else
        %Exactly the same as the above but with two changes to parfors
        if strcmpi(computestr(1:3),'exa')
            parfor j=1:length(n)
                disp(['SLEPTAP computing tapers for time series #' int2str(j) ' of ' int2str(length(n)) '.'])
                if n(j)==Nmax
                    v{j,1}=v1;
                    lambda(:,j)=lambda1;
                elseif n(j)<Nmax
                    v{j,1}=sleptap_one(n(j),p,k);
                    lambda(:,j)=sleptap_lambda_one(n(j),p,k,v{j,1});%wow, need the ",1" for Matlab's parfor to work
                elseif n(j)>Nmax
                    v{j,1}=vnew{j};
                    lambda(:,j)=lambda1;
                end
            end
        else
            v=vnew;
            lambda=lambda1;
        end
    end
end


function[v,lambda]=sleptap_one(n,p,k,nnew)
disp(['SLEPTAP calculating tapers of length ',int2str(n),'.'])

%n,p,k,nnew

w=p./n;

%taper calculation using tridiagonal matrix
%see SAPA Section 8.3 
%you have to calculate the eigenvalues separately 

mat=zeros(n,n);
index=(1:n+1:n*n);
mat(index)=((n-1-2*(0:n-1))./2).^2.*cos(2*pi*w);
index2=index(1:length(index)-1)+1;
index3=index(2:length(index))-1;
mat(index2)=(1:n-1).*(n-(1:n-1))./2;
mat(index3)=(1:n-1).*(n-(1:n-1))./2;
%toc

OPTIONS.disp=0;
OPTIONS.maxit=2000;
%OPTIONS.tol=1e-8;

[v,~]=eigs(mat,min(k,size(mat,1)),'LA',OPTIONS);

for i=size(mat,1)+1:k
    v(:,i)=nan;
end

for i=1:size(v,2)
    if v(floor(end/2),i)<0
        v(:,i)=-v(:,i);
    end
end

lambda=sleptap_lambda_one(n,p,k,v);

if nargin==4
    if n~=nnew
        x=[0:size(v,1)-1]'./(size(v,1)-1);
        xnew=[0:nnew-1]'./(nnew-1);
        vnew=zeros(nnew,k);
        for j=1:size(v,2)
            vnew(:,j)=interp1(x,v(:,j),xnew,'spline');
            vnew(:,j)=vnew(:,j)./sqrt(sum(squared(vnew(:,j))));
        end
        v=vnew;
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

% function[a]=calcdefmat(n,w)
% i=(0:n-1)'*ones(1,n);
% j=i';
% a=sin(2*pi*w*(i-j))./(pi*(i-j));
% a(isnan(a))=2*w;

function[]=sleptap_test

[psi,lambda]=sleptap(200);

tol=1e-6;
bool=false(1,size(psi,2));
for j=1:size(psi,2)       
        bool(1,j)=aresame(vsum(psi(:,j).^2,1),1,tol);
end
reporttest('SLEPTAP unit energy',allall(bool))

%testing bulk calculation 
N=[400:77:800];%just make some weird lengths
clear psi1 lambda1 
for i=1:length(N)
     [psi1{i,1},lambda1{i,1}]=sleptap(N(i),4);
end
[psi2,lambda2]=sleptap(N,4);
[psi3,lambda3]=sleptap(N,4,'exact');
b1=aresame(psi1,psi2,1e-3)&&aresame(psi1,psi3,1e-3);
b2=aresame(lambda1{end},lambda2)&&aresame(lambda1{end},lambda3(:,end));
reporttest('SLEPTAP computing multiple tapers simultaneously',b1&&b2)


function[]=sleptap_interpolation_test

%comparing spline vs. pchip vs. linear interpolation

N1=512;
N2=2*N1;%doubling experiment
N2=N1./2%halving experiment
[psi,lambda]=sleptap(N1); 
[psi2,lambda]=sleptap(N2); 

x1=[0:N1-1]'./(N1-1);
x2=[0:N2-1]'./(N2-1);

psi2l=interp1(x1,psi,x2,'linear');
psi2s=interp1(x1,psi,x2,'spline');
psi2p=interp1(x1,psi,x2,'pchip');

for i=1:size(psi,2)
    psi2(:,i)=psi2(:,i)./sqrt(sum(squared(psi2(:,i))));
    psi2l(:,i)=psi2l(:,i)./sqrt(sum(squared(psi2l(:,i))));
    psi2s(:,i)=psi2s(:,i)./sqrt(sum(squared(psi2s(:,i))));
    psi2p(:,i)=psi2p(:,i)./sqrt(sum(squared(psi2p(:,i))));
end

projl=abs(squeeze(vsum(psi2l.*psi2,1)))
projs=abs(squeeze(vsum(psi2s.*psi2,1)))
projp=abs(squeeze(vsum(psi2p.*psi2,1)))

figure,
subplot(1,3,1),plot(psi2l-psi2)
subplot(1,3,2),plot(psi2s-psi2)
subplot(1,3,3),plot(psi2p-psi2)
%spline interpolation is much better than the other two for doubling
%for halving, it doesn't matter, they are all the same


