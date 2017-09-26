function[varargout] = vmedian(varargin)
%VMEDIAN  Median over finite elements along a specified dimension.
%
%   Y=VMEDIAN(X,DIM) takes the median of all finite elements of X along      
%   dimension DIM. 
%                           
%   [Y1,Y2,...YN]=VMEDIAN(X1,X2,...XN,DIM) also works.
%
%   VMEDIAN(X1,X2,...XN,DIM); with no output arguments overwrites the 
%   original input variables.
%
%   VMEDIAN, like MATLAB's MEDIAN, defines the median over an even 
%   number of values to be the average of the two middle elements.
%
%   VMEDIAN uses a fast algorithm which can be several times faster 
%   than MEDIAN, but unlike MEDIAN, excludes both INFs and NANs.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2015 J.M. Lilly --- type 'help jlab_license' for details        

if strcmpi(varargin{1}, '--t')
  vmedian_test,return
end

dim=varargin{end};

for i=1:length(varargin)-1
   if isreal(varargin{i})
       varargout{i}=vmedian1(varargin{i},dim);
   else
       varargout{i}=vmedian1(real(varargin{i}),dim)+sqrt(-1)*vmedian1(imag(varargin{i}),dim);
   end
end

eval(to_overwrite(nargin-1))

function[med]=vmedian1(data,dim)

numel=sum(isfinite(data),dim);
med=vzeros(size(numel),'nan');

temp=vswap(data,nan,inf);
temp=vswap(temp,-inf,inf);
sorted=sort(temp,dim,'ascend');
sorted=permute(sorted,[1:dim-1 dim+1:ndims(sorted) dim]);
sorted=reshape(sorted,[length(numel(:)),size(sorted,ndims(sorted))]);

ii=(1:length(numel(:)))';
ii=reshape(ii,size(numel));

boolodd=isodd(numel);
booleven=~boolodd&(numel~=0);

indexodd=sub2ind(size(sorted),ii(boolodd),(numel(boolodd)+1)./2);
indexeven1=sub2ind(size(sorted),ii(booleven),numel(booleven)./2);
indexeven2=sub2ind(size(sorted),ii(booleven),numel(booleven)./2+1);

med(boolodd)=sorted(indexodd);
med(booleven)=sorted(indexeven1)./2+sorted(indexeven2)./2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=vmedian_test


N=6;
x1=randn(N,N,N,N,N);
bool=false(5,1);
tol=1e-4;

etime1=0;
etime2=0;
for i=1:5
    tic;y1=vmedian(x1,i);etime1=etime1+toc;
    tic;z1=median(x1,i);etime2=etime2+toc;
    bool(i)=aresame(y1,z1,tol);
end
reporttest('VMEDIAN 5-D with no NANs versus MEDIAN', allall(bool))
disp(['VMEDIAN was ' num2str(etime2./etime1) ' faster than MEDIAN.'])



N=30;
x1=randn(N,N,N);
x1(2:7:end)=nan;

tic
y1=vmedian(x1,1);
y2=vmedian(x1,2);
y3=vmedian(x1,3);
etime1=toc;

z1=zeros(size(y1));
z2=zeros(size(y2));
z3=zeros(size(y3));

tic
for i=1:N
    for j=1:N
        z1(1,i,j)=median(x1(isfinite(x1(:,i,j)),i,j),1);
        z2(i,1,j)=median(x1(i,isfinite(x1(i,:,j)),j),2);
        z3(i,j)=median(x1(i,j,isfinite(x1(i,j,:))),3);
    end
end
etime2=toc;
tol=1e-4;
reporttest('VMEDIAN 3-D including NANs', aresame(y1,z1,tol) && aresame(y2,z2,tol)&&aresame(y3,z3,tol))
disp(['VMEDIAN was ' num2str(etime2./etime1) ' faster than loop with MEDIAN.'])


N=30;
x1=round(randn(N,N,N));
x1(2:7:end)=nan;
x1(1:9:end)=x1(1);

tic
y1=vmedian(x1,1);
y2=vmedian(x1,2);
y3=vmedian(x1,3);
etime1=toc;

z1=zeros(size(y1));
z2=zeros(size(y2));
z3=zeros(size(y3));

tic
for i=1:N
    for j=1:N
        z1(1,i,j)=median(x1(isfinite(x1(:,i,j)),i,j),1);
        z2(i,1,j)=median(x1(i,isfinite(x1(i,:,j)),j),2);
        z3(i,j)=median(x1(i,j,isfinite(x1(i,j,:))),3);
    end
end
etime2=toc;
tol=1e-4;
reporttest('VMEDIAN 3-D including NANs and repeated data', aresame(y1,z1,tol) && aresame(y2,z2,tol)&&aresame(y3,z3,tol))
disp(['VMEDIAN was ' num2str(etime2./etime1) ' faster than loop with MEDIAN.'])



reporttest('VMEDIAN NaNs with row median',aresame([1 nan],vmedian([1 nan; nan nan],1)));
reporttest('VMEDIAN NaNs with column median',aresame([1 nan]',vmedian([1 nan; nan nan],2)));

