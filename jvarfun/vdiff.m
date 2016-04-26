function[varargout]=vdiff(varargin)
%VDIFF	Length-preserving first central difference.
%
%   DX=VDIFF(X,DIM) differentiates X along dimension DIM using the first 
%   central difference; DX is the same size as X.                                 
%                                                                        
%   [D1,D2,...,DN]=VDIFF(X1,X2,...,XN,DIM) for multiple input variables 
%   also works. 
%
%   VDIFF(X1,X2,...,DIM); with no output arguments overwrites the
%   original input variables.
%
%   DXDT=VDIFF(DT,...) optionally uses scalar timestep DT to approximate
%   a time derivative, i.e. DXDT equals DX divided by DT.
%   _____________________________________________________________________
%
%   First and last points
%
%   The first and last points must be treated differently, as the central 
%   difference is not defined there.  Three different methods can be used.
%
%   VDIFF(...,STR) specifies which method to use.
%
%        'endpoint'  uses the first forwards / first backwards difference
%                    at the first and last point, respectively.  
%        'periodic'  treats the array as being periodic along dimension DIM,
%                    so that the central difference is defined at endpoints.
%        'nans'      fills in the first and last values with NANs.
%
%   The default behavior is 'endpoint'.
%   _____________________________________________________________________
%
%   'vdiff --t' runs some tests.
%
%   Usage:  x=vdiff(x,dim);
%           x=vdiff(dt,x,dim);
%           x=vdiff(dt,x,dim,'periodic');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2015 J.M. Lilly --- type 'help jlab_license' for details    
  
%   I am so irritated by diff

 
if strcmpi(varargin{1}, '--t')
  vdiff_test,return
end


if length(varargin{1})==1
    dt=varargin{1};
    varargin=varargin(2:end);
else
    dt=1;
end

if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
else
    str='endpoints';
end
    
n=varargin{end};
varargin=varargin(1:end-1);



for i=1:length(varargin)
    varargout{i}=vdiff1(varargin{i},n,str);
    if dt~=1
       varargout{i}=varargout{i}/dt;
    end
end

eval(to_overwrite(length(varargin)))


function[y]=vdiff1(x,n,str)
  	  
if ~isempty(x)        
    %y=vshift(x,1,n)/2-vshift(x,-1,n)/2;
     
    %/*********************************************************************
    %This block is just a slightly faster way to take the first central difference
    sizex=size(x);
    %Correct for the fact that sizex ignores singleton dimensions
    if n>length(sizex)
        for i=length(sizex)+1:n
            sizex(i)=1;
        end
    end
    sizey=sizex;
    sizey(n)=sizex(n)-2;
    
%     dx=diff(x,1,n)/2;
%     array=zeros(length(size(x)),1);
%     array(n)=-1;
%     y=vindex(dx,1:size(x,n)-2,n)+vindex(dx,2:size(x,n)-1,n);
%  
    y=zeros(sizey);
    i1=1:2:size(x,n);
    if iseven(size(x,n))
        i2=i1+1;
    else
        i2=i1(1:end-1)+1;
    end
  
    %This is just the same as .... y=vshift(x,1,n)/2-vshift(x,-1,n)/2;
    y=vindexinto(y,diff(vindex(x,i1,n),1,n)/2,i1(1:end-1),n)+vindexinto(y,diff(vindex(x,i2,n),1,n)/2,i2(1:end-1),n);
    %\*********************************************************************
        
    if strcmpi(str(1:3),'end')
        dx1=vindex(x,2,n)-vindex(x,1,n);
        dxend=vindex(x,size(x,n),n)-vindex(x,size(x,n)-1,n);
        y=cat(n,dx1,cat(n,y,dxend));
    elseif strcmpi(str(1:3),'nan')
        dx1=nan*vindex(x,1,n);
        y=cat(n,dx1,cat(n,y,dx1));
    elseif strcmpi(str(1:3),'per')
        dx1=(vindex(x,2,n)-vindex(x,1,n))/2+(vindex(x,1,n)-vindex(x,size(x,n),n))/2;
        dxend=(vindex(x,size(x,n),n)-vindex(x,size(x,n)-1,n))/2+(vindex(x,1,n)-vindex(x,size(x,n),n))/2;
        y=cat(n,dx1,cat(n,y,dxend));
    end
else
    y=[];
end
    
function[]=vdiff_test

y1=(1:4)';
y2=2*(1:4)';
[x1,x2]=vdiff(y1,y2,1);
bool=aresame(x1,[1 1 1 1]').*aresame(x2,2*[1 1 1 1]');
reporttest('VDIFF', bool)
vdiff(y1,y2,1);
bool=aresame(y1,[1 1 1 1]').*aresame(y2,2*[1 1 1 1]');
reporttest('VDIFF output overwrite', bool)

dt=pi;
y1=(1:4)';
y2=2*(1:4)';
[x1,x2]=vdiff(pi,y1,y2,1);
bool=aresame(x1,[1 1 1 1]'./dt).*aresame(x2,2*[1 1 1 1]'./dt,1e-10);
reporttest('VDIFF with non-unit time step', bool)

rng(1);
xxx=randn(100000,2);
tic;y=vshift(xxx,1,1)/2-vshift(xxx,-1,1)/2;etime1=toc;
tic;y2=vdiff(xxx,1,'periodic');etime2=toc;
reporttest('VDIFF new and former algorithms match, row dimension, even length', aresame(y,y2,1e-10))
disp(['New algorithm was ' num2str(etime2./etime1) ' times faster than former algorithm.'])

xxx=randn(2,100000);
tic;y=vshift(xxx,1,2)/2-vshift(xxx,-1,2)/2;etime1=toc;
tic;y2=vdiff(xxx,2,'periodic');etime2=toc;
reporttest('VDIFF new and former algorithms match, column dimension, even length', aresame(y,y2,1e-10))
disp(['New algorithm was ' num2str(etime2./etime1) ' times faster than former algorithm.'])

rng(1);
xxx=randn(100001,2);
tic;y=vshift(xxx,1,1)/2-vshift(xxx,-1,1)/2;etime1=toc;
tic;y2=vdiff(xxx,1,'periodic');etime2=toc;
reporttest('VDIFF new and former algorithms match, row dimension, odd length', aresame(y,y2,1e-10))
disp(['New algorithm was ' num2str(etime2./etime1) ' times faster than former algorithm.'])

xxx=randn(2,100001);
tic;y=vshift(xxx,1,2)/2-vshift(xxx,-1,2)/2;etime1=toc;
tic;y2=vdiff(xxx,2,'periodic');etime2=toc;
reporttest('VDIFF new and former algorithms match, column dimension, odd length', aresame(y,y2,1e-10))
disp(['New algorithm was ' num2str(etime2./etime1) ' times faster than former algorithm.'])


