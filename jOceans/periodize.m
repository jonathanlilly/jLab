function[varargout]=periodize(varargin)
%PERIODIZE  Returns a doubly periodic version of an input array.
%
%   [XP,YP,FP]=PERIODIZE(X,Y,F) takes the matrix F and forms a doubly-
%   periodic version FP, which will be 3*SIZE(F). 
%
%   XP and YP are the new X and Y axes extended based on the delta-X and
%   delta-Y values within X and Y.  X and Y must both have uniform spacing.
%   
%   [XE,YE,FP]=PERIODIZE(N,X,Y,F) only extends F periodically by N points, 
%   so that FP is SIZE(F,1)+2*N by SIZE(F,2)+2*N.
%
%   [XE,YE,FP]=PERIODIZE(N,M,X,Y,F) extends F periodically by N points in
%   the X-direction, corresponding to the *columns* of F, and by M points 
%   in the Y-direction, corresponding to the *rows* of F. 
%
%   Thus in this case FP will be SIZE(F,1)+2*M by SIZE(F,2)+2*N.
%
%   F may have multiple higher dimensions, for example a third dimension 
%   corresponding to time.  Only the first two dimensions are periodized.
%
%   For singly periodic versions, use PERIODIZE(0,M,F) or PERIODIZE(N,0,F). 
%
%   [XP,YP,F1P,F2P,...,FKP]=PERIODIZE(...,X,Y,F1,F2,...,FK) also works,
%   and returns periodized versions of all the K arrays F1, F2,...,FK.
%
%   'periodize --t' runs a test.
%   'periodize --f' generates a sample figure.
%
%   Usage: [xp,yp,fp]=periodize(x,y,f);
%          [xp,yp,fp]=periodize(N,x,y,f);
%          [xp,yp,fp]=periodize(N,M,x,y,f);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2015 J.M. Lilly --- type 'help jlab_license' for details
 

if strcmpi(varargin{1}, '--t')
    periodize_test,return
end
if strcmpi(varargin{1}, '--f')
    type makefigs_periodize
    makefigs_periodize;
    return
end

if length(varargin{1})==1
    N=varargin{1};
    varargin=varargin(2:end);
else
    N=[];
end
if length(varargin{1})==1
    M=varargin{1};
    varargin=varargin(2:end);
else
    M=[];
end

x=varargin{1};
y=varargin{2};
varargin=varargin(3:end);


for i=1:length(varargin)
    if i==1
        [varargout{i+2},varargout{1},varargout{2}]=periodize_one(N,M,x,y,varargin{i});
    else
       varargout{i+2}=periodize_one(N,M,x,y,varargin{i});
    end
end


function[zp,xp,yp]=periodize_one(N,M,x,y,z)

if isempty(M)&&isempty(N) 
    N=length(x);
    M=length(y);
elseif isempty(M)&&~isempty(N)
    M=N;
end

zp=zeros(size(z,1)+2*M,size(z,2)+2*N,size(z,3),size(z,4));

%Edges ordered clockwise from top
e=z(1:M,:,:,:);
f=z(:,end-N+1:end,:,:);
g=z(end-M+1:end,:,:,:);
h=z(:,1:N,:,:);

zp(M+1:M+size(z,1),N+1:N+size(z,2),:,:)=z;
zp(1:M,N+1:size(z,2)+N,:,:)=g;
zp(end-M+1:end,N+1:size(z,2)+N,:,:)=e;
zp(M+1:size(z,2)+M,1:N,:,:)=f;
zp(M+1:size(z,2)+M,end-N+1:end,:)=h;

%Corners ordered clockwise from top left
a=z(1:M,1:N,:,:);
b=z(1:M,end-N+1:end,:,:);
c=z(end-M+1:end,end-N+1:end,:,:);
d=z(end-M+1:end,1:N,:,:);

zp(1:M,1:N,:,:)=c;
zp(1:M,end-N+1:end,:,:)=d;
zp(end-M+1:end,end-N+1:end,:,:)=a;
zp(end-M+1:end,1:N,:,:)=b;

%That actually worked the first time

if nargout>1
    xp=zeros(length(x)+2*N,1);
    yp=zeros(length(x)+2*M,1);
    
    xp(N+1:N+size(z,1))=x;
    yp(M+1:M+size(z,2))=y;
    
    dx=x(2)-x(1);
    dy=y(2)-y(1);
    
    xp(1:N)=dx*[1:N]-(N+1)*dx+x(1);
    xp(end-N+1:end)=dx*[1:N]+x(end);
    
    yp(1:M)=dy*[1:M]-(M+1)*dy+y(1);
    yp(end-M+1:end)=dy*[1:M]+y(end);
end

function[]=periodize_test
load qgsnapshot

[xp,yp,fp]=periodize(100,200,qgsnapshot.x,qgsnapshot.y,qgsnapshot.psi);

bool(1)=all(abs(diff(xp)-(xp(2)-xp(1)))./(xp(2)-xp(1))<1e-4);
bool(2)=all(abs(diff(yp)-(yp(2)-yp(1)))./(xp(2)-xp(1))<1e-4);
bool(3)=~anyany((abs(vdiff(fp,1,'nan'))./vstd(qgsnapshot.psi(:),1))>0.12);
bool(4)=~anyany((abs(vdiff(fp,2,'nan'))./vstd(qgsnapshot.psi(:),1))>0.12);
reporttest('PERIODIZE no discontinutities in XE, YE, or FP',allall(bool))

[xp,yp,fp1,fp2]=periodize(100,200,qgsnapshot.x,qgsnapshot.y,qgsnapshot.psi,qgsnapshot.psi);

reporttest('PERIODIZE multiple input/output periodizations',aresame(fp1,fp2))

