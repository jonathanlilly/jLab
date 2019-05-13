function[Q1,Q2,Q3,Q4,Q]=divgeom(varargin)
%DIVGEOM  Geometric decomposition of eddy vorticity flux divergence.
%
%   [F1,F2,F3,F4]=DIVGEOM(DX,DY,K,L,THETA) returns the geometric decomposition
%   of eddy vorticity flux divergence associated with variance ellipses 
%   having kinetic energy K, anisotropy L, and orientation THETA.
%
%   K, L, and THETA are matrices of the same size.  These are defined on an
%   X-Y grid with x oriented in *columns* and y oriented in *rows*.  DX and
%   DY are the sampling intervals in the X and Y directions, respectively.
%
%   Note that K and L are related to ellipse parameters KAPPA and LAMBDA
%   used elsewhere in JLAB by K=KAPPA^2 and L=LAMBDA*KAPPA^2.
%
%   F1, F2, F3, and F4 are four different contributions to the eddy 
%   vorticity flux divergence, as follows:
% 
%       F1     Quadratic variations in the linear energy L
%       F2     Product of linear variations in THETA and L
%       F3     Quadratic variations in the orientation THETA
%       F4     Product of linear variations in orientation THETA
%
%   For details, see Waterman and Lilly (2015), Geometric decomposition of
%   eddy-mean flow feedbacks in barotropic systems, J. Phys. Oceanogr. 
%
%   [F1,F2,F3,F4,F]=DIVGEOM(...) also returns the total eddy flux 
%   divergence F, calculated directly, with F1+F2+F3+F4 = F apart from 
%   numerical error.
%
%   By default, DIVGEOM calculates derivates with repeated applications of
%   a first central difference.  DIVGEOM(...,'arakawa') alternately uses a 
%   modified first central difference appropriate for models that employ an
%   Awakawa advection scheme. 
%   __________________________________________________________________
%
%   Usage: [f1,f2,f3,f4]=divgeom(dx,dy,K,L,theta);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2015 J.M. Lilly --- type 'help jlab_license' for details

if strcmpi(varargin{1}, '--f')
    %   'divgeom --f' generates a sample figure. XX not currently working
    %divgeom_figure,return
end

%Divgeom does run tests, but these are hidden because it involves mat-files
%not distributed as a part of JLAB
%if strcmpi(varargin{1},XXX)
%    divgeom_test,return
%end




%   Do I do this? 
%   Filtering
%  
%   DIVGEOM can optionally filter the second-order derivative terms, F1,
%   F3, and F to reduce small-scale noise.
%
%   [F1,F2,F3,F4,F]=DIVGEOM(...,N) smooths F1, F3, and F with an N point
%   boxcar filter in both the X and Y directions.  N should be odd.
%   __________________________________________________________________
%

diffstr='cartesian';
str='second';

for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'car')||strcmpi(varargin{end}(1:3),'ara')
            diffstr=varargin{end};
        elseif strcmpi(varargin{end}(1:3),'dir')||strcmpi(varargin{end}(1:3),'sec')
            str=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end

dx=varargin{1};
dy=varargin{2};
K=varargin{3};
L=varargin{4};
theta=varargin{5};

if length(varargin)==6
    N=varargin{6};
else
    N=0;
end

%[x,y] = meshgrid([1:size(K,2)]*dx,[1:size(K,1)]*dx);


%/************************************************
%Compute gradients
Lx=divgeom_vdiff(dx,L,2,diffstr);
Ly=divgeom_vdiff(dy,L,1,diffstr);

thetax=divgeom_vdiff(dx,frac(1,2)*unwrap(2*theta,[],2),2,diffstr);
thetay=divgeom_vdiff(dy,frac(1,2)*unwrap(2*theta,[],1),1,diffstr);

cos2x=divgeom_vdiff(dx,cos(2*theta),2,diffstr);
cos2y=divgeom_vdiff(dy,cos(2*theta),1,diffstr);

sin2x=divgeom_vdiff(dx,sin(2*theta),2,diffstr);
sin2y=divgeom_vdiff(dy,sin(2*theta),1,diffstr);

M=L.*cos(2*theta);
N=L.*sin(2*theta);

Nx=divgeom_vdiff(dx,N,2,diffstr);
Ny=divgeom_vdiff(dy,N,1,diffstr);

QM=divgeom_mixederiv(dx,dy,divgeom_vdiff(dx,M,2,diffstr),divgeom_vdiff(dy,M,1,diffstr),diffstr);
QN=divgeom_vdiff(dx,Nx,2,diffstr)-divgeom_vdiff(dy,Ny,1,diffstr);

Q=QN-QM;
%\************************************************

% 
% cx=L;
% cx=vdiff(dx,cx,2)+sqrt(-1)*vdiff(dx,cx,1);
% cx=vdiff(dx,cx,2)+sqrt(-1)*vdiff(dx,cx,1);
% Q1=-imag(rot(-2*theta).*cx);

    
Q1a=-cos(2*theta).*divgeom_mixederiv(dx,dy,Lx,Ly,diffstr);
Q1b= sin(2*theta).*(divgeom_vdiff(dx,Lx,2,diffstr)-divgeom_vdiff(dy,Ly,1,diffstr));
Q1=Q1a+Q1b;

%vsize(L,theta,thetax,thetay)
Q3a= 2*L.*cos(2*theta).*(divgeom_vdiff(dx,thetax,2,diffstr)-divgeom_vdiff(dy,thetay,1,diffstr));
Q3b= 2*L.*sin(2*theta).*divgeom_mixederiv(dx,dy,thetax,thetay,diffstr);
Q3=Q3a+Q3b;


if findstr(str,'dir')
    Q2a= 4*cos(2*theta).*(thetax.*Lx-thetay.*Ly);
    Q2b= 4*sin(2*theta).*(thetax.*Ly+thetay.*Lx);
    Q2=Q2a+Q2b;
    
    Q4a= 4*L.*cos(2*theta).*(2*thetax.*thetay);
    Q4b= -4*L.*sin(2*theta).*(thetax.^2-thetay.^2);
    Q4=Q4a+Q4b;
else
    dkdsin2=L.*(divgeom_vdiff(dx,sin2x,2,diffstr)-divgeom_vdiff(dy,sin2y,1,diffstr));
    dldcos2=L.*divgeom_mixederiv(dx,dy,divgeom_vdiff(dx,cos(2*theta),2,diffstr),...
        divgeom_vdiff(dy,cos(2*theta),1,diffstr),diffstr);
    
    Q2a=QN-Q1b-dkdsin2;
    Q2b=-(QM+Q1a-dldcos2);
    Q2=Q2a+Q2b;
    
    Q4a=-dldcos2-Q3b;
    Q4b=dkdsin2-Q3a;
    Q4=Q4a+Q4b;
end


function[df]=divgeom_vdiff(dx,f,dim,schemestr);
%This is basically just to implement what I think is the way the Arakawa
%advection scheme takes derivatives.  Taken from PSI2FIELDS

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
df=vdiff(dx,f,dim);


function[gxy]=divgeom_mixederiv(dx,dy,fx,fy,diffstr)
%Correction for two possible forms of mixed derivative

gxy= divgeom_vdiff(dy,fx,1,diffstr);
gyx= divgeom_vdiff(dx,fy,2,diffstr);
gxy=gxy+gyx;

%bool=abs(gyx)<abs(gxy);
%gxy(bool)=gyx(bool);
%gxy=2*gxy;

function[]=divgeom_test
load highresjet
use highresjet

[q1,q2,q3,q4,q]=divgeom(dx,K,L,theta);
reporttest('DIVGEOM Q1--Q4 add to total Q',aresame(q1+q2+q3+q4,q,1e-16))

[q1d,q2d,q3d,q4d,qd]=divgeom(dx,K,L,theta,'direct');

rat=abs(frac(q2-q2d,q2+q2d));
rat=sort(rat(:));
index=find(log10(rat)>-1,1,'first');

%No more than 15% of the data has the error ratio > 1/10
reporttest('DIVGEOM Q2 direct and implicit forms agree',index/length(rat)>0.85);

rat=abs(frac(q4-q4d,q4+q4d));
rat=sort(rat(:));
index=find(log10(rat)>-1,1,'first');

%No more than 15% of the data has the error ratio > 1/10
reporttest('DIVGEOM Q4 direct and implicit forms agree',index/length(rat)>0.85);



function[]=divgeom_figure
load jetellipses_highres
use jetellipses
 
[q1,q2,q3,q4,q]=divgeom(x(2)-x(1),x(2)-x(1),kappabar,lambdabar,thetabar,5);

figure
subplot(2,2,1),jpcolor(x,y,q1),axis equal,axis tight
subplot(2,2,2),jpcolor(x,y,q2),axis equal,axis tight
subplot(2,2,3),jpcolor(x,y,q3),axis equal,axis tight
subplot(2,2,4),jpcolor(x,y,q4),axis equal,axis tight
for i=1:4
    subplot(2,2,i),%caxis([-.2 .2]/16)
    hold on,%contour(x,y,lambdabar,[0.1 0.1]);
    %linestyle k
end
packfig(2,2)


%[q1b,q2b,q3b,q4b,qb]=divgeom(x(2)-x(1),kappabar,lambdabar,thetabar,5,'arakawa');


figure
subplot(2,1,1),jpcolor(x,y,q1+q2+q3+q4),axis equal,axis tight
subplot(2,1,2),jpcolor(x,y,q),axis equal,axis tight
for i=1:2
    subplot(2,1,i),caxis([-.2 .2]/16)
    hold on,contour(x,y,lambdabar,[0.1 0.1]);
    linestyle k
end
packfig(2,1)

function[]=divgeom_figure2

filename='/Users/lilly/Data/early/QGBetaPlaneTurbulenceFloats_experiment_04.nc';
ii=425:675;
[u,v]=vzeros(length(ii),1024,1001);
for k=1:1001
    k
    [utemp,vtemp]=FieldsFromTurbulenceFile(filename,k, 'u','v');
    u(:,:,k)=utemp(ii,:);
    v(:,:,k)=vtemp(ii,:);
end

u2=u;
v2=v;

for k=1:251
    k
    u2(k,:,:)=anatrans(u(k,:,:),3,'periodic');
    v2(k,:,:)=anatrans(v(k,:,:),3,'periodic');
end



