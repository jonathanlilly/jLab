function[xc,yc]=closedcurves(varargin)
%CLOSEDCURVES  Locate and interpolate closed curves in a possibly periodic domain.
%   _______________________________________________________________________
%   
%   *|* closedcurves.png --- Figure of CLOSEDCURVES applied to a QG model. 
%   Type 'jhelp closedcurves' to view this image. *|*
%   ______________________________________________________________________
%
%   [XC,YC]=CLOSEDCURVES(F,FO) returns closed curves of the matrix F at 
%   level FO, i.e. curves where F=FO, as the cell arrays XC and YC.  
%
%   If F has N contiguous patches at the contour level FO, then XC and YC
%   will be length N cell arrays.  
%
%   After rounding, the values of XC will be indices into the *columns* 
%   of F, while YC will be indices into the *rows*. 
%
%   [XC,YC]=CLOSEDCURVES(X,Y,F,FO) alternately specifies X and Y axes to go
%   with F.  X is associated with the *columns* of F and Y with its *rows*.
%   The lengths of X and Y must match SIZE(F,2) and SIZE(F,1) respectively.
%
%   Use CELLPLOT(XC,YC) to plot the curves. 
%
%   CLOSEDCURVES works by calling Matlab's CONTOUR routine, re-arranging
%   the contour matrix output, and throwing away non-closed contours. 
%
%   CLOSEDCURVES has options for considering a periodic domain, and for
%   interpolating the curves at a specified resolution, as described below.
%
%   The above figure shows closed curves of the Okubo-Weiss parameter P,
%   specifically the contour level P = -4 times its own standard deviation,
%   in a numerical simulation of quasi-geostrophic (QG) turbulence.
%   __________________________________________________________________
%
%   Periodic extension
%  
%   CLOSEDCURVES can look for curves within a singly or doubly periodic
%   version of the input field F.
%
%   CLOSEDCURVES(...,'periodic') will make F doubly periodic before looking
%   for closed curves.  
%
%   The X and Y axes will be extended based on their regular spacing, and 
%   these extended values will be returned in XC and YC for curves that 
%   fall into the periodically extended domain. 
%
%   CLOSEDCURVES(...,'periodic',N) extends F by N points in all directions.
%
%   CLOSEDCURVES(...,'periodic',N,M) extends F by M points at both the left
%   and the right (rows), and N points at the top and the bottom (columns).
%
%   To apply the periodic extension in only one of the dimesions, use
%   CLOSEDCURVES(...,'periodic',N,0) or CLOSEDCURVES(...,'periodic',0,M).
%
%   The above figure illustrates the difference between the periodized and
%   non-periodized algorithm.  In periodic domain such as this one, closed
%   curves at the region boundary will be missed wihtout periodization.
%
%   See PERIODIZE for details.
%   __________________________________________________________________
%
%   Interpolation
%
%   CLOSEDCURVES can also interpolate the curves to a specified number of
%   points along each curve, in order to increase the spatial resolution.  
%
%   [XC,YC]=CLOSEDCURVES(...,'interpolate',NPOINTS) will spline-interpolate
%   each curve to be length NPOINTS. 
%
%   In this case, XC and YC will be *matrices* instead of cell arrays.  The
%   number of rows of XC and YC is then NPOINTS, with each column 
%   corresponding to a separate curve.
%
%   Choosing a high value of NPOINTS, say NPOINTS=256, will minimize errors
%   in integral calculations based on these curves, such as those carried
%   out by CURVEMOMENTS.
%   __________________________________________________________________
%
%   See also CURVEMOMENTS, CURVEINTERP, PERIODIZE.  
%
%   'closedcurves --t' runs a test.
%   'closedcurves --f' generates the above figure.
%
%   Usage: [xc,yc]=closedcurves(f,fo);
%          [xc,yc]=closedcurves(x,y,f,fo);
%          [xc,yc]=closedcurves(x,y,f,fo,'interpolate',256);
%          [xc,yc]=closedcurves(x,y,f,fo,'periodic',M);
%          [xc,yc]=closedcurves(x,y,f,fo,'periodic',N,M);
%          [xc,yc]=closedcurves(x,y,f,fo,'periodic',N,M,'interpolate',256);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2015 J.M. Lilly --- type 'help jlab_license' for details


if ischar(varargin{1})
    if strcmpi(varargin{1}, '--t')
        closedcurves_test,return
    elseif strcmpi(varargin{1}, '--f')
        type makefigs_closedcurves
        makefigs_closedcurves;
        return
    end
end


M=[];
N=[];
str='noperiodize';
interpstr='nointerp';
na=nargin;

for i=1:3
    if ischar(varargin{end})
        str=varargin{end};
        varargin=varargin(1:end-1);
        na=na-1;
    elseif ischar(varargin{end-1})
        if strcmpi(varargin{end-1}(1:3),'int')
            interpstr=varargin{end-1};
            Mnew=varargin{end};
        else
            str=varargin{end-1};
            N=varargin{end};
            M=N;
        end
        varargin=varargin(1:end-2);
        na=na-2;
    elseif na>2
        if ischar(varargin{end-2})
            str=varargin{end-2};
            N=varargin{end-1};
            M=varargin{end};
            varargin=varargin(1:end-3);
            na=na-3;
        end
    end
end
qo=varargin{end};
varargin=varargin(1:end-1);

na=length(varargin);
if na==3;
    x=varargin{1};
    y=varargin{2};
    q=varargin{3};
else
    x=[];
    y=[];
    q=varargin{1};
end

if strcmpi(str(1:3),'per')
    if isempty(M)
        M=size(q,1);
    end
    if isempty(N)
        N=size(q,2);
    end
    xin=x;
    yin=y;
    [x,y,q]=periodize(N,M,x,y,q);
end
figure('visible','off')

xs=0;
if isempty(x)
    [c,h]=contour(q,[qo qo],'b');
else
    %Note!!! Matlab's contour has a major bug and will work badly if
    %contour(x,y,z) has x containing zero.  To see this try
    %[c,h]=contour(x,y,q,[qo qo],'b');figure,plot(c(1,:),'.'),figure,plot(c(2,:),'.')
    %A value of 0 could *either* mean 'the end of a contour', or else 
    %'this contour crosses x=0'!  To get around this, we have to shift
    %x to not include zero.  
    if minmin(x)<0
        xs=1-minmin(x);
    end
    [c,h]=contour(x+xs,y,q,[qo qo],'b');
end
close


if isempty(c)
    xc=[];
    yc=[];
else
    index=find(c(1,:)==qo);
    index(end+1)=size(c,2)+1;
    for i=1:length(index)-1
        xc{i,1}=c(1,index(i)+1:index(i+1)-1)'-xs;
        yc{i,1}=c(2,index(i)+1:index(i+1)-1)';
    end
end


for i=1:length(xc)
    %Prune those that are not closed
    %Possibly add a tolerance here?  
    bprune=~(xc{i}(1)==xc{i}(end))||~(yc{i}(1)==yc{i}(end));
    
    if strcmpi(str(1:3),'per')
        %Prune those outside principal domain
        if ~bprune
            [xo,yo]=curvemoments(xc{i},yc{i});
            bprune=(xo>max(xin))||(xo<min(xin))||(yo>max(yin))||(yo<min(yin));
        end
    end
    if bprune
        xc{i}=[];
        yc{i}=[];
    end
end

lenx=cellength(xc);
leny=cellength(yc);
bool=(lenx>=3)&(leny>=3);
xc=xc(bool);
yc=yc(bool);

if strcmpi(interpstr(1:3),'int')
    [xcnew,ycnew]=vzeros(Mnew,length(xc));
    for i=1:length(xc)
        t=(1:length(xc{i}))';
        tnew=linspace(1,length(xc{i}),Mnew)';
        xcnew(:,i)=interp1(t,xc{i},tnew,'spline');
        ycnew(:,i)=interp1(t,yc{i},tnew,'spline');
    end
    xc=xcnew;
    yc=ycnew;
end

function[]=closedcurves_test
load qgsnapshot

dx=qgsnapshot.x(2)-qgsnapshot.x(1);
[cv,zeta,N,S,P]=psi2fields(dx,qgsnapshot.psi);
P=frac(P,std(P(:)));

[xc,yc]=closedcurves(qgsnapshot.x,qgsnapshot.y,P,-4);
[ic,jc]=closedcurves(P,-4);

cell2col(xc,yc,ic,jc);
bool=isfinite(xc.*yc.*ic.*jc);
vindex(xc,yc,ic,jc,bool,1);

clear bool
bool(1)=aresame(xc',qgsnapshot.x(round(ic)),dx);
bool(2)=aresame(yc,qgsnapshot.y(round(jc)),dx);

reporttest('CLOSEDCURVES matches for index vs. input axes variations',allall(bool))
%[xc,yc]=closedcurves(qgsnapshot.x,qgsnapshot.y,zeta,0);


% zp=curveinterp(qgsnapshot.x,qgsnapshot.y,cv,xp,yp);
% [xo,yo,L,R,D,a,b,theta]=curvemoments(xp,yp);
% [vort,div,mom,ke,momstd]=curvemoments(xp,yp,zp);
% 
% [xpi,ypi]=closedcurves(qgsnapshot.x,qgsnapshot.y,P,-2,'periodic',100,200,'interpolate',256);
% zpi=curveinterp(qgsnapshot.x,qgsnapshot.y,cv,xpi,ypi);
% [xoi,yoi,Li,Ri,Di,ai,bi,thetai]=curvemoments(xpi,ypi);
% [vorti,divi,momi,kei,momstdi]=curvemoments(xpi,ypi,zpi);
% 
% figure,
% subplot(1,2,1),
% plot(R,div,'ro'),hold on,plot(Ri,divi,'b+'),ylim([-5 5]*1e-6),xlim([0 2.49])
% title('Divergence with (blue) and without (red) interpolation')
% xlabel('Radius R (km)'),ylabel('Spatially-averaged divergence (1/s)')
% subplot(1,2,2),
% plot(R,vort,'ro'),hold on,plot(Ri,vorti,'b+'),ylim([-5 5]*1e-6),xlim([0 2.49])
% title('Vorticity with (blue) and without (red) interpolation')
% xlabel('Radius R (km)'),ylabel('Spatially-averaged vorticity (1/s)')
% packfig(1,2,'columns')
% fontsize 18 14 14 14
