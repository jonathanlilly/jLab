function[ci]=curveinterp(varargin)
%CURVEINTERP  Interpolate a field or its gradient onto a set of curves.
%
%   CURVEINTERP is a low-level function called by CLOSEDCURVES and
%   CURVEMOMENTS.
%
%   FC=CURVEINTERP(X,Y,F,XC,YC) interpolates the field in matrix F, with 
%   axes X and Y, onto the curves specified by XC and YC.  
%
%   X is associated with the *columns* of F and Y with its *rows*.  The 
%   lengths of X and Y must match SIZE(F,2) and SIZE(F,1) respectively.
%
%   Unlike INTERP2, CURVEINTERP permits XC and YC to be cell arrays of 
%   column vectors, in which case FC will be a cell array of the same size.
%   This is useful for working with CLOSEDCURVES and CURVEMOMENTS.
%
%   CURVEINTERP uses a fast loopless algorithm that is an order of 
%   magnitude faster than looping over the individual curves. 
%
%   CURVEINTERP works by calling INTERP2 with the 'linear' algorithm.
%   __________________________________________________________________
%
%   Gradient interpolation
%
%   GRADFC=CURVEINTERP(X,Y,F,XC,YC,'gradient') alternately interpolates not
%   F, but instead its *gradient*, and returns the result in the complex-
%   valued quantity GRADFC=dF/dX+sqrt(-1)*dF/dY.
%  
%   This can be very useful when F is large, because it is not necessary
%   to compute the gradient of F everywhere.  It is only computed locally
%   in near the curves, which can lead to a considerable speed improvment.
%   __________________________________________________________________
%
%   See also CLOSEDCURVES, CURVEMOMENTS.
%
%   'curveinterp --t' runs a test.
%
%   Usage: fc=curveinterp(x,y,f,xc,yc);
%          gradfc=curveinterp(x,y,f,xc,yc,'gradient');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2015 J.M. Lilly --- type 'help jlab_license' for details
 
 
if strcmpi(varargin{1}, '--t')
    curveinterp_test,return
end

str='field';
algstr='fast';
for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'fie')||strcmpi(varargin{end}(1:3),'gra')
             str=varargin{end};
        else
             algstr=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end

x=varargin{1};
y=varargin{2};
varargin=varargin(3:end);

if aresame(size(x),size(varargin{2}))||aresame(length(x),size(varargin{2},2))
    z=varargin{1};
    c=varargin{2};
    varargin=varargin(3:end);
else
    c=varargin{1};
    z=[];
    varargin=varargin(2:end);
end

xc=varargin{1};
yc=varargin{2};
varargin=varargin(3:end);

if ~isempty(z)
    zc=varargin{1};
    varargin=varargin(2:end);
else
    zc=[];
end 

if ~iscell(xc)
    bcellinput=false;
    sizexc=size(xc);
    xtemp=xc;clear xc
    ytemp=yc;clear yc
    xc{1}=xtemp(:);
    yc{1}=ytemp(:);
    if ~isempty(zc)
        ztemp=zc;clear zc
        zc{1}=ztemp(:);
    end
else
    bcellinput=true;
end

interpstr='linear';
%interpstr='spline';  %Seems to make little difference

if strcmpi(str(1:3),'gra')
    dx=frac(1,2)*(x(2)-x(1));  %Experimentation shows this is a good choice
                               %Better than dx and much better than dx/4
    if isempty(z)
        cxp=curveinterp(x,y,c,celladd(dx,xc),yc,algstr);
        cxm=curveinterp(x,y,c,celladd(-dx,xc),yc,algstr);
        cyp=curveinterp(x,y,c,xc,celladd(dx,yc),algstr);
        cym=curveinterp(x,y,c,xc,celladd(-dx,yc),algstr);        
    else
        cxp=curveinterp(x,y,z,c,celladd(dx,xc),yc,zc,algstr);
        cxm=curveinterp(x,y,z,c,celladd(-dx,xc),yc,zc,algstr);
        cyp=curveinterp(x,y,z,c,xc,celladd(dx,yc),zc,algstr);
        cym=curveinterp(x,y,z,c,xc,celladd(-dx,yc),zc,algstr);
    end   
    cir=celladd(cellmult(frac(1,2*dx),cxp),cellmult(-frac(1,2*dx),cxm));%figure,cellplot(cir)
    cii=celladd(cellmult(frac(1,2*dx),cyp),cellmult(-frac(1,2*dx),cym));%figure,cellplot(cii)
    ci=celladd(cir,cellmult(sqrt(-1),cii));
else
    ci=cellmult(nan,xc);

    if strcmpi(algstr(1:3),'loo')
        for i=1:length(xc);
            if isempty(z)
                %vsize(x,y,c,xc{i},yc{i})
                ci{i}=interp2(x,y,c,xc{i},yc{i},interpstr);
            else
                ci{i}=interp3(x,y,z,c,xc{i},yc{i},zc{i},interpstr);
            end
        end
    elseif strcmpi(algstr(1:3),'fas')
        cell2col(xc,yc);
        ci=nan*xc;
        %For boundary condition treatment, fix later
        %[xc,yc]=curveinterp_boundaries(x,y,xc,yc);
        xc(xc<minmin(x))=minmin(x);  
        xc(xc>maxmax(x))=maxmax(x);  
        yc(yc<minmin(y))=minmin(y);  
        yc(yc>maxmax(y))=maxmax(y);  
        
        if isempty(z)
            bool=isfinite(xc)&isfinite(yc);
            %vsize(x,y,c,xc(bool),yc(bool))
            ci(bool)=interp2(x,y,c,xc(bool),yc(bool),interpstr);
        else
            cell2col(zc);
            bool=isfinite(xc)&isfinite(yc)&isfinite(zc);
            ci(bool)=interp3(x,y,c,xc(bool),yc(bool),zc(bool),interpstr);
        end
        col2cell(ci);
    end
end

if ~bcellinput
    ci=ci{1};
    ci=reshape(ci,sizexc);
    %ci=reshape(ci,sizexc(2),sizexc(1));
    %ci=conj(ci');
end


function[xc,yc]=curveinterp_boundaries(x,y,xc,yc)
%any([minmin(xc)<minmin(x) maxmax(xc)>maxmax(x) minmin(yc)<minmin(y) maxmax(yc)>maxmax(y)])

%figure,plot(xc,yc,'b.')
ddx=x(2)-x(1);
dx=minmin(x)-xc;
xc(dx>0)=maxmax(x)-dx(dx>0);

dx=maxmax(x)-xc;
xc(dx<0)=minmin(x)-dx(dx<0);

dy=minmin(y)-yc;
yc(dy>0)=maxmax(y)-dy(dy>0);

dy=maxmax(y)-yc;
yc(dy<0)=minmin(y)-dy(dy<0);

%hold on,plot(xc,yc,'ro')

%any([minmin(xc)<minmin(x) maxmax(xc)>maxmax(x) minmin(yc)<minmin(y) maxmax(yc)>maxmax(y)])

function[]=curveinterp_test

load qgsnapshot
dx=qgsnapshot.x(2)-qgsnapshot.x(1);
tic;[cv,zeta,N,S,P]=psi2fields(dx,qgsnapshot.psi);t0=toc;
P=frac(P,std(P(:)));
[xc,yc]=closedcurves(qgsnapshot.x,qgsnapshot.y,P,-5);

tic;cvc=curveinterp(qgsnapshot.x,qgsnapshot.y,cv,xc,yc,'loop');t1=toc;
tic;cvc2=curveinterp(qgsnapshot.x,qgsnapshot.y,cv,xc,yc,'fast');t2=toc;

reporttest('CURVEINTERP fast and slow algorithms match',aresame(cell2col(cvc),cell2col(cvc2)))
disp(['CURVEINTERP fast algorithm was ' int2str(t1./t2) ' times faster than loops.'])

cvc3=vempty;
for i=1:length(xc);
    tempu=interp2(qgsnapshot.x,qgsnapshot.y,real(cv),xc{i},yc{i},'linear');
    tempv=interp2(qgsnapshot.x,qgsnapshot.y,imag(cv),xc{i},yc{i},'linear');
    cvc3{i,1}=tempu+sqrt(-1)*tempv;
end
reporttest('CURVEINTERP matches manual interpolation',aresame(cell2col(cvc),cell2col(cvc3)))

tic;gradpsi=curveinterp(qgsnapshot.x,qgsnapshot.y,qgsnapshot.psi,xc,yc,'gradient');t3=toc;
x1=cell2col(cvc2);
x2=sqrt(-1)*cell2col(gradpsi)*100/1000;  %Convert psi from m^2 to cm * km
err=frac(vsum(abs(x1-x2).^2,1),vsum(abs(x1).^2,1));

reporttest('CURVEINTERP gradient interpolation',err<4e-5)
disp(['CURVEINTERP gradient interpolation was ' int2str((t0+t2)./t3) ' times faster than direct calculation.'])

