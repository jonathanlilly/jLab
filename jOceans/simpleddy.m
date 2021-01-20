function[v,psi,zeta,psif]=simpleddy(varargin)
%SIMPLEDDY  Streamfunction, velocity, and vorticity for various eddy profiles. 
%
%   [V,PSI,ZETA]=SIMPLEDDY(R,RO,VO) returns the streamfunction PSI,
%   azimuthal velocity V, and relative vorticity ZETA at radial locations R
%   for the Rankine vortex.
%
%   The eddy core radius is RO, in km, and at this radius the eddy 
%   currents achieve their maximum value of |VO| cm/s.  The sign of VO
%   gives the rotation sense of the eddy, positive for counterclockwise.
%
%   PSI has units of km^2/s, V has units of cm/s, and ZETA of 1/s. 
%
%   RO and VO should either be scalars or an arrays of the same size as R. 
%   The output arrays will have size LENGTH(R) x N.
%   _______________________________________________________________________
%
%   Eddy models
%
%   SIMPLEDDY can also return profiles associated with various other types
%   of model eddies.  The following forms are available:
%
%         SIMPLEDDY(R,RO,VO,'rankine')        -- Rankine vortex   [default]
%         SIMPLEDDY(R,RO,VO,'gaussian')       -- Gaussian vortex
%
%   RO, VO, and the shape parameter ALPHA may then each be either scalars
%   or arrays of the same length. 
%   _______________________________________________________________________
% 
%   Eddy 'slices'
%
%   SIMPLEDDY can also be used to return an eddy 'slice', that is, the
%   eddy streamfunction, currents, and vorticity observed along a curve.
%  
%   CV=SIMPLEDDY(CX,RO,VO,...) where CX=X+iY is an array of locations on 
%   the complex plane, with units of kilometers, returns the complex-valued
%   currents CV=U+iV due to an eddy located at the origin.
%
%   [CV,PSI,ZETA]=SIMPLEDDY(CX,RO,VO,...) also returns the streamfunction
%   and vorticity observed along this curve.  
%   _______________________________________________________________________
%
%   'simpleddy --f' makes some sample figures.
%
%   Usage: [v,psi,zeta]=simpleddy(r,RO,VO);
%          [v,psi,zeta]=simpleddy(r,RO,VO,'gaussian');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015--2020 J.M. Lilly --- type 'help jlab_license' for details
 

%         SIMPLEDDY(R,RO,VO,'matern',ALPHA)   -- Matern vortex
%
%   The Matern family is described by an additional shape parameter, 
%   termed ALPHA.  This parameter must be ALPHA>1/2, with the Matern family
%   becoming increasingly similar to a Gaussian as ALPHA increases.

%   See also MATERNFUN, EDDYCURVE, EDDYFIT, EDDYGUESS.



%   XXX what if X has more than one column.  
%          [v,psi,zeta]=simpleddy(r,RO,VO,'matern',ALPHA);


if strcmpi(varargin{1}, '--f')
    rankineddy_fig
    gaussian_fig
    %materneddy_figures
    return
end

if strcmpi(varargin{1}, '--t')
    simpleddy_test;
    return
end

x=varargin{1};
if size(x,1)==1
    disp('Transposing the first argument to SIMPLEDDY to make it a column vector.')
    x=conj(x');
end

ro=varargin{2};
vo=varargin{3};

eddystr='gaussian';
alpha=[];
        
if nargin>3
    if ~ischar(varargin{4})
        alpha=varargin{4};
        if alpha==0
             eddystr='rankine';
             alpha=1;
        elseif alpha==inf
            eddystr='gaussian';
        elseif alpha>1/2
            eddystr='matern';
        else
            error(['Sorry, the value ALPHA=' num2str(alpha) ' is not supported.'])
        end
    else
        if ischar(varargin{end})
            eddystr=varargin{end};
            varargin=varargin(1:end-1);
        elseif ischar(varargin{end-1})
            eddystr=varargin{end-1};
            alpha=varargin{end};
            varargin=varargin(1:end-2);
        end
    end
end

vo=vo/100/1000;
if strcmpi(eddystr(1:3),'mat')
    [v,psi,zeta]=materneddy(alpha,ro,vo,abs(x));
elseif strcmpi(eddystr(1:3),'gau')
    [v,psi,zeta]=gaussianeddy(ro,vo,abs(x));
elseif strcmpi(eddystr(1:3),'ran')
    [v,psi,zeta]=rankineddy(alpha,ro,vo,abs(x));
end


if size(x,2)==1
    x=vrep(x,size(v,2),2);
end
if ~isreal(x)||anyany(x<0)
    v=sqrt(-1)*frac(x,abs(x)).*v;
end

v(~isfinite(v))=0;
psi(~isfinite(psi))=0;
zeta(~isfinite(zeta))=0;



function[v,psi,zeta]=rankineddy(alpha,ro,vo,r)
%RANKINEDDY Velocity and streamfunction for a generalized Rankine vortex.
%
%   RANKINEDDY computes the velocity and streamfunction observed at
%   the origin due to a Rankine vortex located at a specified point,
%   as described in Lilly and Rhines (2002).
%
%   [v,psi]=rankineddy(eta,ro,vo);
%
%     Input
%	eta: complex-valued position of eddy (km)
%	ro:  size of eddy (km)
%	vo:  velocity at edge, negative for anticyclone (cm/s) 
%
%     Output
%	v,psi:	azimuthal velocity and streamfunction at origin
%
%   'rankineddy --f' generates a sample figure
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2013 J.M. Lilly --- type 'help jlab_license' for details    

if isempty(alpha)||~isfinite(alpha)
    alpha=1;
end

iO=find(r>ro);
v=vo.*r./ro;
vO=vo.*(ro.^alpha)./(r.^alpha);
if ~isempty(iO)
 %size(r),size(v), size(ro),size(iO), size(vo)
    v(iO)=vO(iO);
end
v=v*100*1000;

psi=(1/2).*(r./ro).^2;
if alpha==1
    psiO=log(r./ro)+1/2;
else
    psiO=1/2+frac(1,alpha-1)-frac(1,alpha-1).*frac(ro.^(alpha-1),r.^(alpha-1));
end

if ~isempty(iO)
    psi(iO)=psiO(iO);
end
psi=psi.*vo.*ro;

zeta=2*vo./ro+0*r;
zetaO=frac(vo.*ro,r.^(alpha+1)).*(1-alpha);
if ~isempty(iO)
    zeta(iO)=zetaO(iO);
end

function[v,psi,zeta]=gaussianeddy(ro,vo,r)
%GAUSSIANEDDY Velocity and streamfunction for a Gaussian vortex.
%
%   GAUSSIANEDDY computes the velocity and streamfunction observed at
%   the origin due to an eddy having a Gaussian streamfunction profile
%   located at a specified point.
%
%   [v,psi]=gaussianeddy(eta,ro,vo);
%
%     Input
%	eta: complex-valued position of eddy (km)
%	ro:  size of eddy (km)
%	vo:  velocity at edge, negative for anticyclone (cm/s) 
%
%     Output
%	v,psi:	azimuthal velocity and streamfunction at origin
%
%   'gaussianeddy --f' generates a sample figure
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2012 J.M. Lilly --- type 'help jlab_license' for details    

psi=-1.*vo.*ro.*exp(1/2.*(1-r.^2./ro.^2));
v=100*1000*r.*(vo./ro).*exp(1/2.*(1-r.^2./ro.^2));
zeta=(vo./ro).*(2-r.^2./ro.^2).*exp(1/2.*(1-r.^2./ro.^2));



function[v,psi,zeta]=materneddy(alpha,R,V,r)
%MATERNEDDY
%
%   MATERNEDDY
%
%   Usage: []=materneddy();
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2012 J.M. Lilly --- type 'help jlab_license' for details
 
%vsize(R,V,alpha,r)
%vsize(R,V,alpha,r)
%Halpha=gamma(alpha).*(2.^(alpha-1));  
%Don't need this because I only need the ratio H/Halpha

%figure,plot(abs(r(:)))
alpha
R
V
maternfun(alpha-1,maternrad(alpha))
maternrad(alpha)
alpha=alpha(:);
R

%vcolon(R,V,alpha);
R=R(:)';
V=V(:)';
alpha=alpha(:)';

c=-frac(maternrad(alpha).^2,R).*maternfun(alpha-1,maternrad(alpha))';
%vsize(V,c)
Hnorm=frac(V,c);
%xx=maternrad(alpha);
%figure,plot(xx(:));
Rnorm=R./maternrad(alpha);
%vsize(Rnorm,alpha,r)

size(Rnorm)
size(Hnorm)


if ~aresame(size(alpha),size(r))
    if size(r,2)==1&&(length(alpha)~=1)
        r=vrep(r,length(alpha),2); 
    end
    if size(r,1)~=1&&(length(alpha)~=1)
        Rnorm=vrep(Rnorm,size(r,1),1);
        Hnorm=vrep(Hnorm,size(r,1),1);   
    end
end

%figure,plot(abs(Rnorm(:)),'.')
%figure,plot(abs(Hnorm(:)))
%figure,plot(abs(alpha(:)))
vsize(r,Rnorm,alpha)

v=-100*1000*Hnorm.*r.*(1./Rnorm.^2).*maternfun(alpha-1,r./Rnorm);

if nargout >=2
    psi=Hnorm.*maternfun(alpha,r./Rnorm);
end
if nargout ==3
    zeta= Hnorm.*(1./Rnorm.^2).*((r./Rnorm).^2.*maternfun(alpha-2,r./Rnorm)-2*maternfun(alpha-1,r./Rnorm));
end


function[]=materneddy_figures
dr=0.1;
r=[0:dr:500]';
%alpha=[1/2:1/4:10];
%alpha=logspace(log10(1.05),log10(25),10);
alpha=logspace(log10(3/4),log10(50),10);
alpha=([1.2:.2:2 3:2:11]);
[v,psi,zeta]=simpleddy(r,10,15,'matern',alpha);
[vo,psio,zetao]=simpleddy(r,10,15,'gaussian');
%[vo,psio,zetao]=simpleddy(r,10,15,'rankine');

psi=[psio psi];
v=[vo v];
zeta=[zetao zeta];

c=frac(maternfun(alpha-1,maternrad(alpha)),maternfun(alpha-1,0));
c=[1 c.*sqrt(exp(1))];
for i=1:length(c)
    psi(:,i)=psi(:,i).*c(i);
    v(:,i)=v(:,i).*c(i);
    zeta(:,i)=zeta(:,i).*c(i);
end

psitwosides=[flipud(psi);psi(2:end,:)];
Psi=abs(fft(psitwosides));
%[mu,sigma,skew,kurt]=pdfprops(1:length(Psi),psitwosides(:,[2:end 1]));
%[muf,sigmaf,skewf,kurtf]=pdfprops(1:length(Psi),fftshift(Psi(:,[2:end 1])));

f=fourier(length(Psi));
Psi=Psi(1:length(r),:);

figure
subplot(2,2,1),h=plot(r,psi);vlines(10,'k:'), xlim([0 50]),
h1=h(1);h=h(2:end);title('Streamfunction'),xlabel('Distance')
linestyle -h h1 -4K
linestyle -h h default
subplot(2,2,2),h=plot(r,v);vlines(10,'k:'),ylim([0 15.5]), xlim([0 50])
h1=h(1);h=h(2:end);title('Velocity'),xlabel('Distance')
linestyle -h h1 -4K 
linestyle -h h default
subplot(2,2,3),h=plot(r,zeta);vlines(10,'k:'),hlines(0,'k:'),ylim([-1 5.1]/1e5), xlim([0 50])
h1=h(1);h=h(2:end);title('Vorticity'),xlabel('Distance')
linestyle -h h1 -4K
linestyle -h h default
subplot(2,2,4),h=plot(f./dr,Psi);xlog,ylog,ylim(10.^[-6 1]),vlines(1./10,'k:'),xlim([minmin(f./dr) maxmax(f./dr)])
h1=h(1);h=h(2:end);title('Spectrum'),xlabel('Radian Wavenumber')
linestyle -h h1 -4K
linestyle -h h default
%letterlabels(2)

fontsize 14 12 12 12
orient landscape
%print -depsc materneddies_core.eps


alpha=logspace(log10(3/4),log10(25),100);
figure,
plot(alpha,frac(maternfun(alpha-1,maternrad(alpha)),maternfun(alpha-1,0)))
linestyle -2k
hlines(1./sqrt(exp(1)),'k:')
title('Rossby Number Ratio'),xlabel('Shape Parameter \alpha'),xlim([0 25])
fontsize 14 12 12 12
orient portrait
%print -depsc rossbyratio.eps

%A la McWilliams
figure
h=plot((r./10).^2,abs(v));vlines(10,'k:'),hlines(0,'k:'),xlim([0 25]),ylog
h1=h(1);h=h(2:end);title('Velocity'),xlabel('Distance Squared'),ylim([.1 15])
orient portrait


function[]=simpleddy_test

x=[-100:.1:100]';
y=[0:2:16];
clear cx
for i=1:length(y)
    cx(:,i)=x+sqrt(-1)*y(i);
end

Ro=10;
Vo=15;

bool(1)=aresame(simpleddy(cx,Ro,Vo,'matern',.75),simpleddy(cx,Ro,Vo,.75));
bool(2)=aresame(simpleddy(cx,Ro,Vo,'gaussian'),simpleddy(cx,Ro,Vo,inf));
bool(3)=aresame(simpleddy(cx,Ro,Vo,'rankine'),simpleddy(cx,Ro,Vo,0));
reporttest('SIMPLEDDY alternate input format matches string format',all(bool))

function[]=eddyslice_figures

x=[-100:.1:100]';
y=-[0:2:16];
clear cx
for i=1:length(y)
    cx(:,i)=x+sqrt(-1)*y(i);
end

Ro=10;
Vo=15;


%[cv,psi,zeta]=simpleddy(cx,Ro,Vo,'gaussian');
%[cv,psi,zeta]=simpleddy(cx,Ro,Vo,inf);

[cv,psi,zeta]=simpleddy(cx,Ro,Vo,'matern',.75);

%[cv,psi,zeta]=simpleddy(cx,Ro,Vo,'matern',1/2);
%[cv,psi,zeta]=simpleddy(cx,Ro,Vo,'matern',1);
%[cv,psi,zeta]=simpleddy(cx,Ro,Vo,'matern',2);
%[cv,psi,zeta]=simpleddy(cx,Ro,Vo);

figure,
subplot(2,2,1),plot(cx),hold on,axis([-1 1 -1 1]*Ro*2),axis equal
h=plot(Ro*rot(0:.01:2*pi+.1));
linestyle -h h 2G
subplot(2,2,2),hodograph(cv)
subplot(2,2,3),plot(x,real(cv)),vlines([-Ro Ro],'k:');hlines(0,'k:');
title('Along-Stream Currents')
subplot(2,2,4),plot(x,imag(cv)),vlines([-Ro Ro],'k:');hlines(0,'k:');
title('Cross-Stream Currents')

% 
% x=[-20:.1:20]';
% [xg,yg]=meshgrid(x,x);
% [cv,psi,zeta]=simpleddy(xg+sqrt(-1)*yg,Ro,Vo,'matern',1);

function[]=manyeddy_fig
ro=10;
vo=-10;
dx=0.1
eta=osum((-500:dx:500)', sqrt(-1)*(0:1:50)');
[taper,lambda]=sleptap(length(eta),4); 


N=size(taper,1);
[v,psi]=simpleddy(eta,ro,vo,'matern',1);
[f,s]=mspec(.1,psi,ones(N,1)./sqrt(N)); 
[v,psi]=simpleddy(eta,ro,vo,'matern',2);
[f(:,:,2),s(:,:,2)]=mspec(.1,psi,ones(N,1)./sqrt(N)); 
[v,psi]=simpleddy(eta,ro,vo,'matern',4);
[f(:,:,3),s(:,:,3)]=mspec(.1,psi,ones(N,1)./sqrt(N)); 
[v,psi]=simpleddy(eta,ro,vo,'matern',8);
[f(:,:,4),s(:,:,4)]=mspec(.1,psi,ones(N,1)./sqrt(N)); 
[v,psi]=simpleddy(eta,ro,vo,'matern',16);
[f(:,:,5),s(:,:,5)]=mspec(.1,psi,ones(N,1)./sqrt(N)); 
[v,psi]=simpleddy(eta,ro,vo,'gaussian');
[f(:,:,6),s(:,:,6)]=mspec(.1,psi,ones(N,1)./sqrt(N)); 


figure
for i=1:6
subplot(3,2,i)
plot(f(:,:,i),s(:,:,i)),xlog,ylog,hold on,linestyle D,
xlim([2*pi/1000 2*pi/1]),ylim([10^(-30) 10^(-4)]),
h=plot(f(:,:,i),vmean(s(:,:,i),2));linestyle -h h 2k
vlines(2*pi/10,'g')
end

packfig(3,2)


function[]=gaussian_fig
ro=10;
vo=-10;
eta=-sqrt(-1)*5+(-50:.1:50)';

%[v,psi]=simpleddy(eta,ro,vo,inf);
[v,psi]=simpleddy(eta,ro,vo,'gaussian');

figure,
subplot(121),
uvplot(real(eta),v),
title('Anticyclone sliced north of center, halfway to edge')
subplot(122)
polar(angle(v),abs(v))
title('Hodograph of Gaussian eddy')

function[]=rankineddy_fig
ro=10;
vo=-10;
eta=-sqrt(-1)*5+(-50:.1:50)';
[v,psi]=simpleddy(eta,ro,vo,'rankine',2);

figure,
subplot(121),
uvplot(real(eta),v),
title('Anticyclone sliced north of center, halfway to edge')
subplot(122)
polar(angle(v),abs(v))
title('Hodograph of Rankine eddy')
