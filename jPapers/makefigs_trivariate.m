function[varargout]=makefigs_trivariate(str)
%MAKEFIGS_TRIVARIATE  Make figures for Lilly (2011).
% 
%   This function makes all figures for the paper 
%
%   Lilly, J. M. (2011). Modulated oscillations in three dimensions. IEEE
%      Transactions on Signal Processing, 59 (12), 5930--5943.
%
%   Usage: makefigs_trivariate
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2023 J.M. Lilly --- type 'help jlab_license' for details

if nargin==0
    str='noprint';
end

%/************************************************************************
%Trivariate Example
 
%t=(0:.1:8000)';
t=(0:800)';

betao=pi/6;
alphao=pi/6;
thetao=pi/4;

phi=vrep((t/1000*5)*2*pi,5,2);
kappa=2.5+0*phi;
lambda=0.4+0*phi;
theta=thetao+0*phi;
alpha=alphao+0*phi;
beta=betao+0*phi;

%Column one is amplitude modulation
kappa(:,1)=3*exp(2*0.393*(t/1000-1));

%Column 2 is nutation
beta(:,2)=beta(:,2)+phi(:,2)./14.45*lambda(1,2)*sqrt(2)*0.904;

%Column three is precession
phi(:,3)=phi(:,3)*0.941;
%phi(:,3)=phi(:,1);  %Not sure why I'm not using this one right now
theta(:,3)=pi/4+phi(:,3)./14.45*.96;

%Column 4 is external precession
lambda(:,4)=0;
phi(:,4)=phi(:,4)*1.062;%+(t-mean(t)).^2*.00000004;
alpha(:,4)=alpha(:,4)-phi(:,4)./14.45*2*lambda(1,3)*sqrt(2)*.851;

%Column five is deformation
rho=1/636*(t+2);
lambda(:,5)=sin(rho);
lambda(lambda>1)=1;

[x,y,z]=ellsig(kappa,lambda,theta,phi,alpha,beta,'real');
[a,b,c,d,e]=ellband(kappa,lambda,theta,phi,alpha,beta); 

om=vdiff(phi,1)+sqrt(1-lambda.^2).*(vdiff(theta,1)+cos(beta).*vdiff(alpha,1));
ombar=vmean(om,1,squared(kappa));
sig=sqrt(abs(a).^2+abs(b).^2+abs(d).^2+abs(e).^2+(om-vrep(ombar,length(om),1)).^2);
sigbar=sqrt(vmean(sig.^2,1,squared(kappa)));


%Tests
%aresame(ombar,0.0314*[1 1 1 1 1],1e-3)
%b1=a(:,1);b2=abs(e(:,2));b3=c(:,3);b4=sqrt(d(:,4).^2+abs(e(:,4)).^2);b5=abs(b(:,5));
%aresame(vmean([b1 b2 b3 b4 b5],1),0.786*[1 1 1 1 1]./1000,1e-5)


%plot([b1 b2 b3 b4 b5])
%aresame(vmean(a,1),0.0314*[1 1 1 1 1],1e-3)
%plot(sig)


textstr{1}='Amplitude Modulation';
textstr{2}='Nutation';
textstr{3}='Internal Precession';
textstr{4}='External Precession';
textstr{5}='Deformation';

%/****************
psi=sleptap(length(t),2);
%psi=ones(size(t))./length(t);
clear s
for i=1:5
    for j=1:3
        [f,s(:,j,i)]=mspec(x(:,i),psi);
        [f,s(:,j,i)]=mspec(y(:,i),psi); 
        [f,s(:,j,i)]=mspec(z(:,i),psi);   
    end
end
sbar=squeeze(vsum(s,2));
for i=1:size(sbar,2)
    sbar(:,i)=sbar(:,i)./vsum(sbar(:,i),1);
end
%\****************

% 
% clear s
% for i=1:5
%     for j=1:3
%         s(:,j,i)=abs(fft(detrend(x(:,i))));
%         s(:,j,i)=abs(fft(detrend(y(:,i))));
%         s(:,j,i)=abs(fft(detrend(z(:,i))));
%     end
% end
% sbar=squeeze(vsum(s,2));
% sbar=sbar(end/2+1:end,:);
% for i=1:size(sbar,2)
%     sbar(:,i)=sbar(:,i)./vsum(sbar(:,i),1);
% end
% 
[mu,sig]=pdfprops(vrep(f,5,2),sbar);

%Since the plot will read up-down not left-right
numplot=[1 4 2 5 3 6];

figure

for i=1:5
    subplot(3,2,i)
    scatter3(x(:,i),y(:,i),z(:,i),5+0*y(:,i),t,'filled'),colormap gray,hold on
    scatter3(x(:,i),y(:,i),-3.5+0*z(:,i),5+0*y(:,i),t,'filled'),colormap gray
    set(gca,'dataaspectratio',[ 1 1 1])
    axis([-1 1 -1 1 -1 1]*3.5)
    xtick([-4:2:4]),ytick([-4:2:4]),ztick([-4:2:4])
    set(gca,'box','on'),grid off
    hlines([-3.5 3.5],'k:'),vlines([-3.5 3.5],'k:') 
    set(gca,'xticklabel',[]),set(gca,'yticklabel',[]),set(gca,'zticklabel',[])
    view(-25,10)
    plot3([0 0],[0 0 ],[-3.5 3.5],'k--')
    plot3(0,0,0,'ko','markersize',5,'markerfacecolor','k')
    if i==1||i==3||i==5
        [xs,ys,zs]=vectmult(jmat3(alphao,3)*jmat3(betao,1),[-1 -1 1 1 -1]*3.2,[-1 1 1 -1 -1]*3.2,[0 0 0 0 0]);
        plot3(xs,ys,zs,'k')
    else
        sty{1}='K';sty{2}='G';sty{3}='C';
        index=[1 floor(length(t)/2) length(t)];
        for j=1:length(index)
            J=jmat3(alpha(index(j),i),3)*jmat3(beta(index(j),i),1);
            [xs,ys,zs]=vectmult(J,[-1 -1 1 1 -1]*3.2,[-1 1 1 -1 -1]*3.2,[0 0 0 0 0]);
            h=plot3(xs,ys,zs);linestyle(h,sty{j})
        end 
    end
    title(textstr{i}) 
    text(-3,3,-2.8,['(' char(96+numplot(i)) ')'])
    
    if ~verLessThan('matlab','8.4.0')
        set(gca,'BoxStyle','full')
    end
end


subplot(3,2,6),cla
plot(f,sbar,'k'),hold on,plot(f,sbar,'k.'),xlog,ylog
xlim([10^(-3) max(f)]),ylim([10^(-7) 1])
title('Average Spectra')
xlabel('Frequency (cycles/point)')
text(10^-2.8,10^-6.5,'(f)')
vlines(1/2/100,'k:')


orient tall
fontsize jpofigure 

if strcmpi(str,'print')
    print -depsc trivariate-example.eps
end
%\************************************************************************

%/************************************************************************
%Schematic
figure
a=3;b=2;
theta=pi/3;
beta=pi/4;alpha=pi/6;
[kappa,lambda]=ab2kl(a,b);

J1=[1 0 0 ; 0 cos(beta) -sin(beta); 0 sin(beta) cos(beta)];
J3=[cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0 ; 0 0 1];

[x,y,z]=ellsig(kappa,lambda,theta,5*pi/6+[0:.05:2*pi],alpha,beta,'real');plot3(x,y,z,'k','linewidth',2),hold on

plot3([0 x(1)],[0 y(1)],[0 z(1)],'k--','linewidth',2)
plot3(x(1),y(1),z(1),'ko','markersize',10)

[x,y,z]=vectmult(J3*J1,[-1 -1 1 1 -1]*3,[-1 1 1 -1 -1]*3,[0 0 0 0 0]);plot3(x,y,z,'k:')
[x,y,z]=vectmult(J3*J1,[-1 1]*3,[0 0]*3,[0 0]);plot3(x,y,z,'k--')
[x,y,z]=vectmult(J3*J1,[0 0]*3,[-1 1]*3,[0 0]);plot3(x,y,z,'k--')
[x,y,z]=vectmult(J3*J1,[0 cos(theta)]*3,[0 sin(theta)]*3,[0 0]);
h=plot3(x,y,z);linestyle -h h 2E--

[x,y,z]=ellsig(1,0,theta,[0:.05:5*pi/6*1.05],alpha,beta,'real');plot3(x,y,z,'k'),hold on
[x,y,z]=ellsig(1.5,0,theta,[-pi/3:0.05:0],alpha,beta,'real');plot3(x,y,z,'k'),hold on

cxe=1.5*rot(3*pi/2:0.05:2*pi-pi/3);plot3(real(cxe),imag(cxe),3+0*cxe,'k');

cxe=2.5*rot(0:0.05:pi/4);
[x,y,z]=vectmult(J3,0*cxe,-imag(cxe),real(cxe));
plot3(x,y,z,'k')

plot3([0 0],[-3 3],[3 3],'k','linewidth',1)
plot3([-3 3],[0 0 ],[3 3],'k','linewidth',1)
set(gca,'dataAspectRatio',[1 1 1]),set(gca,'plotboxAspectRatio',[1 1 1]),axis equal
axis([-3 3 -3 3 -3 3]),boxon,view(75,15)

xlabel('X-Axis'),ylabel('Y-Axis'),zlabel('Z-Axis')

[x,y,z]=vectmult(J3*J1,[0 0],[0 0],[0 4.24]);plot3(x,y,z,'k','linewidth',2)
h=plot3(x,y,[3 3]);linestyle -h h 2E
plot3([0 0],[0 0 ],[0 3],'k','linewidth',1)
plot3(x(end),y(end),3,'ko','markersize',10,'markerfacecolor','k')

title('Ellipse Schematic in 3D')
text(0.35,-2,3,'$\alpha$')
text(-0.5,-0.7,1.6,'$\beta$')
text(.6,.6,.35,'$\theta$')
text(-1,-1.2,-.25,'$\phi$')

fontsize 12 10 10 10
set(gcf,'paperposition',[1 1 4 4])
  
if ~verLessThan('matlab','8.4.0')
    set(gca,'BoxStyle','full')
end
    
if strcmpi(str,'print')
   print -depsc trivariate-schematic.eps
end
%\************************************************************************


%/************************************************************************
load solomon 
use solomon

vindex(x,y,z,num,5000:12000,1);

%Unit signal vector 
normx=sqrt(x.^2+y.^2+z.^2);
xu=frac(x,normx);
yu=frac(y,normx);
zu=frac(z,normx);

[x,y,z]=anatrans(x,y,z);
t=(num-num(1))*24*3600;

[kappa,lambda,theta,phi,alpha,beta]=ellparams(x,y,z);
[abar,ombar,upbar]=instmom(t(2)-t(1),[x y z],1,2);

%Unit normal vector
mat=pagemtimes(jmat3(alpha,3),jmat3(beta,1));
[xn,yn,zn]=vectmult(mat,0*x,0*x,1+0*x);

%Sphere
[xs,ys,zs]=sphere(20);

ii1=1200:2150;
ii2=2151:6000;

figure
subplot(1,2,1)
h=mesh(xs,ys,zs,1+0*zs,'facealpha',0.7,'markeredgecolor','k');
colormap gray,hold on,axis equal
plot3(xu(ii1),yu(ii1),zu(ii1),'k+')
plot3(xu(ii2),yu(ii2),zu(ii2),'.','color',[1 1 1]*0.4)
view([90-37.5-10, 20])
xlabel('X-Axis'),ylabel('Y-Axis'),zlabel('Z-Axis')
title('Unit Signal Vector')
set(gca,'dataAspectRatio',[1 1 1]),set(gca,'plotboxAspectRatio',[1 1 1]),axis equal
set(gca,'xtick',[-1 -1/2 0 1/2 1])
set(gca,'ytick',[-1 -1/2 0 1/2 1])
set(gca,'ztick',[-1 -1/2 0 1/2 1])
plot3([0 0],[0 0],[0 1],'k','linewidth',3)
plot3([0 0],[0 1],[0 0],'k','linewidth',3)
plot3([0 1],[0 0],[0 0],'k','linewidth',3)
plot3([0 0],[0 0],[0 -1],'k--','linewidth',3)
plot3([0 0],[0 -1],[0 0],'k--','linewidth',3)
plot3([0 -1],[0 0],[0 0],'k--','linewidth',3)
text(-.95,-.95,-.9,'(a)')
plot3(1.05*real(rot(jdeg2rad(12.3-90))),1.05*imag(rot(jdeg2rad(12.3-90))),0,'ko','markersize',12,'markerfacecolor','w')
plot3(1.05*real(rot(jdeg2rad(12.3-90))),1.05*imag(rot(jdeg2rad(12.3-90))),0,'wo','markersize',8,'markerfacecolor','k')

subplot(1,2,2)
h=mesh(xs,ys,zs,1+0*zs,'facealpha',0.7,'markeredgecolor','k');
colormap gray,hold on,axis equal
plot3(xn(ii1),yn(ii1),zn(ii1),'k+')
plot3(xn(ii2),yn(ii2),zn(ii2),'k.','color',[1 1 1]*0.4)
view([90-37.5-10, 20])
xlabel('X-Axis'),ylabel('Y-Axis'),zlabel('Z-Axis')
title('Unit Normal Vector')
set(gca,'dataAspectRatio',[1 1 1]),set(gca,'plotboxAspectRatio',[1 1 1]),axis equal
set(gca,'xtick',[-1 -1/2 0 1/2 1])
set(gca,'ytick',[-1 -1/2 0 1/2 1])
set(gca,'ztick',[-1 -1/2 0 1/2 1])
plot3([0 0],[0 0],[0 1],'k','linewidth',3)
plot3([0 0],[0 1],[0 0],'k','linewidth',3)
plot3([0 1],[0 0],[0 0],'k','linewidth',3)
plot3([0 0],[0 0],[0 -1],'k--','linewidth',3)
plot3([0 0],[0 -1],[0 0],'k--','linewidth',3)
plot3([0 -1],[0 0],[0 0],'k--','linewidth',3)
text(-.95,-.95,-.9,'(b)')
plot3(1.05*real(rot(jdeg2rad(12.3-90))),1.05*imag(rot(jdeg2rad(12.3-90))),0,'ko','markersize',12,'markerfacecolor','w')
plot3(1.05*real(rot(jdeg2rad(12.3-90))),1.05*imag(rot(jdeg2rad(12.3-90))),0,'wo','markersize',8,'markerfacecolor','k')

orient landscape
if strcmpi(str,'print')
    print -depsc trivariate_sphere.eps
end
%\************************************************************************




%/************************************************************************
%Transverse and radial portions
xr=real(real(x)*rot(-jdeg2rad(12.3))+sqrt(-1)*real(y)*rot(-jdeg2rad(12.3)));
xt=imag(real(x)*rot(-jdeg2rad(12.3))+sqrt(-1)*real(y)*rot(-jdeg2rad(12.3)));
xax=[0 1700];

figure,
subplot(4,1,1),plot(t,real([xr xt z])./1e4),yoffset 5,xlim(xax),ylim([-10 37]*.4)
ylabel('Signal x(t) ($\times 10^4$)'),linestyle k k k
title('Characteristics of a Seismic Record')
text(30,2,'r'),text(30,7,'t'),text(30,12,'v')
subplot(4,1,2),plot(t,kappa/1e4,'k'),xlim(xax)
ylabel('Amplitude $\kappa(t)$ ($\times 10^4$)'),ytick([1:5]),ylim([0 22]*2e3/1e4)
linestyle 2k k-. k
subplot(4,1,3),plot(t,lambda,'k'),xlim(xax),ylim([0 1]),ytick([0:.2:1])
ylabel('Linearity $\lambda(t)$')
subplot(4,1,4),plot(t,ombar,'k'),xlim(xax),ylim([0 .7])
for i=1:4
    subplot(4,1,i),vlines(t([1200 2150]),'k:')
    xtick([0:400:2000])
end
ylabel('Frequency $\omega(t)$ (rad/sec)')
letterlabels(2)
packfig(4,1,'rows')
set(gcf,'paperposition',[1 1 4 8])
fontsize 12 10 10 10
if strcmpi(str,'print')
    print -depsc trivariate_seismic.eps
end
%\************************************************************************
