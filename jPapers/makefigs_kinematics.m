function[varargout]=makefigs_kinematics(str)
%MAKEFIGS_KINEMATICS  Make figures for Lilly (2018).
%
%   This function makes all figures for the paper 
%
%   Lilly, J. M. (2018) Kinematics of a fluid ellipse in a linear flow. 
%       Fluids, 3 (1) 16: 1--29.
%
%   'makefigs_kinematics' generates all figures.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2018--2023 J.M. Lilly --- type 'help jlab_license' for details

if nargin==0
    str='noprint';
end

%cd /Users/lilly/Desktop/Dropbox/Projects/kinematics
%/************************************************************************
%Ellipse & Annulus Schematic
figure
subplot(1,2,1)
a=3;b=2;
theta=pi/3;
alpha=1.2;
[kappa,lambda]=ab2kl(a,b);

ellipseplot(kappa,lambda,theta,0+sqrt(-1)*0,'npoints',64,'k');hold on
linestyle 2k
[x,y]=ellsig(kappa*alpha,lambda,theta,[0:pi/12:2*pi]','real');

phi=[0:pi/48:2*pi-1e-10]';
[x,y]=ellsig(kappa,lambda,theta,phi,'real');
r=sqrt(x.^2+y.^2);
z=rot(theta).*(a*cos(phi(1:4:end))+1i*b*sin(phi(1:4:end)));
for i=1:length(z)
    h1=plot(real(z(i))*[1 alpha],imag(z(i))*[1 alpha],'k');
end
linestyle -h h1 G
h1=plot(real(z(1))*[0 alpha],imag(z(1))*[0 alpha],'k');
h1=plot(real(z(5))*[0 alpha],imag(z(5))*[0 alpha],'k');
h1=plot(real(z(1))*[0 1],imag(z(1))*[0 1]);linestyle -h h1 3G
h1=plot(real(z(7))*[0 1],imag(z(7))*[0 1]);linestyle -h h1 3G--
plot(z(5),'wo','markerfacecolor','w','markersize',9)
plot(z(5),'ko','markerfacecolor','k','markersize',6)

% vartheta=atan2(imag(z),real(z));  %Azimuth angle
% [vartheta,index]=sort(vartheta);
% R=interp1([vartheta;vartheta(1)+2*pi],[abs(z(index));abs(z(index(1)))],angle(rot(phi(1:4:end))));
%
% z2=rot(phi(1:4:end)).*R;
% for i=1:length(z2)
%     if i~=1&&i~=7&&i~=13&&i~=19 %Omit lines that vlines+hlines will plot
%         h1=plot(real(z2(i))*[0 1],imag(z2(i))*[0 1],':');
%         linestyle -h h1 D:
%     end
% end
hlines(0,'k:'),vlines(0,'k:')

plot(rot(0:pi/48:angle(z(1))),'k')
phi=0:pi/48:angle(z(1));
plot(rot(theta)*(a.*cos(phi)+1i*b.*sin(phi))/2,'k')
text(.5,.25,'\theta','interpreter','tex')
text(.1,.75,'\phi','interpreter','tex')
text(1.15,1.75,'a')
text(-1.15,0.4,'b')
axis([-3.25 3.25 -3.75 3.75])
title('Ellipse Schematic')
axis off
%************************************************************************
%Annulus Schematic
subplot(1,2,2)
a=3;b=2;
theta=pi/3;
alpha=1.2;
[kappa,lambda]=ab2kl(a,b);

ellipseplot(kappa,lambda,theta,0+sqrt(-1)*0,'npoints',64,'k');hold on
ellipseplot(kappa*alpha,lambda,theta,0+sqrt(-1)*0,'npoints',64,'k');
linestyle 2k
[x,y]=ellsig(kappa*alpha,lambda,theta,[0:pi/12:2*pi]','real');

phi=[0:pi/48:2*pi-1e-10]';
[x,y]=ellsig(kappa,lambda,theta,phi,'real');
r=sqrt(x.^2+y.^2);
z=rot(theta).*(a*cos(phi(1:4:end))+1i*b*sin(phi(1:4:end)));

for i=1:length(z)
    h1=plot(real(z(i))*[1 alpha],imag(z(i))*[1 alpha],'k');
end
vartheta=atan2(imag(z),real(z));  %Azimuth angle
[vartheta,index]=sort(vartheta);
R=interp1([vartheta;vartheta(1)+2*pi],[abs(z(index));abs(z(index(1)))],angle(rot(phi(1:4:end))));

for ii=17:20;
    patch([x(ii) alpha*x(ii) alpha*x(ii+1) x(ii+1)],...
        [y(ii) alpha*y(ii) alpha*y(ii+1) y(ii+1)],'k')
end

% z2=rot(phi(1:4:end)).*R;
% for i=1:length(z2)
%     if i~=1&&i~=7&&i~=13&&i~=19 %Omit lines that vlines+hlines will plot
%         h1=plot(real(z2(i))*[0 1],imag(z2(i))*[0 1],':');
%         linestyle -h h1 D:
%     end
% end
hlines(0,'k:'),vlines(0,'k:')

h1=plot(real(z(1))*[0 alpha],imag(z(1))*[0 alpha],'k');
h1=plot(real(z(5))*[0 alpha],imag(z(5))*[0 alpha],'k');
h1=plot(real(z(1))*[0 1],imag(z(1))*[0 1]);linestyle -h h1 3G
h1=plot(real(z(7))*[0 1],imag(z(7))*[0 1]);linestyle -h h1 2w
h1=plot(real(z(7))*[0 1],imag(z(7))*[0 1]);linestyle -h h1 3G--

%hlines(0,'k:'),vlines(0,'k:')
plot(rot(0:pi/48:angle(z(1))),'k')
phi=0:pi/48:angle(z(1));
plot(rot(theta)*(a.*cos(phi)+1i*b.*sin(phi))/2,'k')
text(.5,.25,'\theta','interpreter','tex')
text(.1,.75,'\phi','interpreter','tex')
text(1.15,1.75,'a')
text(-1.15,0.4,'b')

%legend([h1 h2 h3],'Time','Phase','Azimuth')
axis([-3.25 3.25 -3.75 3.75])
title('Elliptical Ring Schematic')
fontsize 14 12 12 12
axis off
letterlabels(4)

packfig(1,2,'columns')

set(gcf,'paperposition',[1 1 12 8])
set(gcf,'color','w')
set(gcf,'inverthardcopy','off')

if strcmp(str,'print')
    print -deps kinematics-schematic.eps
end
%\************************************************************************


%/************************************************************************
figure,
x=[-7:1:7]';y=x;
[xg,yg]=meshgrid(x,y);

% subplot(1,4,1);quiver(x,y,xg,yg,'k','linewidth',1),title('Divergence')
% subplot(1,4,2);quiver(x,y,-yg,xg,'k','linewidth',1),title('Vorticity')
% subplot(1,4,3);quiver(x,y,xg,-yg,'k','linewidth',1),title(' Normal Strain')
% subplot(1,4,4);quiver(x,y,yg,xg,'k','linewidth',1),title('Shear Strain')

for i=1:4,subplot(1,4,i), axis equal,axis([-1 1 -1 1]*max(x)),
    %vlines(0,'k:'),hlines(0,'k:'),
    noxlabels,noylabels,hold on,boxon,end

subplot(1,4,1);quiver(x,y,xg,yg,'k','linewidth',0.5),title('Ix','FontWeight','bold')
subplot(1,4,2);quiver(x,y,-yg,xg,'k','linewidth',0.5),title('Jx','FontWeight','bold')
subplot(1,4,3);quiver(x,y,xg,-yg,'k','linewidth',0.5),title('Kx','FontWeight','bold')
subplot(1,4,4);quiver(x,y,yg,xg,'k','linewidth',0.5),title('Lx','FontWeight','bold')

packfig(1,4,'columns')

orient landscape
fontsize 9 9 9 9
if strcmp(str,'print')
    print -depsc kinematics-ijkl-vectors.eps
end

figure
x=[-10:.1:10]';y=x;
[xg,yg]=meshgrid(x,y);
subplot(1,4,1);contour(x,y,xg.^2+yg.^2,10,'k'),title('x$^T$Ix','fontweight','bold'),
subplot(1,4,2);title('x$^T$Jx','fontweight','bold')
subplot(1,4,3);
z=xg.^2-yg.^2;zp=z;zn=z;zp(z<0)=nan;zn(z>0)=nan;
contour(x,y,zp,10,'k'),title('x$^T$Kx','fontweight','bold'),
hold on, contour(x,y,zn,10,'-.','color',0.6*[1 1 1]),
subplot(1,4,4);
z=2*xg.*yg;zp=z;zn=z;zp(z<0)=nan;zn(z>0)=nan;
contour(x,y,zp,10,'k'),title('x$^T$Lx','fontweight','bold'),
hold on,contour(x,y,zn,10,'-.','color',0.6*[1 1 1]),
for i=1:4,subplot(1,4,i), axis equal,axis([-10 10 -10 10]),boxon,
    %vlines(0,'D:'),hlines(0,'D:'),
    noxlabels,noylabels,end
packfig(1,4,'columns')


orient landscape
fontsize 9 9 9 9
if strcmp(str,'print')
    print -depsc kinematics-ijkl-quadratic.eps
end
%\************************************************************************


