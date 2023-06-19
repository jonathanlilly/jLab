function[varargout]=makefigs_transfer(str)
%MAKEFIGS_TRANSFER Makes all figures for Lilly and Elipot (2021).
%
%   This function makes all figures for 
%
%        Lilly, J. M. and S. Elipot (2021). A unifying perspective on 
%            transfer function solutions to the unsteady Ekman problem. 
%            Fluids, 6 (2): 85, 1--36. 
%
%   This may take a while as some of the plots are computationally
%   intensive.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2021--2023 J.M. Lilly --- type 'help jlab_license' for details


if nargin==0
  str='noprint';
end

if strcmp(str,'--f')
     makefigs_transfer('noprint');return
end

%This is where you want things to be printed to
dirname='/Users/lilly/Desktop/Dropbox/Projects/transfer/theory/figures';

%/*************************************************************************
%basic figure of transfer function
rho=1027;
fc=1e-4;
omega=[-8.5:.01:8.5]*fc;
z=[0:.1:50]';
zo=20;
delta=20;zo=20;h=50;mu=delta.^2./zo;
G1=rho.*fc.*windtrans(omega*24*3600,z,fc*24*3600,delta,mu,inf);
G2=rho.*fc.*windtrans(omega*24*3600,z,fc*24*3600,delta,mu,50);
G1(~isfinite(G1))=1e6;
[ommat,zmat]=meshgrid(omega,z);
zeta=2*sqrt(2)*frac(zo,delta).*sqrt((1+zmat./zo).*abs(1+ommat./fc));

%str='Transfer Function Magnitude (m^2s/kg)';
labelstr='Log10 Transfer Function Magnitude (m^{-1})';

figure,
clf
subplot(1,2,1),contourf(omega./fc,z,log10(abs(G1))',-10:.1:6),nocontours,flipy,hold on
text(-8,48,'(a)','color','w')
subplot(1,2,2),contourf(omega./fc,z,log10(abs(G2))',-10:.1:6),nocontours,flipy,hold on
text(-8,48,'(b)','color','w')

for i=1:2
    subplot(1,2,i)
    caxis([-4 -1]),hlines(20,'0.5k:'),%vlines(-1,'w'),
    %hlines([10 40],'2w')
    if i==1
        hlines([10 40],'0.5D')
    else
        hlines([40],'0.5D')
    end
    xtick([-10:2:10])
    ylabel('Depth (m)')
    xlabel('Nondimensional Frequency $\omega/|f|$','interpreter','latex')
    %contour(omega./fc,z,zeta,[4 7],'w','linewidth',2)
    contour(omega./fc,z,zeta,[4 7.5],'w')
end
subplot(1,2,1)
text(-1,5,'I','color',[1 1 1],'HorizontalAlignment','center')
text(-1,25,'IV','color',[1 1 1],'HorizontalAlignment','center')
text(-1,45,'VII','color',[1 1 1],'HorizontalAlignment','center')
text(2,5,'II','color',[1 1 1]*0,'HorizontalAlignment','center')
text(-4,5,'II','color',[1 1 1]*0,'HorizontalAlignment','center')
text(1,25,'V','color',[1 1 1]*0,'HorizontalAlignment','center')
text(-3,25,'V','color',[1 1 1]*0,'HorizontalAlignment','center')
text(0.4,45,'VIII','color',[1 1 1]*0,'HorizontalAlignment','center')
text(-2.4,45,'VIII','color',[1 1 1]*0,'HorizontalAlignment','center')
text(5.5,5,'III','color',[1 1 1]*0,'HorizontalAlignment','center')
text(-7.5,5,'III','color',[1 1 1]*0,'HorizontalAlignment','center')
text(5.5,25,'VI','color',[1 1 1]*0,'HorizontalAlignment','center')
text(-7.5,25,'VI','color',[1 1 1]*0,'HorizontalAlignment','center')
text(5.5,45,'IX','color',[1 1 1],'HorizontalAlignment','center')
text(-7.5,45,'IX','color',[1 1 1],'HorizontalAlignment','center')

subplot(1,2,2)
%text(-1,5,'I-$h$','color',[1 1 1],'HorizontalAlignment','center')
text(-1,25,'IV-$h$','color',[1 1 1],'HorizontalAlignment','center')
text(-1,45,'VII-$h$','color',[1 1 1]*0,'HorizontalAlignment','center')
%text(2,5,'II-$h$','color',[1 1 1]*0,'HorizontalAlignment','center')
%text(-4,5,'II-$h$','color',[1 1 1]*0,'HorizontalAlignment','center')
text(1,25,'V-$h$','color',[1 1 1]*0,'HorizontalAlignment','center')
text(-3,25,'V-$h$','color',[1 1 1]*0,'HorizontalAlignment','center')
text(0.4,45,'VIII-$h$','color',[1 1 1]*0,'HorizontalAlignment','center')
text(-2.4,45,'VIII-$h$','color',[1 1 1]*0,'HorizontalAlignment','center')
%text(5.5,5,'III-$h$','color',[1 1 1]*0,'HorizontalAlignment','center')
%text(-7.5,5,'III-$h$','color',[1 1 1]*0,'HorizontalAlignment','center')
text(5.5,25,'VI-$h$','color',[1 1 1]*0,'HorizontalAlignment','center')
text(-7.5,25,'VI-$h$','color',[1 1 1]*0,'HorizontalAlignment','center')
text(5.5,45,'IX-$h$','color',[1 1 1],'HorizontalAlignment','center')
text(-7.5,45,'IX-$h$','color',[1 1 1],'HorizontalAlignment','center')

ha=packfig(1,2,'col');
axes(ha(1)),%ylim([0 62]),ytick([0:10:50])
hc=colorbar('SouthOutside');
hc.Label.String=labelstr;
pos1=get(ha(1),'position');
pos2=get(ha(2),'position');
pos=get(hc,'position');%pos(4)=pos(4)/2;
%set(hc,'position',[pos(1)+0.195 pos(2) pos(3) pos(4)])
set(hc,'position',[pos(1)+0.195 pos(2)+.02 pos(3) pos(4)/2])
set(ha(2),'position',[pos2(1) pos1(2) pos2(3) pos1(4)])
set(ha(1),'position',pos1)
fontsize 10 10 9 8
set(gcf,'paperposition',[1 1 10 6])
if strcmp(str,'print')
    jprint(dirname,'transferfunctionschematic')
%    jprint(dirname,'transferfunctionschematic-highres','-r500')
end
% %--------------------------------------------------------------------------
% %Green's function 
% rho=1027;
% fc=1e-4;
% %omega=[-100:.01:100]*fc;
% omega=[-1000:.01:1000]*fc;
% z=[0:.1:49.9]';
% zo=20;
% delta=20;zo=20;h=50;mu=delta.^2./zo;
% %G1=rho.*fc.*windtrans(omega*24*3600,z,fc*24*3600,delta,mu,inf);
% G2=rho.*fc.*windtrans(omega*24*3600,z,fc*24*3600,delta,mu,50);
% g1=ifft(ifftshift(G1,1));
% g2=ifft(ifftshift(G2,1));
% 
% figure,
% jpcolor(1:1e4,z,log10(abs(g2(1:1e4,:)))'),caxis([-7 -3.5]),flipy,hold on
% contour(1:1e4,z,real(g2(1:1e4,:))',[0 0],'k')
%\*************************************************************************
  
%/*************************************************************************
%basic figure of transfer function, flower version
rho=1027;
fc=1e-4;
omega=[-8.5:.01:8.5]*fc;
z=[0:5:45]';
zo=20;
delta=20;zo=20;h=50;mu=delta.^2./zo;
lnomega=logspace(2,100,10000);
omega=[-fliplr(lnomega) -100:.1:-8.6 -8.5:.001:-1-0.005 -1-0.005:1e-7:-1+0.005  -1+0.005:0.001 :8.5 8.6:.1:100 lnomega]*fc;

G1=rho.*abs(fc).*windtrans(omega*24*3600,z,fc*24*3600,delta,mu,inf);
G2=rho.*abs(fc).*windtrans(omega*24*3600,z,fc*24*3600,delta,mu,50);
figure
clf,
subplot(1,2,1)
hl=plot(G1);hold on,linestyle(hl,'0.5D'),
hl=plot(G1(:,1)); linestyle(hl,'2k')
hl=plot(G1(:,end));linestyle(hl,'2k--')
subplot(1,2,2)
hl=plot(G2);hold on,linestyle(hl,'0.5D'),
hl=plot(G2(:,1)); linestyle(hl,'2k')
hl=plot(G2(:,end));linestyle(hl,'2k--')
for i=1:2
    subplot(1,2,i)
    axis square, axis equal
    if i==1
        axis([-.1 1.5 -.8 .8]),
        xtick([-0.8:.2:1.6])
        text(1.35,0.7,['(' setstr(96+i) ')'])
    elseif i==2
        axis([-.01 .13 -.07 .07])
        text(0.117,0.06,['(' setstr(96+i) ')'])
    end
          
    vlines(0,'0.5k:'),hlines(0,'0.5k:')
    xlabel('Real Part of $\widetilde G(\omega,z)$ (m$^{-1}$)','interpreter','latex')%(Dimensionless)')
    ylabel('Imaginary Part of $\widetilde G(\omega,z)$ (m$^{-1}$)','interpreter','latex')% (Dimensionless)')
    %xtick([-.4:.2:1.2]),ytick([-.6:.2:.6])
end
set(gcf,'paperposition',[1 1 10 5])
if strcmp(str,'print')
    jprint(dirname,'schematicflowers')
end
%\*************************************************************************


%/*************************************************************************
%green's functions for schematic
rho=1027;
fc=1e-4;
%omega=[-100:.0001:100]*fc;
omega=[-10000:.01:10000]*fc;
z=[0:1/2:49.5]';
zo=20;
delta=20;zo=20;h=50;mu=delta.^2./zo;

[G1,G2]=vzeros(length(omega),length(z));

%z=[0:5:45]';
parfor i=1:length(z)
    i
%   G1(:,i)=rho.*abs(fc).*windtrans(omega*24*3600,z(i),fc*24*3600,delta,mu,inf);
    G2(:,i)=rho.*abs(fc).*windtrans(omega*24*3600,z(i),fc*24*3600,delta,mu,50);
end
G2(:,end+1)=1e-16+0*G2(:,end);
z=[z;50];
%g1=ifft(ifftshift(G1,1));
g2=ifft(ifftshift(G2,1));
g2=g2./maxmax(abs(g2));


%fc=2 pi rad/ 10^4 seconds
%fN= 2 pi rad
%dt= 1/2 seconds
%Tf = 2 pi / fc = 10^4 seconds 
t=[1:size(g2,1)]'*frac(1,2)./10^4;

figure
ii=1:40000;
subplot(3,1,1)
contourf(t(ii),z,log10(abs(g2(ii,:)))',[-10:.1:0]);nocontours,hold on,flipy
contour(t(ii),z,log10(abs(g2(ii,:)))',[-12:1:0],'w');
hlines(z(1:10:end),'0.2D')
hlines(z(1),'2k')
hlines(z(end-10),'2k--')
%contour(t(ii),z,real(g2(ii,:))',[0 0],'k:');
caxis([-5.5 -1.5]),ylim([-.01 49.999]),xlim([0 2]),xtick([0:.25:2])
vlines(0.25005+[0:1/2:2],'0.25k:')
ylabel('Depth (m)')
text(1.875,2.5,'(a)','color','k')
subplot(3,1,2)
h=plot(t(ii),real(g2(ii,1:10:end))*1000);ylim(1000*[-.009 .009]+1000*0.005)
linestyle(h,'0.5D'),linestyle(h(1),'2k'),linestyle(h(end-1),'2k--')
ylabel('Real Part of $g(t,z)$ ($\times 1000 / \mathrm{max}\{g\}$)','interpreter','latex')
text(1.875,12.7,'(b)','color','k')
hold on,vlines(0.25005+[0:1/2:2],'0.25k:'),xtick([0:.25:2])
%h=plot(t(ii),abs(g2(ii,1:10:end))*1000);
subplot(3,1,3)
h=plot(t(ii),imag(g2(ii,1:10:end))*1000);ylim(1000*[-.009 .009])
linestyle(h,'0.5D'),linestyle(h(1),'2k'),linestyle(h(end-1),'2k--')
ylabel('Imaginary Part of $g(t,z)$ ($\times 1000 / \mathrm{max}\{g\}$)','interpreter','latex')
xlabel('Time (Inertial Periods $2\pi / f$)','interpreter','latex')
text(1.875,8,'(c)','color','k')
hold on,vlines(0.25005+[0:1/2:2],'0.25k:'),xtick([0:.25:2])
packfig(3,1)

set(gcf,'paperposition',[1 1 5 10])
if strcmp(str,'print')
    jprint(dirname,'greensfunction')
end
%\*************************************************************************

if 0
%/*************************************************************************
%green's functions for schematic
rho=1027;
fc=1e-4;
%omega=[-100:.0001:100]*fc;
omega=[-10000:.01:10000]*fc;
z=[0:1/2:49.5]';
delta=20;zo=20;h=50;mu=delta.^2./zo;
delta=10;zo=10;h=50;mu=delta.^2./zo;

[G1,G2]=vzeros(length(omega),length(z));

%z=[0:5:45]';
parfor i=1:length(z)
    i
   % G1(:,i)=rho.*abs(fc).*windtrans(omega*24*3600,z(i),fc*24*3600,delta,mu,inf);
    G2(:,i)=rho.*abs(fc).*windtrans(omega*24*3600,z(i),fc*24*3600,delta,mu,50);
end
%G1(~isfinite(G1))=0;
G2(:,end+1)=1e-16+0*G2(:,end);
z=[z;50];
%g1=ifft(ifftshift(G1,1));
g2=ifft(ifftshift(G2,1));
g2=g2./maxmax(abs(g2));


%fc=2 pi rad/ 10^4 seconds
%fN= 2 pi rad
%dt= 1/2 seconds
%Tf = 2 pi / fc = 10^4 seconds 
t=[1:size(g2,1)]'*frac(1,2)./10^4;

clf
ii=1:40000;
subplot(3,1,1)
contourf(t(ii),z,log10(abs(g2(ii,:)))',[-12:.1:0]);nocontours,hold on,flipy
contour(t(ii),z,log10(abs(g2(ii,:)))',[-12:1/2:0],'w');
hlines(z(1:10:end),'0.2D'),hlines(z(1),'2k'),hlines(z(end-10),'2k--')
%contour(t(ii),z,real(g2(ii,:))',[0 0],'k:');
caxis([-5.5 -1.5]),ylim([-.01 49.999]),xlim([0 2]),xtick([0:.25:2])
vlines(0.25005+[0:1/2:2],'0.25k:')
ylabel('Depth (m)')
text(1.875,2.5,'(a)','color','k')
subplot(3,1,2)
h=plot(t(ii),real(g2(ii,1:10:end))*1000);ylim(1000*[-.009 .009]+1000*0.005-1.5)
linestyle(h,'0.5D'),linestyle(h(1),'2k'),linestyle(h(end-1),'2k--')
ylabel('Real Part of $g(t,z)$ ($\times 10^3 / \mathrm{max}\{g\}$)','interpreter','latex')
text(1.875,12.7-1.5,'(b)','color','k')
hold on,vlines(0.25005+[0:1/2:2],'0.25k:'),xtick([0:.25:2])
%h=plot(t(ii),abs(g2(ii,1:10:end))*1000);
subplot(3,1,3)
h=plot(t(ii),imag(g2(ii,1:10:end))*1000);ylim(1000*[-.009 .009])
linestyle(h,'0.5D'),linestyle(h(1),'2k'),linestyle(h(end-1),'2k--')
ylabel('Imaginary Part of $g(t,z)$ ($\times 10^3 / \mathrm{max}\{g\}$)','interpreter','latex')
text(1.875,8,'(c)','color','k')
hold on,vlines(0.25005+[0:1/2:2],'0.25k:'),xtick([0:.25:2])
xlabel('Time (Inertial Periods $2\pi / f$)','interpreter','latex')
% subplot(4,1,4)
% %gke=5*fliplr(log10(cumsum(squared(fliplr(g2)),2)));
% gke=5*log10(cumsum(squared(g2),2));
% contourf(t(ii),z,gke(ii,:)',[-50:.5:0]);nocontours,hold on,flipy
% contour(t(ii),z,gke(ii,:)',[-50:5:0],'w');
% hlines(z(1:10:end),'0.2D'),hlines(z(1),'2k'),hlines(z(end-10),'2k--')
% caxis([-40 -5]),ylim([-.01 49.999]),xlim([0 2]),xtick([0:.25:2])
% vlines(0.25005+[0:1/2:2],'0.25k:')
% text(1.875,2.5,'(d)','color','k')
% %ylabel('Kinetic Energy $\int_{-h}^z |g(t,x)|^2 \mathrm{d}x$')
% ylabel('Kinetic Energy $\int_{0}^z |g(t,x)|^2 \mathrm{d}x$')
% xlabel('Time (Inertial Periods $2\pi / f$)','interpreter','latex')
% packfig(4,1)
%set(gcf,'paperposition',[1 1 5 13])

packfig(3,1)
set(gcf,'paperposition',[1 1 5 10])
if strcmp(str,'print')
    jprint(dirname,'greensfunction_broader')
end
%\*************************************************************************
end


%/*************************************************************************
%illustration of problems with overflow
log10delta=linspace(-1,4,100);
log10mu=linspace(-7,4,100);
h=100;
%h=1000;
%h=1e6;
%h=1e8;
fc=1e-4;

Gi1=vzeros(length(log10mu),length(log10delta));
Gi2=vzeros(length(log10mu),length(log10delta));
Gi3=vzeros(length(log10mu),length(log10delta));
Gi4=vzeros(length(log10mu),length(log10delta));
%f=0;
f=2.*fc*24*3600;
%f=fc*24*3600*0.99;
%f=fc*24*3600*8;
tic
for i=1:length(log10mu)
    i
    for j=1:length(log10delta)
     Gi1(i,j)=windtrans(f,15,fc*24*3600,10.^log10delta(j),10.^log10mu(i),h,'two');
     Gi2(i,j)=windtrans(f,15,fc*24*3600,10.^log10delta(j),10.^log10mu(i),h,'far');
     Gi3(i,j)=windtrans(f,15,fc*24*3600,10.^log10delta(j),10.^log10mu(i),h,'general');
     Gi4(i,j)=windtrans(f,15,fc*24*3600,10.^log10delta(j),10.^log10mu(i),h,'expansion');
    end
end
%etime1=toc
%etime2=toc

figure,
clf
subplot(1,3,1),contourf(log10delta,log10mu,log10(abs(Gi3)),[-200:.1:0]),nocontours,hold on
contour(log10delta,log10mu,log10(abs(Gi3)),[-200:10:0],'w'),caxis([-12 maxmax(log10(abs(Gi1)))])
hc=colorbar('South');
hc.Label.String='Log10 Magnitude at $\omega = -2f$';
hc.Label.Interpreter='latex';
pos=get(hc,'position');
set(hc,'position',[pos(1)+pos(1)/10 pos(2) pos(3)+pos(1)/10 pos(4)/2])
subplot(1,3,2),contourf(log10delta,log10mu,log10(abs(Gi1)),[-200:.1:0]),nocontours,hold on
contour(log10delta,log10mu,log10(abs(Gi1)),[-200:10:0],'w'),caxis([-12 maxmax(log10(abs(Gi1)))])
subplot(1,3,3),contourf(log10delta,log10mu,log10(abs(Gi2-Gi4)./abs(Gi4)),[-14:.1:0]),nocontours,hold on
caxis([-14 -3])
hc2=colorbar('North');
hc2.Label.String='Log10 Fractional Error at $\omega = -2f$';
hc2.Label.Interpreter='latex';
pos2=get(hc2,'position');
set(hc2,'position',[pos2(1)-pos2(1)/28 pos2(2)-0.03 pos2(3)+pos(1)/10 pos2(4)/2])

[demat,mumat]=meshgrid(10.^log10delta,10.^log10mu);
for i=1:3
    subplot(1,3,i),%contour(log10delta,log10mu,log10(abs(Gi1-Gi2)./abs(Gi1)),[-3.4 -3.4],'k','linewidth',2)
    text(3.5,3.7,['(' setstr(96+i) ')'],'color','k')
    xlim([-1 3.99])
    axis equal,xlim([-1 3.99]),ylim([-5 4])
    ylabel('Log10 Madsen Depth $\mu$ (m)','interpreter','latex')
    xlabel('Log10 Ekman Depth $\delta$ (m)','interpreter','latex')
    %contour(log10delta,log10mu,log10(demat./mumat*2*sqrt(2)*sqrt(3)),[1 1]*3,'color',[1 1 1]*0.6,'linewidth',2)
    contour(log10delta,log10mu,log10(demat./mumat*2*sqrt(2)*sqrt(3)),[1 1]*2.9,'color',[1 1 1]*0.6,'linewidth',2)
end
packfig(1,3,'columns')
set(hc,'position',[pos(1)+pos(1)/10 pos(2)+.1 pos(3)+pos(1)/10 pos(4)/2])
set(hc2,'position',[pos2(1)-pos2(1)/28 pos2(2)-0.03-.07 pos2(3)+pos(1)/10 pos2(4)/2])
%set(gcf,'paperposition',[1 1 11 5])
set(gcf,'paperposition',[1 1 11 8])
if strcmp(str,'print')
    jprint(dirname,'overflowproblem')
end
%\*************************************************************************


%/*************************************************************************
%symmetry plots
rho=1027;
%log10delta=linspace(-1,4,101);
%log10mu=linspace(-7,4,100);
log10delta=linspace(-4,4,201);
log10mu=linspace(-8,4,200);
%h=10000;
h=[16 1e2 1e3 1e4 1e5 1e6];
%h=[15.01 15.1 16 26 1e2 1e3 1e4 1e5 1e6];
fc=1e-4;


[demat,mumat]=meshgrid(10.^log10delta,10.^log10mu);
zomat=demat.^2./mumat;

A=zeros(size(demat,1),size(demat,2),length(h));
A0=zeros(size(demat,1),size(demat,2),length(h));
for i=1:length(h)
    A(:,:,i)=frac(2,mumat).*log(frac(1+h(i)./zomat,1+15./zomat));
    A0(:,:,i)=frac(2,mumat).*log(1+h(i)./zomat);
end    

%Gf=vzeros(length(log10mu),length(log10delta),length(h));
%Gf0=vzeros(length(log10mu),length(log10delta),length(h));
% for i=1:length(log10mu)
%     i
%     parfor j=1:length(log10delta)
%         for k=1:6
%           Gf(i,j,k)=rho.*fc.*windtrans(-fc*24*3600,15,fc*24*3600,10.^log10delta(j),10.^log10mu(i),h(k));
%          Gf0(i,j,k)=rho.*fc.*windtrans(-fc*24*3600,0,fc*24*3600,10.^log10delta(j),10.^log10mu(i),h(k));
%         end
%     end
% end
% maxmax(abs(A-Gf))
% maxmax(abs(A0-Gf0))
% jpcolor(log10(abs(A(:,:,1)-Gf(:,:,1))./A(:,:,1)))
% %these agree to the numerical noise level

omega=[-3:.001:1]*fc;
log10zo=[-12:.1:12]';
ii=1:10:length(log10zo);

clear dei mui dei0 mui0 dei1 mui1
zo=10.^log10zo;
for j=1:length(h)
    dei{j,1}=sqrt(frac(2*zo,1).*log(frac(1+h(j)./zo,1+15./zo)));
    mui{j,1}=dei{j}.^2./zo;
    dei1{j,1}=sqrt(frac(2*zo,10).*log(frac(1+h(j)./zo,1+15./zo)));
    mui1{j,1}=dei1{j}.^2./zo;
    dei0{j,1}=sqrt(frac(2*zo,1).*log(frac(1+h(j)./zo,1)));
    mui0{j,1}=dei0{j}.^2./zo;
end
%--------------------------------------------------------------------------
figure,
clf
contourf(log10delta,log10mu,log10(A(:,:,3)),[-10:.1:10]),nocontours,hold on
contour(log10delta,log10mu,log10(A(:,:,3)),[0 0 ],'k','linewidth',1.5)
contour(log10delta,log10mu,log10(A(:,:,3)),[-10:1:10],'k','linewidth',0.5)
contour(log10delta,log10mu,log10(zomat),[-20:20],'w','linewidth',0.5)
contour(log10delta,log10mu,log10(zomat),[0 0],'w','linewidth',1.5)
plot(log10(dei{3}(1:10:end)),log10(mui{3}(1:10:end)),'ko','markerfacecolor','k','markersize',4)
caxis([-4.7 7.3]),axis equal,ylim([-7 4]),xlim([-4 4]),
ylabel('Log10 Madsen Depth $\mu$ (m)','interpreter','latex')
xlabel('Log10 Ekman Depth $\delta$ (m)','interpreter','latex')
hc=colorbar('SouthOutside');
hc.Label.String='Log10 Inertial Amplitude A (m$^{-1}$)';
hc.Label.Interpreter='latex';
hc.Ticks=[-4:1:8];
set(gcf,'paperposition',[1 1 4*1.5 5*1.5])
%text(-0.92,3.7,'M'),text(3.8,-6.8,'E')
text(-3.925,3.7,'M'),text(3.7,-6.8,'E','color','w')
if strcmp(str,'print')
    jprint(dirname,'deltamuplane')
end
%--------------------------------------------------------------------------
%figure,cellplot(dec,muc)
%figure,cellplot(dec0,muc0)
omega=[-100:.1:-8.6 -8.5:.001:-1-0.005 -1-0.005:1e-7:-1+0.005  -1+0.005:0.001:8.5 8.6:.1:100]*fc;
%omega=[-1e4:10:-1000 -1000:1:-100 -100:.1:-8.6 -8.5:.001:-1-0.005 -1-0.005:1e-7:-1+0.005  -1+0.005:0.001 :8.5 8.6:.1:100 100:1:1000 1000:10:1e4]*fc;
G=nan*zeros(length(omega),length(ii),length(h));
[Gm,Ge,Gminf,Geinf]=vzeros(length(omega),length(h));
%G0=nan*zeros(length(omega),length(ii),length(h));
%G1=nan*zeros(length(omega),length(ii),length(h));
for j=1:size(G,3)
    j
   Gm(:,j)=rho.*abs(fc).*windtrans(omega*24*3600,15,fc*24*3600,0,mui{j}(1),h(j));
   Ge(:,j)=rho.*abs(fc).*windtrans(omega*24*3600,15,fc*24*3600,dei{j}(end),0,h(j));
   parfor k=1:length(ii)
        k
        G(:,k,j)=rho.*abs(fc).*windtrans(omega*24*3600,15,fc*24*3600,dei{j}(ii(k)),mui{j}(ii(k)),h(j));
         %G1(:,k,j)=rho.*abs(fc).*windtrans(omega*24*3600,15,fc*24*3600,dei1{j}(ii(k)),mui1{j}(ii(k)),h(j));
        %G0(:,k,j)=rho.*abs(fc).*windtrans(omega*24*3600,0,fc*24*3600,dei0{j}(ii(k)),mui0{j}(ii(k)),h(j));
    end
end
%--------------------------------------------------------------------------
figure
clf
%[~,jj]=min(abs(omega));
for j=1:size(G,3)
    subplot(2,3,j)
    hl=plot(G(:,:,j));axis equal, axis square
    axis([-.3 1 -.65 .65]*1.05),vlines(0,'0.5k:'),hlines(0,'0.5k:')
    linestyle(hl,'0.5D')
    hl=plot(G(:,1,j));linestyle(hl,'2k--')
    hl=plot(G(:,end,j));linestyle(hl,'2k')
    hl=plot(Ge(:,j));linestyle(hl,'0.5w')
    hl=plot(Gm(:,j));linestyle(hl,'0.5k')
    text(0.9,0.59,['(' setstr(96+j) ')'])
%    xlabel('Real Part of Transfer Function (m^2s / kg)')%(Dimensionless)')
%    ylabel('Imaginary Part of Transfer Function (m^2s / kg)')% (Dimensionless)') 
    xlabel('Real Part of $\widetilde G(\omega,15$ m$)$ (m$^{-1}$)','interpreter','latex')%(Dimensionless)')
    ylabel('Imaginary Part of $\widetilde G(\omega,15$ m$)$ (m$^{-1}$)','interpreter','latex')% (Dimensionless)')
    xtick([-.4:.2:1.2]),ytick([-.6:.2:.6])
    %plot(G(jj,:,j),'ko','markerfacecolor','k','markersize',4)
end
packfig(2,3)
set(gcf,'paperposition',[1 1 9.4 6])
fontsize 10 10 9 9
if strcmp(str,'print')
    jprint(dirname,'similarflowers')
end
%%check invariance
%figure,plot(G(:,:,1)),hold on,plot(G1(:,:,1)/10,'r')
%figure,plot(G(:,:,end)),hold on,plot(G1(:,:,end)/10,'r')
% %--------------------------------------------------------------------------
% figure
% clf
% for j=1:size(Gf,3)
%     subplot(2,3,j)
%     hl=plot(omega./fc,abs(G(:,2:end-1,j)));hold on
%     axis([-1.99 0 1.001*1e-4 1.25]),ylog 
%     linestyle(hl,'0.5D')
%     hl=plot(omega./fc,abs(G(:,1,j)));
%     linestyle(hl,'2k')
%     hl=plot(omega./fc,abs(G(:,end,j)));
%     linestyle(hl,'2k--')
%     ylabel('Log10 Magnitude of $G(\omega,15$ m$)$','interpreter','latex')
%     xlabel('Nondimensional Frequency $\omega/f$','interpreter','latex')
%     text(-.17,0.8,['(' setstr(96+j) ')'])
%     xtick([-2:.5:0])
%     vlines(-1,'0.5k:')
% end
% packfig(2,3)
% set(gcf,'paperposition',[1 1 9*1.25 5*1.25])
% jprint(dirname,'similarmagnitudes')
% %--------------------------------------------------------------------------
% figure
% clf
% for j=1:size(Gf,3)
%     subplot(2,3,j)
%     hl=plot(omega./fc,(360/2/pi)*unwrap(angle(G(:,2:end-1,j))));hold on
%     axis([-1.99 0 -179.9 180]),
%     linestyle(hl,'0.5D')
%     hl=plot(omega./fc,(360/2/pi)*unwrap(angle(G(:,1,j))));
%     linestyle(hl,'2k')
%     hl=plot(omega./fc,(360/2/pi)*unwrap(angle(G(:,end,j))));
%     linestyle(hl,'2k--')
%     ylabel('Angle of $G(\omega,15$ m$) (degrees) $','interpreter','latex')
%     xlabel('Nondimensional Frequency $\omega/f$','interpreter','latex')
%     %text(-.17,0.8,['(' setstr(96+j) ')'])
%     xtick([-2:.5:0])
%     vlines(-1,'0.5k:')
% end
% packfig(2,3)
% set(gcf,'paperposition',[1 1 9*1.25 5*1.25])
% jprint(dirname,'similarphases')
%--------------------------------------------------------------------------
% figure
% clf
% for j=1:size(Gf0,3)
%     subplot(2,3,j)
%     hl=plot(G0(:,2:end-1,j));axis equal, axis square
%     axis([-.05 1.05 -.505 .505]),vlines(0,'0.5k:'),hlines(0,'0.5k:')
%     linestyle(hl,'0.5D'),
%     hl=plot(G0(:,end,j)); linestyle(hl,'2k')
%     hl=plot(G0(:,1,j));linestyle(hl,'2k--')
%     text(0.9,0.44,['(' setstr(96+j) ')'])
%     xlabel('Real Part of $G(\omega,0$ m$)$','interpreter','latex')%(Dimensionless)')
%     ylabel('Imaginary Part of $G(\omega,0$ m$)$','interpreter','latex')% (Dimensionless)') 
%     xtick([-.4:.2:1.2]),ytick([-.6:.2:.6])
% end
% packfig(2,3)
% set(gcf,'paperposition',[1 1 9.4 6])
% jprint(dirname,'similarflowers0')
lnomega=logspace(2,100,10000);
omega=[-fliplr(lnomega) -100:.1:-8.6 -8.5:.001:-1-0.005 -1-0.005:1e-7:-1+0.005  -1+0.005:0.001 :8.5 8.6:.1:100 lnomega]*fc;
G0=vzeros(length(omega),length(ii));
Gm0=rho.*abs(fc).*windtrans(omega*24*3600,0,fc*24*3600,0,mui0{3}(1),h(3));
Ge0=rho.*abs(fc).*windtrans(omega*24*3600,0,fc*24*3600,dei0{3}(end),0,h(3));
%Gm0=rho.*abs(fc).*windtrans(omega*24*3600,1e-10,fc*24*3600,0,mui0{3}(ii(end)),h(3));
%Ge0=rho.*abs(fc).*windtrans(omega*24*3600,0,fc*24*3600,dei0{3}(ii(1)),0,h(3));
%zeta=2*sqrt(2)*frac(zo,delta).*sqrt(abs(1+ommat./fc));

for k=1:length(ii)
    k
    G0(:,k)=rho.*abs(fc).*windtrans(omega*24*3600,0,fc*24*3600,dei0{3}(ii(k)),mui0{3}(ii(k)),h(3));
 %   G1(:,k)=rho.*abs(fc).*windtrans(omega*24*3600,0,fc*24*3600,dei0{3}(ii(k)),mui0{3}(ii(k)),h(1));
end
figure
clf
hl=plot(G0);axis equal, axis square
axis([-.05 1.05 -.505 .505]),vlines(0,'0.5k:'),hlines(0,'0.5k:')
linestyle(hl,'0.5D'),
hl=plot(G0(:,1)); linestyle(hl,'2k--')
hl=plot(G0(:,end)); linestyle(hl,'2k')
hl=plot(Ge0); linestyle(hl,'0.5w')
hl=plot(Gm0);linestyle(hl,'0.5k')
%text(0.9,0.44,['(' setstr(96+j) ')'])
xlabel('Real Part of $\widetilde G(\omega,0$ m$)$ (m$^{-1}$)','interpreter','latex')%(Dimensionless)')
ylabel('Imaginary Part of $\widetilde G(\omega,0$ m$)$ (m$^{-1}$)','interpreter','latex')% (Dimensionless)')
xtick([-.4:.2:1.2]),ytick([-.6:.2:.6])
set(gcf,'paperposition',[1 1 4 4])
plot(2*h(3)/squared(dei0{3}(end)),0,'wo','markerfacecolor','k')
fontsize 10 10 9 9
if strcmp(str,'print')
    jprint(dirname,'similarflowers0')
end


