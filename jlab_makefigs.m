function[varargout]=jlab_makefigs(namestr,str)
%JLAB_MAKEFIGS  Makes figures for papers by J.M. Lilly.
%
%   JLAB_MAKEFIGS NAME makes all figures for the publication NAME, as follows. 
%
%  'jlab_makefigs gulfcensus':
%   Lilly, J. M. and P. Perez-Brunius (2021).  Extracting statistically 
%       significant eddy signals from large Lagrangian datasets using 
%       wavelet ridge analysis, with application to the Gulf of Mexico.
%       Submitted to Nonlinear Processes in Geophysics.
%
%  'jlab_makefigs gulfdrifters':
%   Lilly, J. M. and P. Perez-Brunius (2021). A gridded surface current
%       product for the Gulf of Mexico from consolidated  drifter
%       measurements.  Accepted to Earth Science System Data.
%
%  'jlab_makefigs kinematics':
%   Lilly, J. M. (2018) Kinematics of a fluid ellipse in a linear flow. 
%       Fluids, 3 (1) 16: 1--29.
%
%  'jlab_makefigs matern':
%   Lilly, J. M., A. M. Sykulski, J. J. Early, and S. C. Olhede (2017).
%       Fractional Brownian motion, the Matern process, and stochastic 
%       modeling of turbulent dispersion.  Nonlinear Processes in
%       Geophysics, 24: 481--514.
%
%   'jlab_makefigs element':
%   Lilly, J. M.  (2017).  Element analysis: a wavelet-based method for
%       analyzing time-localized events in noisy time series.  Proceedings 
%       of the Royal Society of London, Series A, 473 (2200), 1--28.
%
%   'jlab_makefigs superfamily':
%   Lilly, J. M., and S. C. Olhede (2012). Generalized Morse wavelets as a
%      superfamily of analytic wavelets.  IEEE Transactions on Signal
%      Processing 60 (11), 6036--6041.
%
%   'jlab_makefigs multivariate':
%   Lilly, J. M., and S. C. Olhede (2012). Analysis of modulated 
%      multivariate oscillations. IEEE Transactions on Signal Processing, 
%      60 (2), 600--612.
%
%   'jlab_makefigs vortex': 
%   Lilly, J. M., R. K. Scott, and S. C. Olhede (2011). Extracting waves 
%      and vortices from Lagrangian trajectories. Geophysical Research 
%      Letters, 38, L23605, 1--5.
%
%   'jlab_makefigs trivariate':
%   Lilly, J. M. (2011). Modulated oscillations in three dimensions. IEEE
%      Transactions on Signal Processing, 59 (12), 5930--5943.
%
%   'jlab_makefigs analytic':
%   Lilly, J. M., and S. C. Olhede (2010). On the analytic wavelet 
%      transform. IEEE Transactions on Information Theory, 56 (8),
%      4135--4156.
%
%   'jlab_makefigs bandwidth':
%   Lilly, J. M., and S. C. Olhede (2010). Bivariate instantaneous 
%      frequency and bandwidth. IEEE Transactions on Signal Processing,
%      58 (2), 591--603.
%
%   'jlab_makefigs asilomar':
%   Lilly, J. M., and S. C. Olhede (2009). Wavelet ridge estimation of 
%      jointly modulated multivariate oscillations. IEEE 43rd Asilomar 
%      Conference on Signals, Systems, and Computers, 452--456. 
%
%   'jlab_makefigs morsies':
%   Lilly, J. M., and S. C. Olhede, (2009). Higher-order properties of 
%      analytic wavelets. IEEE Transactions on Signal Processing, 57 (1),
%      146--160.
%   
%   'jlab_makefigs ridges':
%   Lilly, J. M., and J.-C. Gascard (2006). Wavelet ridge diagnosis of 
%      time-varying elliptical signals with application to an oceanic eddy. 
%      Nonlinear Processes in Geophysics 13, 467--483.
%   
%   'jlab_makefigs all' makes the figures for all papers.
%  
%   'jlab_makefigs NAME print' will  print all figure for paper NAME as
%   .eps files in the current diretory, e.g. 'jlab_makefigs ridges print'. 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2020 J.M. Lilly --- type 'help jlab_license' for details     
 

if nargin==1
    str='noprint';
end

jj=0;
jj=jj+1;names{jj}='ridges';
jj=jj+1;names{jj}='morsies';
jj=jj+1;names{jj}='asilomar';
jj=jj+1;names{jj}='bandwidth';
jj=jj+1;names{jj}='analytic';
jj=jj+1;names{jj}='trivariate';
jj=jj+1;names{jj}='vortex';
jj=jj+1;names{jj}='multivariate';
jj=jj+1;names{jj}='superfamily';
jj=jj+1;names{jj}='element';
jj=jj+1;names{jj}='matern';
jj=jj+1;names{jj}='kinematics';
jj=jj+1;names{jj}='gulfdrifters';
jj=jj+1;names{jj}='gulfcensus';

%cd(jlab_settings('dirnames.figures'))

dti=get(0,'defaultTextInterpreter');
dli=get(0,'defaultTextInterpreter');
set(0,'defaultTextInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')


if strcmpi(namestr(1:3),'all')
    for i=1:length(names)
        eval(['jlab_makefigs_' names{i} '(''' str ''');'])
    end
else
    if exist(['jlab_makefigs_' namestr])
        eval(['jlab_makefigs_' namestr '(''' str ''');'])
    else
        disp(['Sorry, formula for paper ''' namestr ''' not found.'])
    end
end

set(0,'defaultTextInterpreter',dti)
set(0,'defaultLegendInterpreter',dli)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[varargout]=jlab_makefigs_gulfdrifters(str) 
makefigs_gulfdrifters(str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=jlab_makefigs_kinematics(str) 
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

%END of jlab_makefigs_kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[varargout]=jlab_makefigs_matern(str) 
makefigs_matern(str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[varargout]=jlab_makefigs_element(str) 
makefigs_element(str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[varargout]=jlab_makefigs_superfamily(str)

%/************************************************************
figure

N=2^16;
fs=1/512;
t=1:N;t=t-mean(t);t=t.*fs;
clear psi psif 
P=zeros(4,5);
ga=[1/3 1 3 9 27]';
be=flipud(ga(1:end));
for i=1:5
    for j=1:5
        P(i,j)=sqrt(be(j)*ga(i));
        [psi(:,i,j) psif(:,i,j)]=morsewave(N,1,ga(i),be(j),2*pi*fs,'bandpass');  
    end
end
for i=1:5
    for j=1:5
        subplot(5,5,i+(j-1)*5)
        x=[squeeze(real(psi(:,i,j))) squeeze(imag(psi(:,i,j))) squeeze(abs(psi(:,i,j)))]./max(abs(psi(:,i,j)));
        plot(t./(P(i,j)./pi),x),xlim([-2 2]),ylim([-0.85 1.1])
        linestyle k k-- 1.5k
        hlines(0,'k:'),vlines([-1/2 1/2],'k:')
        if iseven(i+(5-j+1)*5),set(gca,'color',[1 1 1]*0.9),end
        if j==1,title(['$\gamma=$' num2str(ga(i),3)]),end
        if i==1,ylabel(['$\beta=$' num2str(be(j),3)]),end
        set(gca,'xtick',[]),set(gca,'ytick',[]),set(gca,'xcolor','w'),set(gca,'ycolor','w')
        set(get(gca,'ylabel'),'color','k')
        if j==5&&i==1,text(-1.8,-0.6,'(a)'),end
    end
end
packfig(5,5)

fontsize 14 14 14 14
orient tall
set(gcf,'paperposition',[1 1 8 8])

if strcmpi(str,'print')
    print -deps morsie_families_new.eps
end
%\*********************************************************


%/************************************************************
figure

N=2^15;
fs=1/512;
t=1:N;t=t-mean(t);t=t.*fs;
clear psi psif 
P=zeros(4,5);
ga=[1/3 1 3 9 27]';
be=flipud(ga(1:end));
for i=1:5
    for j=1:5
        P(i,j)=sqrt(be(j)*ga(i));
        [psi(:,i,j) psif(:,i,j)]=morsewave(N,1,ga(i),be(j),2*pi*fs,'bandpass');  
    end
end


f=frac(1,2*pi)*fourier(N);
for i=1:5
    for j=1:5
        subplot(5,5,i+(j-1)*5)
        psigauss=simplepdf(f,fs,fs./P(i,j),'gaussian');
        psigauss=psigauss./maxmax(psigauss).*2;
        fbg=-frac(1,6)*P(i,j).^2.*(f./fs-1).^3.*(ga(i)-3)-frac(1,24)*P(i,j).^2.*(f./fs-1).^4.*((ga(i)-3).^2+2);
        plot(f,[psif(1:length(f),i,j) psigauss psigauss.*exp(fbg)]),xlim([0 f(129+64)]),ylim([0 2.1])
        linestyle 2k k-.  k 
        vlines(1/512,'k:')
        if iseven(i+(5-j+1)*5),set(gca,'color',[1 1 1]*0.9),end
        if j==1,title(['$\gamma=$' num2str(ga(i),3)]),end
        if i==5,ylabel(['$\beta=$' num2str(be(j),3)]),end
        %if i==1,vlines(0,'k'),end
        if j==5&&i==1,text(.0005,0.25,'(b)'),end
    end
end
h=packfig(5,5);
h=reshape(h,5,5);
for i=1:5
    for j=1:5
        axes(h(i,j))
        if i==5,ylabel(['$\beta=$' num2str(be(j),3)]),end
        set(gca,'xtick',[]),set(gca,'ytick',[]),set(gca,'xcolor','w'),set(gca,'ycolor','w')
        set(get(gca,'ylabel'),'color','k'),set(gca,'yaxislocation','right')
    end
end

fontsize 14 14 14 14
orient tall
set(gcf,'paperposition',[1 1 8 8])

if strcmpi(str,'print')
    print -deps morsie_families_fourier.eps
end
%\*********************************************************




%/*********************************************************************
ga1=logspace(log10(.1),2,100);
be1=logspace(log10(1/4+.001),2,101);

[ga,be]=meshgrid(ga1,be1);
[fm,fe,fi,cf] = morsefreq(ga,be);
[a,sigt,sigo,skew]=morsebox(ga,be);
p=sqrt(be.*ga);

figure
contourf(ga1,be1,a,[1/2:.05:1.5]),xlog,ylog,nocontours
hold on, contour(ga1,be1,a,(.500:.005:.55),'k'),colormap gray,flipmap
hold on, contour(ga1,be1,skew,[0 0],'k--'),colormap gray,flipmap

vlines([1 2 3],'k:'),
plot(ga1,frac(ga1-1,2),'k','linewidth',2)
for i=-1:3
    h=plot(ga1,(3.^(2*i))./ga1);linestyle -h h 2w
    h=plot(ga1,(3.^(2*i))./ga1);linestyle -h h E
end



%contour(ga1,be1,(k3./k2.^(3/2)),[0 0],'k')
caxis([1/2 1.5])
axis([1/10 100 1/4 100])

xtick([1/2 1 2 3 4 5 6 8 16 32])
ytick([1/2 1 2 4 8 16 32 64 128])

if verLessThan('matlab','8.4.0')
    xtl=get(gca,'xticklabel');xtl(1,:)='1/2';set(gca,'xticklabel',xtl)
    ytl=get(gca,'yticklabel');ytl(1,:)='1/2';set(gca,'yticklabel',ytl)
else
    xtl=get(gca,'xticklabel');xtl{1}='1/2';set(gca,'xticklabel',xtl)
    ytl=get(gca,'yticklabel');ytl{1}='1/2';set(gca,'yticklabel',ytl)
end
outticks

title('Generalized Morse Wavelet Phase Diagram')
xlabel('Gamma Parameter'),ylabel('Beta Parameter')

text(.12+.25,17.86,'B')
text(1.25,90,'C')
text(2.25,90,'G')
text(3.25,90,'A')
text(1/2,85,'e')
text(64,0.28,'S')
text(.30,.28,'a')

ga=[1/3 1 3 9 27]';
for i=1:5,for j=1:5,plot(ga(i),ga(j),'wo','markersize',4,'markerfacecolor','k'),end,end
plot(1/10,22,'ks','markersize',5,'markerfacecolor','k')

hc=colorbar;hc=discretecolorbar(hc,[.5 1.5],(1/2:.05:1.5));nocontours
axes(hc),hlines(.500:.005:.55,'k')
ylabel('Heisenberg Area'),ytick(1/2:.1:1.5),fixlabels([0 -1])



orient tall
fontsize 12 12 12 12
set(gcf,'paperposition',[1 1 5 5])


if strcmpi(str,'print')
   print -depsc analytic-morsephase.eps
end
%/************************************************************************


%/********************************************************************
%Compute Morlet area and projection onto Gaussian
disp('Sorry, this next bit takes a while....')

om=logspace(log10(1.1),log10(2.5*pi),100);
nu=morlfreq(om);
p_morlet=om.*sqrt(om.*(om-nu)+1);
N=2^13;
t=[0:N-1]';t=t-mean(t);
L=50;

clear sigmat sigmao psi g
for i=1:length(om)
        [psii,psifi]=morlwave(N,om(i),1./L,'bandpass');
        [mut,sigmat(i)]=pdfprops(t,abs(psii).^2); 
        [muo,sigmao(i)]=pdfprops(2*pi*fftshift(t./length(t)),abs(psifi).^2); 
        %Demodulated version
        psi(:,i)=morlwave(N,om(i),1/L,'energy').*rot(-1/L.*t);
        g(:,i)=simplepdf(t,0,p_morlet(i)*L,'gaussian');
        g(:,i)=g(:,i)./sqrt(sum(abs(g(:,i)).^2));
end
proj_morlet=abs(squeeze(vsum(conj(psi).*g,1)));


%ga1=[2.99 3 3.01];
ga1=(1:1:6);
be1=logspace(-10,2,200);
[ga,be]=meshgrid(ga1,be1);
[fm,fe,fi,cf] = morsefreq(ga,be);
a=morsebox(ga,be);
[p,skew,kurt]=morseprops(ga,be);

clear psi_morse g
for i=1:length(be1)
    for j=1:length(ga1)
        psi_morse(:,i,j)=morsewave(N,1,ga1(j),be1(i),1/L,'energy').*rot(-1/L.*t);
        g(:,i,j)=simplepdf(t,0,p(i,j)*L,'gaussian');
        g(:,i,j)=g(:,i,j)./sqrt(sum(abs(g(:,i,j)).^2));
    end
end
proj_morse=abs(squeeze(vsum(conj(psi_morse).*g,1)));

figure, 
%subplot(2,1,1)
plot(p/pi,1./vswap(a,inf,nan))
hold on, plot(p_morlet/pi,1./real(sigmat.*sigmao))
linestyle k k-- 2k 2k-- k-- k-. 4D
xlim([0 2.5]),ylim([0 2])
vlines(sqrt([1:6]/2)/pi,'D:')
fixlabels(-1)

hlines(1./0.55,'D')
title('Wavelet Inverse Heisenberg Area')
xlabel('Duration $P_{\beta,\gamma}/\pi$')
ylabel('Inverse Area $1/A_{\beta,\gamma}$')
text(0.07,0.1,'(a)')
legend('$\gamma=1$','$\gamma=2$','$\gamma=3$','$\gamma=4$','$\gamma=5$','$\gamma=6$','Morlet','location','southeast');

set(gcf,'paperposition',[1 1 5 4])
fontsize 12 12 12 12

if strcmpi(str,'print')
     print -depsc wavelets_heisenberg_area.eps
end

figure
plot(p./pi,proj_morse.^2),xlim([0 2.5]),ylim([.92 1])
hold on, plot(p_morlet./pi,proj_morlet.^2)
linestyle k k-- 2k 2k-- k-- k-. 4D

title('Wavelet Envelope Versus Gaussian')
xlabel('Duration $P_{\beta,\gamma}/\pi$')
ylabel('Squared Inner Product')
ytick([.92:.01:1]),fixlabels([-1 -2])

text(0.1,0.925,'(b)')

legend('$\gamma=1$','$\gamma=2$','$\gamma=3$','$\gamma=4$','$\gamma=5$','$\gamma=6$','Morlet','location','southeast');

set(gcf,'paperposition',[1 1 5 4])
fontsize 12 12 12 12

if strcmpi(str,'print')
     print -depsc wavelets_versus_gaussian.eps
end
%\********************************************************************


%END of jlab_makefigs_superfamily
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[varargout]=jlab_makefigs_multivariate(str)
%Note: the last panel of Figure 1 had incorrect error ellipses for the top 
%two floats in the published version of the paper, due to an incorrect 
%calculation of the deviation vectors.  These estimates are actually 
%quite biased because I am using long wavelets to compensate for the
%fact that these two signals are near the Nyquist.  The signals do not
%appear to be well resolved, try 'cellplot(kap(25:26))' and note the strong
%oscillations in the estimated amplitude.

%/******************************************************************
load ebasnfloats
use ebasnfloats

len=cellength(lat);
index=find(len>200);
lato=30;
lono=-25;

vindex(id,num,lat,lon,p,t,index,1);


clear ga be fs p
for i=1:length(lat)
    if i==25||i==26
        %Special treatment for 36 (148000, now 25) and 37 (149000, now 26), 
        %       which are not getting ridges otherwise
        ga(i)=3;
        be(i)=8;
        fs{i}=morsespace(ga(i),be(i),{0.05,pi},2*pi/100,8);
    else
        ga(i)=3;
        be(i)=3;
        fs{i}=morsespace(ga(i),be(i),{0.05,2*pi/3},2*pi/100,8);
    end
    ppsi(i)=morseprops(ga(i),be(i));
end


%Compute wavelet transforms using generalized Morse wavelets
clear wx wy cx ir jr wxr wyr fxr fyr
for i=1:length(lat)
    i
    cx{i}=fillbad(latlon2xy(lat{i},lon{i},lato,lono));   
    [wx{i},wy{i}]=wavetrans(real(cx{i}),imag(cx{i}),{ga(i),be(i),fs{i},'bandpass'},'mirror');
    dt=num{i}(2)-num{i}(1);
    [wxr{i},wyr{i},ir{i},jr{i},fr{i},er,br,cr]=ridgewalk(dt,wx{i},wy{i},fs{i},ppsi(i),3); 
    
    if length(ir{i})>0
        bxr{i}=br(:,1);byr{i}=br(:,2);
        cxr{i}=cr(:,1);cyr{i}=cr(:,2);
    else
        [wxr{i},wyr{i},fxr{i},fyr{i},bxr{i},byr{i},cxr{i},cyr{i}]=vempty;
    end
end


clear xhatr yhatr xres lathat lonhat latres lonres kap lam the phi
clear fxhat fyhat cxhat cyhat fbar 

for i=1:length(lat)
    %Small overlap of two weak ridges for i=9 and i=23 
    [xhatr{i},yhatr{i},fbar{i},cxhat{i},cyhat{i}]=...
         ridgemap(length(cx{i}),wxr{i},wyr{i},fr{i},cxr{i},cyr{i},ir{i},'collapse');
    xhatmag=sqrt(squared(xhatr{i})+squared(yhatr{i}));
    cxhat{i}=frac(1,2)*squared(frac(ppsi(i),fbar{i})).*cxhat{i}.*xhatmag;
    cyhat{i}=frac(1,2)*squared(frac(ppsi(i),fbar{i})).*cyhat{i}.*xhatmag;
    [kap{i},lam{i},the{i},phi{i}]=ellparams(xhatr{i},yhatr{i});  
    xhat{i}=real(xhatr{i})+sqrt(-1)*real(yhatr{i});
    xres{i}=cx{i}-vswap(xhat{i},nan,0);
    [kaperr{i},lamerr{i},theerr{i},phierr{i}]=ellparams(cxhat{i},cyhat{i});  
    [latres{i},lonres{i}]=xy2latlon(xres{i},lato,lono);
end

clear cv cvres cvhat
for i=1:length(lat)
    cv{i}=latlon2uv(num{i},lat{i},lon{i});
    cv{i}=fillbad(cv{i});
    cvres{i}=latlon2uv(num{i},latres{i},lonres{i});
    cvres{i}=fillbad(cvres{i});
    cvhat{i}=cv{i}-cvres{i};
end
jj=find(cellength(ir)>0);
%\******************************************************************

%/******************************************************************
figure
ax=[-30.5 -19.5 18 36.5];

subplot(141),h=cellplot(lon,lat);linestyle -h h D
hold on,h=cellplot(lon(jj),lat(jj));linestyle -h h k
h1=plot(lon{24},lat{24});linestyle -h h1 5C
h1=plot(lon{24},lat{24});linestyle -h h1 1k


ylabel('Latitude')
xlabel('West Longitude'),latratio(30),axis(ax)
title('Float Trajectories')

subplot(142)
xlabel('West Longitude'),title('Signal Ellipses')

for i=1:length(lat)
    if anyany(~isnan(kap{i}))
        dt=round(frac(2*pi,fbar{i}));
        clear index
        index(1)=dt(find(~isnan(dt),1,'first'));

        while index(end)<length(kap{i})-dt(find(~isnan(dt),1,'last'))
            index(end+1)=index(end)+dt(index(end));
        end
        index=nonnan(index);
        index=index(~isnan(kap{i}(index)));
        if ~isempty(index)
            ar=[frac(360,2*pi*radearth) frac(360,2*pi*radearth*cosd(30))];
            h=ellipseplot(2*kap{i},lam{i},the{i},lonres{i}+sqrt(-1)*latres{i},ar,'index',index);hold on
        end
    end
end
latratio(30),axis(ax)
linestyle k D

subplot(143),h=cellplot(lonres,latres);linestyle -h h D
hold on,h=cellplot(lonres(jj),latres(jj));linestyle -h h k
xlabel('West Longitude'),latratio(30),axis(ax),title('Residuals')
h1=plot(lonres{24},latres{24});linestyle -h h1 5C
h1=plot(lonres{24},latres{24});linestyle -h h1 1k

subplot(144)
xlabel('West Longitude'),title('Bias Ellipses')

for i=1:length(lat)
    if anyany(~isnan(kaperr{i}))
        dt=round(frac(2*pi,fbar{i}));
        clear index
        index(1)=dt(find(~isnan(dt),1,'first'));

        while index(end)<length(kaperr{i})-dt(find(~isnan(dt),1,'last'))
            index(end+1)=index(end)+dt(index(end));
        end
        index=nonnan(index);
        index=index(~isnan(kaperr{i}(index)));
        if ~isempty(index)
            ar=[frac(360,2*pi*radearth*cosd(30)) frac(360,2*pi*radearth)];
            h=ellipseplot(2*kaperr{i},lamerr{i},theerr{i},lonres{i}+sqrt(-1)*latres{i},ar,'index',index);hold on
        end
    end
end
latratio(30),axis(ax)
linestyle k D


for i=1:4
    subplot(1,4,i)
    xtick([-30:2.5:-20])
    xtl=get(gca,'xticklabel');
    set(gca,'xticklabel',xtl(:,2:end))
    %fixlabels([-1,0])
    boxon
end
letterlabels(4)
packfig(1,4,'columns')


orient landscape
fontsize 12 10 10 10
set(gcf,'paperposition',[1 1 11 7])

if strcmpi(str,'print')
   print -depsc ebasn-trajectories.eps
end
%\******************************************************************


%/******************************************************************
figure
ci=(0:5:65);
ii=24;
numo=datenum(1986,1,1)-1;
[h,hl]=wavespecplot(num{ii}-numo,cv{ii},2*pi./fs{ii},sqrt(abs(wx{ii}).^2+abs(wy{ii}).^2),1,ci);
linestyle -h hl k k--
axes(h(1)),ylim([-18 18]),ylabel('Current Velocity (cm/s)'),title('Multivariate Ridge Method Example')
text(-90,15,'(a)')

axes(h(2)),caxis([0 40]),colormap gray,flipmap,ylim([3.6 60]),hold on
plot(num{ii}-numo,2*pi./fbar{ii},'w','linewidth',4)
plot(num{ii}-numo,2*pi./fbar{ii},'k','linewidth',2)

xlabel('Day of Year 1986'),ylabel('Period in Days')
set(gca,'ytick',2.^(2:.5:5.5))
set(gca,'yticklabel',[' 4';'  ';' 8';'  ';'16';'  ';'32';'  '])
inticks
text(-90,4.5,'(b)')


orient landscape
fontsize 16 14 14 14
set(gcf,'paperposition',[1 1 9 5])
if strcmpi(str,'print')
    print -deps ebasn-example.eps
end
%\******************************************************************


%END of jlab_makefigs_multivariate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[varargout]=jlab_makefigs_vortex(str)

disp('Making figures for Lilly, Scott, and Olhede (2011).')
disp('This may take a while.')

load vortex
%/************************************************************************
use vortex.drifters

disp('Computing wavelet transforms using generalized Morse wavelets.')
ga=3;be=2; 
fs=morsespace(ga,be,{.2 pi},pi./300,8); 
%[psi,psif]=morsewave(length(x),1,ga,be,fs);
[wx,wy]=wavetrans(unwrap(x/L)*L,y,{ga,be,fs,'bandpass'},'mirror');

make vortex.drifters wx wy
clear ir jr kr xr yr fbar kappa lambda theta phi R V
for k=1:size(x,2)
    disp(['Wavelet ridge analysis for drifter #' int2str(k) ' of 101.'])
    %The physical cutoff is 400/1000 km = 0.4 km
    [xr{k},yr{k},ir{k},jr{k},fbar{k}]=ridgewalk(dt,wx(:,:,k),wy(:,:,k),fs,morseprops(ga,be),1);
  
    %Keep track of which number drifter we're on, for future reference
    kr{k}=k+0*ir{k};
    
    %This is the joint instantaneous frequency, see Lilly and Olhede (2010)
    %Calculate ellipse parameters from ridges
    [kappa{k},lambda{k},theta{k},phi{k}]=ellparams(xr{k},yr{k});
    R{k}=ellrad(kappa{k},lambda{k},phi{k});
    V{k}=ellvel(dt*24*3600,kappa{k},lambda{k},theta{k},phi{k},1e5);
end
make vortex.ridges ir jr kr xr yr fbar kappa lambda theta phi R V
%\************************************************************************
%vsize(ir,jr,kr,xr,yr,fbar,kappa)


%/************************************************************************
disp('Mapping ridges back into time series.')

xhat=nan*x;
yhat=nan*y;

clear fhat
for k=1:size(x,2)    
    %RIDGEMAP sums over all ridges and forms properties the same size as the original time series
    [xhat(:,k),yhat(:,k)]=ridgemap(length(x),xr{k},yr{k},ir{k},'collapse');
    
    %This is the estimate of the aggregate oscillatory signals
    xhat(:,k)=real(xhat(:,k));
    yhat(:,k)=real(yhat(:,k));
end

%Careful to form x-residual from unwrapped x, then re-wrap
xres=L*angle(rot((L*unwrap(x/L)-vswap(xhat,nan,0))/L));
yres=y-vswap(yhat,nan,0);

%Calculate residual value appropriate for each ridge
xresr=ir;
yresr=ir;
for k=1:size(x,2)
    xresr{k}=xres(ir{k}(~isnan(ir{k})),k);
    yresr{k}=yres(ir{k}(~isnan(ir{k})),k);
end

x(abs(x-vshift(x,1,1))>L*pi)=nan;
xres(abs(xres-vshift(xres,1,1))>L*pi)=nan;

make vortex.ridges ir jr kr xr yr fbar kappa lambda theta phi R V xresr yresr
make vortex.drifters x y xres yres xhat yhat
%\************************************************************************

%That's the end of the processing, on to the figure making
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%/************************************************************************

figure
use vortex.snapshot
pcolor(x,y,abs(q)),shading interp,axis equal,axis square,axis tight,caxis([1 9]),hold on
colormap(squared(flipud(colormap('gray'))))
use vortex.ridges
cell2col(ir,jr,kr,fbar,R,V,kappa,theta,lambda,xresr,yresr);
kk=68; %That's the time index into the snapshot
index=find(ir==(kk-1)*20+1); 
h=ellipseplot(kappa,lambda,theta,xresr+sqrt(-1)*yresr,[1 1],'index',index,'3w'); 
index=find(ir==(kk-1)*20+1&abs(lambda)>0.8223); 
h=ellipseplot(kappa,lambda,theta,xresr+sqrt(-1)*yresr,[1 1],'index',index,'1.5g'); 
index=find(ir==(kk-1)*20+1&abs(lambda)<0.8223&V>0); 
h=ellipseplot(kappa,lambda,theta,xresr+sqrt(-1)*yresr,[1 1],'index',index,'1.5r'); 
index=find(ir==(kk-1)*20+1&abs(lambda)<0.8223&V<0); 
h=ellipseplot(kappa,lambda,theta,xresr+sqrt(-1)*yresr,[1 1],'index',index,'1.5b'); 

    
use vortex.drifters
plot(x((kk-1)*20+1,:),y((kk-1)*20+1,:),'o','markerfacecolor','w','markeredgecolor','k','markersize',4)
index=find(~isnan(xhat((kk-1)*20+1,:)));
plot(x((kk-1)*20+1,index),y((kk-1)*20+1,index),'o','markerfacecolor','k','markeredgecolor','w','markersize',4)
%text(-550*2,550*2,'Day 240','color','w')
axis off,xlabel('Distance East (km)'),ylabel('Distance North (km)')
%title('Snapshot of Unstable Jet with Lagrangian Ellipses')
h=colorbar('South','xcolor','w','ycolor','w');
if verLessThan('matlab','8.4.0')
    axes(h)
    xlabel('Potential Vorticity Magnitude')
else
    h.Label.String='Potential Vorticity Magnitude';
    pos=h.Position;
    h.Position=[0.2 pos(2) 0.6 pos(4)];
end

orient portrait
fontsize 14 12 12 12

if strcmpi(str,'print')
   print -depsc vortex-snapshot.eps
end
%\************************************************************************
   
%/************************************************************************
%Wavelet transform, for supplementary materials

jj=80;  %This is the drifter I'm using 
use vortex.drifters
figure
h=wavespecplot(num,x(:,jj)+sqrt(-1)*y(:,jj),2*pi./fs.*dt,log10(abs(wx(:,:,jj)).^2+abs(wy(:,:,jj)).^2));
axes(h(1)),hold on,linestyle 2b 2r
ylim([-900 900]),xlim([0 250]),ytick([-750:250:750])
uvplot(num,xres(:,jj)+sqrt(-1)*yres(:,jj),'g')
uvplot(num,vswap(xhat(:,jj),nan,0)+sqrt(-1)*vswap(yhat(:,jj),nan,0));hlines(0,'k:')
title('Example of Wavelet Ridge Analysis'),ylabel('Position (km)')
axes(h(1))
vlines([80 125],'k:')
axes(h(2)),hold on,xlim([0 250]),caxis([-2.75 4.25]),ytick([1 10 100])

%Time index
use vortex.ridges
numr=nan*ir{jj};numr(~isnan(ir{jj}))=num(ir{jj}(~isnan(ir{jj})));
plot(col2mat(numr),2*pi./col2mat(fbar{jj}));
linestyle b g c m
vlines([80 125],'k:')

ylabel('Period (days)'),xlabel('Time (days)')
letterlabels(4)

hc=colorbar('eastoutside');
if verLessThan('matlab','8.4.0')
    axes(hc),ylabel('Signal Stength (Log10 km)')
else
    hc.Label.String='Signal Stength (Log10 km)';
    colormap jet
end
pos1=get(h(1),'position');pos2=get(h(2),'position');
set(h(1),'position',[pos1(1:2) pos2(3) pos1(4)])

orient landscape
fontsize 22 18 18 18
if strcmpi(str,'print')
   print -depsc vortex-transform.eps
end
%\************************************************************************




%/************************************************************************
use vortex.ridges
%Make the ridges into one long column
cell2col(ir,jr,kr,fbar,R,V,kappa,theta,lambda,xresr,yresr);

figure
xbin=[0:2:160];ybin=[-130:2:130];

subplot(1,2,1)
[matm,xmid,ymid]=twodhist(R,V,xbin,ybin);
pcolor(xmid,ymid,log10(matm)),shading flat
colormap jet,map=colormap;map(1,:)=[1 1 1];colormap(map);
caxis([0 3]),hold on
xlabel('Radius (km)'),ylabel('Velocity (cm/s)')

%text(45,140,'Distributions of Estimated Oscillations in Unstable Jet Trajectories')

subplot(1,2,2)
[matm,xmid,ymid]=twodmed(R,V,ecconv(abs(lambda),'lin2ecc').^2,xbin,ybin);
pcolor(xmid,ymid,matm),shading flat
caxis([0 1]),hold on


for i=1:2
    subplot(1,2,i)
    axis([0 145 -180 130]),xtick([0:25:125]),ytick([-125:25:125]),%hlines(-32.5,'k')
    xlabel('Radius (km)'),ylabel('Velocity (cm/s)')
        
    %Ro= 2 V/Rf   so V/R =Ro f/2 vm/rm = (1/2)*0.8*
    %2*sind(45)*7.292e-5*100*1000  = 4.12
    
    plot([0+sqrt(-1)*0;10.5+10.5*sqrt(-1)*maxmax(fs)*100*1000/24/3600/dt],'k','linewidth',2)
    plot([0+sqrt(-1)*0;10.5-10.5*sqrt(-1)*maxmax(fs)*100*1000/24/3600/dt],'k','linewidth',2)   
    plot([0+sqrt(-1)*0;150+150*sqrt(-1)*minmin(fs)*100*1000/24/3600/dt],'k','linewidth',2)
    plot([0+sqrt(-1)*0;150-150*sqrt(-1)*minmin(fs)*100*1000/24/3600/dt],'k','linewidth',2)
    
    h1=plot([0+sqrt(-1)*0;30+sqrt(-1)*4.16*30]);linestyle -h h1 2D
    h1=plot([0+sqrt(-1)*0;30-sqrt(-1)*4.16*30]);linestyle -h h1 2D
end
letterlabels(2)
ha=packfig(1,2,'columns');

for i=1:2
    axes(ha(i))
    h=colorbar('South');
    if verLessThan('matlab','8.4.0')
        axes(h)
        switch i
            case 1
                xlabel('Log10 Number of Ridge Points')
            case 2
                xlabel('Median Squared Eccentricity \epsilon^2')
                xtick([0:0.2:1])
        end
    else
        switch i
            case 1
                h.Label.String='Log10 Number of Ridge Points';
            case 2
                h.Label.String='Median Squared Eccentricity \epsilon^2';
                h.Ticks=[0:0.2:1];
        end
    end

    pos=get(h,'position');
    set(h,'position',[pos(1) pos(2)+0.01 pos(3) pos(4)/2]);
end
orient tall
fontsize 12 10 10 12
set(gcf,'paperposition',[0.5 1 8 6])
 
%Seriously that is not cool Matlab set(gcf,'renderer','opengl')
set(gcf,'renderer','zbuffer')
%set(gcf,'renderer','painters')
 

if strcmpi(str,'print')
   print -depsc vortex-distributions.eps
end
%\************************************************************************

%/************************************************************************
figure

use vortex.drifters
ha(1)=subplot(2,2,1);
plot(x+sqrt(-1)*y),axis([-1 1 -1 1]*pi*L),set(gca,'dataaspectratio',[1 1 1])
linestyle default
ha(2)=subplot(2,2,2);
plot(xres+sqrt(-1)*yres),axis([-1 1 -1 1]*pi*L),set(gca,'dataaspectratio',[1 1 1])
linestyle default

use vortex.ridges
cell2col(ir,jr,kr,fbar,R,V,kappa,theta,lambda,xresr,yresr);
bool=~isinf(ir);
vindex(ir,jr,kr,fbar,R,V,kappa,theta,lambda,xresr,yresr,bool,1);
%Make the ridges into one long column

%0.8223=ecconv(0.95,'ecc2lin');l . %FYI

ha(3)=subplot(2,2,3);
index=periodindex(dt,fbar,3);
ii=find(lambda(index)<0.8223&V(index)<0);
h=ellipseplot(kappa,lambda,theta,xresr+sqrt(-1)*yresr,[1 1],'index',index(ii),'b');hold on;
ii=find(lambda(index)<0.8223&V(index)>0);
h=ellipseplot(kappa,lambda,theta,xresr+sqrt(-1)*yresr,[1 1],'index',index(ii),'r');hold on;
index=periodindex(dt,fbar,1/2)';
ii=find(lambda(index)>0.8223);
h=ellipseplot(kappa,lambda,theta,xresr+sqrt(-1)*yresr,[1 1],'index',index(ii),'2w');hold on
h=ellipseplot(kappa,lambda,theta,xresr+sqrt(-1)*yresr,[1 1],'index',index(ii),'g');hold on
axis([-1 1 -1 1]*pi*L)

ha(4)=subplot(2,2,4);
index=periodindex(dt,fbar,3)';
ii=find(lambda(index)<0.8223&V(index)<0);
h=ellipseplot(kappa,lambda,theta,xresr+sqrt(-1)*yresr,[1 1],'index',index(ii));hold on;
linecolor(h,log10(2*pi./fbar(index(ii))),.275,1.65)
ii=find(lambda(index)<0.8223&V(index)>0);
h=ellipseplot(kappa,lambda,theta,xresr+sqrt(-1)*yresr,[1 1],'index',index(ii));hold on;
linecolor(h,log10(2*pi./fbar(index(ii))),.275,1.65)
index=periodindex(dt,fbar,1/2)';
ii=find(lambda(index)>0.8223);
h=ellipseplot(kappa,lambda,theta,xresr+sqrt(-1)*yresr,[1 1],'index',index(ii),'2w');hold on
h=ellipseplot(kappa,lambda,theta,xresr+sqrt(-1)*yresr,[1 1],'index',index(ii));hold on
linecolor(h,log10(2*pi./fbar(index(ii))),.275,1.65)
axis([-1 1 -1 1]*pi*L)

for i=1:4
    set(ha(i),'xtick',[-1000:500:1000])
    set(ha(i),'ytick',[-1000:500:1000])
end
letterlabels(1);
packfig(2,2)

axes(ha(4))
h=colorbar('south');
if verLessThan('matlab','8.4.0')
    axes(h)
    xtick([ceil([log10([2.5 5 10 20 40 80])-.25]*64/1.4)]/64)
    set(gca,'xticklabel',['2.5';'5  ';'10 ';'20 ';'40 ';'80 '])
    xlabel('Log10 Oscillation Period')
else
    h.Label.String='Log10 Oscillation Period';
    h.Ticks=[ceil([log10([2.5 5 10 20 40 80])-.25]*64/1.4)]/64;
    h.TickLabels={'2.5','5  ','10 ','20 ','40 ','80 '};
end
pos=get(h,'position');
set(h,'position',[pos(1) pos(2) pos(3) pos(4)/2]);

orient tall
fontsize 12 10 10 10
set(gcf,'paperposition',[1 1 8.3 8.5])


if strcmpi(str,'print')
   print -depsc vortex-decomposition.eps
end
%\************************************************************************


if false  %To make the movie, change this line to 'if true'
    %/**********************************************************
    figure
    if ~exist('vortexpv.mat')==2
        disp('Sorry, MAKEFIGS_VORTEX can''t find the file VORTEXPV.MAT')
        disp('You need to download this via anonymous ftp from ftp.nwra.com/outgoing/lilly/vortex/.')
    end
    
    load vortexpv
    
    use vortex.ridges
    %Make the ridges into one long column
    cell2col(ir,jr,kr,fbar,R,V,kappa,theta,lambda,xresr,yresr);
    
    n=0;%
    
    cd vortexmovie
    colors=[0 0 1;0 1/2 0;1 0 0; 0 .75 .75; .75 0 .75;.75 .75 0];
    for k=[(1:101) 101-1:-1:2]
        clf
        use vortexpv
        pcolor(x,y,abs(q(:,:,k))),shading interp,axis equal,axis square,axis tight,caxis([1 9]),hold on
        colormap(squared(flipud(colormap('gray'))))
        index=find(ir==(k-1)*20+1);
        h=ellipseplot(kappa,lambda,theta,xresr+sqrt(-1)*yresr,[1 1],'index',index);linestyle -h h 3w
        h=ellipseplot(kappa,lambda,theta,xresr+sqrt(-1)*yresr,[1 1],'index',index);%linestyle -h h colors2
        
        for i=1:length(index)
            set(h(i),'color',colors(mod(kr(index(i))-1,6)+1,:),'linewidth',2)
        end
        
        index=find(ir==(k-1)*20+1&kr==80);   %... the one I use later
        h=ellipseplot(kappa,lambda,theta,xresr+sqrt(-1)*yresr,[1 1],'index',index);
        if ~isempty(h), linestyle -h h 2k, end
        
        use vortex.drifters
        plot(x((k-1)*20+1,:),y((k-1)*20+1,:),'o','markerfacecolor','w','markeredgecolor','k','markersize',4)
        index=find(~isnan(xhat((k-1)*20+1,:)));
        plot(x((k-1)*20+1,index),y((k-1)*20+1,index),'o','markerfacecolor','k','markeredgecolor','w','markersize',4)
        plot(x((k-1)*20+1,80),y((k-1)*20+1,80),'o','markerfacecolor','m','markeredgecolor','k','markersize',4)
        
        use vortexpv
        text(-550*2,550*2,['Day ' int2str(floor(num(k))) '.' int2str(floor(10*(num(k)-floor(num(k)))))],'color','w')
        text(-550*2,-550*2,'Vortex extraction from Lagrangian trajectories. Lilly, Scott, & Olhede (2011).','color','w')
        %axis off,
        xlabel('Distance East (km)'),ylabel('Distance North (km)')
        title('Movie of Unstable Jet with Lagrangian Ellipses')
        h=gca;
        axis([min(x) max(x) min(y) max(y)])
        %hc=colorbar('South','xcolor','w','ycolor','w');
        %hc.Label.String=('Potential Vorticity Magnitude');
        axis off,title('')
        set(gcf,'paperposition',[0 0 5 5])
        set(gcf,'papersize',[5 5])
        set(gca,'position',[0 0 1 1])
        print('-djpeg',['movieframe' int2str(n)])
        n=n+1;
    end
end
%\**********************************************************


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if false
%To double-check parameter values
load vortexpv
use vortexpv
figure
plot(y,detrend(q(:,:,1))),hold on
plot(y,-vdiff(squared(cos(y./40*pi/2)),1)*130,'g')
hold on
vlines([-40 40])


load vortex
use  vortex.drifters
cv=vdiff(unwrap(x/L)*L,1)./dt+sqrt(-1)*vdiff(y,1)./dt;
figure,plot(abs(cv))

zeta=pi*2.08/(80*1000)
2*2*omega*sind(45)./radearth*80/zeta
1/10/piS
%ok
%1.6*L*1000/(T*3600*24)=2.08
end


%END of jlab_makefigs_vortex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[varargout]=jlab_makefigs_trivariate(str)
 
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
mat=matmult(jmat3(alpha,3),jmat3(beta,1),1);
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
print -depsc trivariate_sphere.eps
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
print -depsc trivariate_seismic.eps
%\************************************************************************


%END of jlab_makefigs_trivariate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[]=jlab_makefigs_analytic(str)
%/*************************************************
figure
ga=[3 3 3 3 nan 1/3.5 1 6 18 nan];be=[1.5 3 7 81 9./ga(5:end)];
[p,skew,kurt]=morseprops(ga,be);

t=1:5000;t=t-mean(t);
[psi,psif]=vzeros(5000,5,'nan');
for i=1:length(be)
    [psi(:,i),psif(:,i)]=morsewave(5000,1,ga(i),be(i),2*pi/40);
end

textstr{1}='(1.5,3)';
textstr{2}='(3,3)';
textstr{3}='(7,3)';
textstr{4}='(81,3)';
textstr{5}='';
textstr{6}='(31.5,0.29)';
textstr{7}='(9,1)';
textstr{8}='(1.5,6)';
textstr{9}='(0.5,18)';
textstr{10}='';

clear ha
for i=1:10
    ha(i)=subplot(2,5,i);
    if i~=5&&i~=10
        psi(:,i)=psi(:,i)./maxmax(abs(psi(:,i)));
        uvplot(t./p(i).*(2*pi/40),psi(:,i)),hold on,plot(t./p(i).*(2*pi/40),abs(psi(:,i))),linestyle k E-- 2k
        ylim([-1 1.05]),xlim([-4.5 4.5]),
        if i==3
            %text(2,1.5,'Examples of Generalized Morse Wavelets','fontsize',14)
            title('Examples of Generalized Morse Wavelets')
        %    pos=get(get(gca,'title'),'position');
        %    set(get(gca,'title'),'position',pos-[0 .01 0])
        end
        text(-4.2,-.9,textstr{i})
        fixlabels([0 -1]),hlines(0,'k:'),
        if i<4,hlines([-1 1.05],'k'),end
        if i==1,vlines(-4.5,'k'),end
        if i==3,vlines(4.5,'k'),end
        set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
    else
        ylim([0 2.02]),xlim([0 2.5]),hold on
    end
   %title(textstr{i})
end

f=(0:5000-1)./126;
h=subplot(2,5,5);plot(f,psif(:,1:4)),ylim([0 2.02]),xlim([0 2.5]), linestyle E-- 2k k 2E--
ytick(0:.5:2),fixlabels([-1 -1])
boxoff,set(gca,'xtick',[],'ytick',[]),vlines(1,'k:'),
%title('Frequency Domain')
legend('a','b','c','d')

h=subplot(2,5,10);plot(f,psif(:,6:10)),ylim([0 2.02]),xlim([0 2.5]), linestyle E-- 2k k 2E--
ytick(0:.5:2),xlabel('Frequency / Peak Frequency'),fixlabels([-1 -1])
boxoff,set(gca,'xtick',[],'ytick',[]),vlines(1,'k:')
%title('Frequency Domain')
legend('f','g','h','i')

letterlabels(ha,1);

packfig(2,5)

orient landscape
fontsize jpofigure
set(gcf,'paperposition',[1 1 9 4.5])

if strcmpi(str,'print')
   %
   print -depsc analytic-morsies.eps
end
%\******************************


%/************************************************************
N=20;
ga1=[(1:1:21)];
be1=[(1:1:21)];
[ga,be]=meshgrid(ga1,be1);
[p,skew,kurt]=morseprops(ga,be);
dcell=morsederiv(N,ga,be);

x=zeros([size(ga) N]);
for n=1:length(dcell);
   if iseven(n)
       fact=n/2;
   else
       fact=(n-1)/2;
   end
   x(:,:,n)=dcell{n}./((p.^2).^fact)./factorial(n);
end

figure
for i=1:size(x,1)
    for j=1:size(x,2)
           if ga1(j)<6
                plot(2:size(x,3),squeeze(abs(x(i,j,2:end))),'ko','markersize',5,'markerfacecolor','k'),hold on,ylog
           else 
                plot(2:size(x,3),squeeze(abs(x(i,j,2:end))),'k.','color',[.6 .6 .6]),hold on,ylog
           end
           ylim([10^-6 5]),xlim([1.8 20.5]),xtick([2 3 4 5 10 15 20]),xlog
    end
end

for i=1:size(x,1)
    for j=1:size(x,2)  
           if ga1(j)==3
%                 h=plot(4:size(x,3),squeeze(abs(x(i,j,4:end)))); linestyle -h h 2w
                 h=plot(4:size(x,3),squeeze(abs(x(i,j,4:end)))); linestyle -h h F

           end
    end
end

hlines(1,'k:'),title('Morse Wavelet Decay with $\beta>1$, $1<\gamma< 6$')
xlabel('Derivative Number at Peak Frequency'),ylabel('Normalized Magnitude')

fontsize jpofigure
set(gcf,'paperposition',[1 1 3.5 3.5])

if strcmpi(str,'print')
   print -depsc analytic-decay.eps
end
%\********************************************************  



%/*************************************************************************
load ebasnfloats
use ebasnfloats
num=num{33};lat=lat{33};lon=lon{33};p=p{33};t=t{33};
vindex(num,lat,lon,1:549,1);

num=num-datenum(1986,1,1)+1;

%vindex(num,lat,lon,1:400,1);
cx=latlon2xy(lat,lon,vmean(lat(:),1),vmean(lon(:),1));
cv=latlon2uv(num,lat,lon);

dt=(num(2)-num(1));
mlat=vmean(lat(:),1);

fmax=abs(corfreq(mlat))*frac(dt*24,2); %One cycle per 2 inertial periods = 1 cycle per 2.6 day
fmin=abs(corfreq(mlat))*frac(dt*24,40);%One cycle per 40 inertial periods = 1 cycle per 53 days
fs=morsespace(3,3,fmax,fmin,8);

ga=[3 3 3];
be=[1.5 3 7];
p=morseprops(ga, be);
c=squared(2./squared(p));

clear ir jr xr fr ar br cr xra fra bra cra wx
    
for i=1:3
      wx(:,:,i)=wavetrans(real(cx),{1,ga(i),be(i),fs,'bandpass'},'mirror');
    [xr{i},ir{i},jr{i},fr{i},br{i},cr{i}]=ridgewalk(dt,wx(:,:,i),fs,2*morseprops(ga(i),be(i)));  
    [xra(:,i),fra(:,i),bra(:,i),cra(:,i)]=ridgemap(length(cx),xr{i},fr{i},br{i},cr{i},ir{i},'collapse');
end

if 0
    x=[real(xra(:,1)) imag(xra(:,1))];
    for i=1:8
        x=[x;flipud(x)];
    end
    audiowrite('eddysound.wav',x,8192/2);
end
    
figure
M=4;N=3;
ci=(0:5:65)/2;
  
for i=1:size(xra,2)
    subplot(M,N,i)
    plot(num,vdiff([  real(cx-xra(:,i)) real(xra(:,i)) real(cx)],1))
    yoffset 25,
    linestyle  1.5k  G k
    axis([-95 num(end) -10 70 ])
    if i==1,ylabel('Velocity (cm/s)'),end
    title(['Ridge Analysis with $\gamma=$' int2str(ga(i)) ', $\beta=$' num2str(be(i))])
    
    subplot(M,N,i+N)
    contourf(num,2*pi./fs,abs(wx(:,:,i)'),ci),caxis([1 25]),colormap gray,flipmap,flipy,ylog,hold on,nocontours
    plot(num(ir{i}(~isnan(ir{i}))),2*pi./fs(jr{i}(~isnan(ir{i}))),'w','linewidth',4)
    plot(num(ir{i}(~isnan(ir{i}))),2*pi./fs(jr{i}(~isnan(ir{i}))),'k','linewidth',2)
    ytick([2 4 8 16 32])
    if i==1,ylabel('Period in Days'),end
    
    subplot(M,N,i+2*N)
    plot(num,fra(:,i),'k'),hold on,plot(num,bra(:,i),'k--')
    axis([-95  num(end) -.5 1.9]),ytick([0:.5:1.5])
    if i==1,ylabel('$\omega(t)$ $\&$ $\upsilon(t)$'),fixlabels([0 -1]),end
    hlines(0,'k:')
    
    subplot(M,N,i+3*N)
    uvplot(num,frac(cra(:,i),fra(:,i).^2)),hold on,ytick([-.2:.1:.2]),fixlabels([0 -1])
    plot(num,frac(bra(:,i),fra(:,i)).^2),%plot(num,[abs(frac(cra(:,i),fra(:,i).^2)) -abs(frac(cra(:,i),fra(:,i).^2))])
    linestyle k G 2k     
    axis([-95  num(end) -.225 .225]),hlines([-1 1]*c(i),'k:')
    if i==1,ylabel('$\rho_2(t)$ $\&$ $\rho_1^2(t)$'),end  
    
    xlabel('Day of Year 1986')
end

for i=1:M*N
    h(i)=subplot(M,N,i);%text(,['(' setstr(i+96) ')' ])
    axes(h(i))
    xtick([-100:100:500])
    set(gca,'xticklabelmode','auto')
end
letterlabels(h,1);


packfig(4,3)
set(gcf,'paperposition',[1 1 9 4.5])
fontsize jpofigure
if strcmpi(str,'print')
   print -depsc analytic-transforms.eps
end

% %Reconstructing mean and median estimates
% for i=1:3
%      wx(:,:,i)=wavetrans(real(xra(:,i)),{1,ga(i),be(i),fs,'bandpass'},'mirror');
%     [ir{i},jr{i},xr{i},fr{i},br{i},cr{i}]=ridgewalk(dt,wx(:,:,i),fs,2*morseprops(ga(i),be(i)));  
%     [xra2(:,i),fra(:,i),bra(:,i),cra(:,i)]=ridgemap(length(cx),xr{i},fr{i},br{i},cr{i},ir{i},'collapse');
% end
% 
% for i=1:3
%     vmean(squared(xra2(:,i)-xra(:,i))./squared(xra(:,i)),1)
%     vmedian(squared(xra2(:,i)-xra(:,i))./squared(xra(:,i)),1)
% end

%END of jlab_makefigs_analytic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=jlab_makefigs_bandwidth(str)
%/****************************************************************
%Schematic of nonzero bandwidth

%This figure is in makefigs_ellband
makefigs_ellband


if strcmpi(str,'print')
    orient landscape
    set(gcf,'paperposition',[1 8 1 4])
    fontsize 14 12 12 12
    print -dpsc stability_schematic.eps
end
%\****************************************************************


%/*************************************************
load ebasnfloats

use ebasnfloats
num=num{33};lat=lat{33};lon=lon{33};p=p{33};t=t{33};
vindex(num,lat,lon,1:549,1);

num=num-datenum(1986,1,1)+1;

cx=latlon2xy(lat,lon,vmean(lat(:),1),vmean(lon(:),1));
cv=latlon2uv(num,lat,lon);

ga=3;be=3;

dt=(num(2)-num(1));
mlat=vmean(lat(:),1);

 
fmax=abs(corfreq(mlat))*frac(dt*24,2); %One cycle per 2 inertial periods = 1 cycle per 2.6 day
fmin=abs(corfreq(mlat))*frac(dt*24,40);%One cycle per 40 inertial periods = 1 cycle per 53 days
fs=morsespace(ga,be,{0.2,fmax},fmin,8);


%Compute wavelet transforms using generalized Morse wavelets
[wx,wy]=wavetrans(real(cx),imag(cx),{ga,be,fs,'bandpass'},'mirror');
%[wp,wn]=vectmult(tmat,wx,wy);

%Form ridges of component time series
[wrx,wry,ir,~,fr,er,br]=ridgewalk(dt,wx,wy,fs,3,2);   


%Map into time series locations
[wrx,wry,fr]=ridgemap(length(cx),wrx,wry,fr,ir);

%Convert xy transforms to ellipse formsdf
[kappa,lambda,theta,phi]=ellparams(wrx,wry);

%Other ellipse properties from xy transforms
rm=ellrad(kappa,lambda);
vm=ellvel(24*3600,kappa,lambda,theta,phi,1e5);

%Other frequencies
[fphi,fth]=vdiff(2*pi*dt,phi,theta,1,'nans');

z=ellsig(kappa,lambda,theta,phi);
fr=fr./(2*pi);

[upa,upb,upc]=ellband(dt,kappa,lambda,theta,phi);
%\*************************************************

%Bandwidth from ridgewalk not agreeing with ellband,
%also figure looks sightly different for bandwith
%up=sqrt(squared(wrx)).*br(:,1))+squared(wry.*br(:,2)))./sqrt(squared(wrx)+squared(wry));
%plot(up),hold on
%plot(sqrt(squared(upa)+squared(upb)+squared(upc)))


%/*************************************************
figure
subplot(511)
uvplot(num,z),hold on, plot(num,[kappa -kappa]),linestyle k k-- 2k 2k
axis tight,ylim([-30 30]),hlines(0,'k:'),ytick(-20:10:30),xlim([-100 440])
ylabel('Displacement')
title('Vortex Motion as a Modulated Elliptical Signal')

subplot(512)
plot(num,-lambda),
axis tight,ylim([0 .55]),boxon,ytick(0:.10:.5),fixlabels([0 -2]),xlim([-100 440])
linestyle k,hlines(0,'k:');
ylabel('Linearity Magnitude')

subplot(513)
plot(num,[fr fphi fth]),linestyle 2k k k-- 2k
axis tight,ylim([-0.1 0.35]),hlines(0,'k:'),ytick(-.1:.1:.3),xlim([-100 440])
ylabel('Cyclic Frequency')
fixlabels([0 -1])

subplot(514)
plot(num,vfilt(abs([upa upb upc]),20,'mirror'))
axis tight,ylim([0 .09]),boxon,ytick(0:.02:.09),fixlabels([0 -2]),xlim([-100 440])
linestyle 2k k k-- k-.
hlines(0,'k:');
ylabel('Bandwidth')


subplot(515)
plot(num,vfilt(abs([upa upb upc]),20,'mirror')./[fr fr fr]./(2*pi))
axis tight,ylim([0 0.17]),boxon,ytick(0:.03:.15),fixlabels([0 -2]),xlim([-100 440])
linestyle 2k k k-- k-.
hlines(0,'k:');
ylabel('Relative Bandwidth')
xlabel('Day of Year 1986')
letterlabels(2)

packfig(5,1,'rows')


if strcmpi(str,'print')
    fontsize 14 12 12 12
    orient tall
    print -deps bandwidth_signal.eps
end
%\*************************************************


%/*************************************************
cxr=cx-z;

cx_nan=cx+0*vswap(z,0,nan);
cxr_nan=cxr+0*vswap(z,0,nan);

figure,
subplot(131);h=plot(cx,'k');hold on
axis equal,axis([-250 280  -280 350]),
title('Bivariate Data')
xtick(-1200:100:400),ytick(-1200:100:600)
ylabel('Displacement North (km)')
xlabel('Displacement East (km)')
plot(cx(1,:),'k^','markersize',5, 'MarkerFaceColor','k')

subplot(132),
%h=plot(cxr); linestyle -h h C
clear index
index(1)=6*4;
while index(end)<length(kappa)-6*2
    index=[index;round(index(end)+(1./2)./fr(index(end)))];
end
%index=(6*4:12:length(a)-6*4);
h=ellipseplot(kappa,lambda,theta,cxr,'index',index);hold on,linestyle D K
%linestyle C E G I K
%linestyle(h(~isnan(h)),'k')

axis equal,axis([-250 280 -280 350]),
xtick(-1200:100:400),ytick(-1200:100:600)
title('Estimated Local Ellipses')
xlabel('Displacement East (km)')
plot(cx(1,:),'k^','markersize',5, 'MarkerFaceColor','k')


subplot(133),
h=plot(cxr,'k');hold on
axis equal,axis([-250 280  -280 350]),
title('Residual')
xtick(-1200:100:400),ytick(-1200:100:600)
xlabel('Displacement East (km)')
plot(cx(1,:),'k^','markersize',5, 'MarkerFaceColor','k')

letterlabels(4)
packfig(1,3,'columns')


if strcmpi(str,'print')
    fontsize 10 8 8 8
    set(gcf,'paperposition',[1 1 6 4])
    print -depsc bandwidth_looper.eps
end
%\*************************************************


%/********************************************************
%Ellipsesketch

a=3.5;
b=2;

phi=pi/3;
th=pi/3;

[k,l]=ab2kl(a,b);
figure
ellipseplot(k,l,th,'npoints',64,'phase',phi),hold on,linestyle k
plot(rot(th+pi/2)*[0 1]*b,'k--')
plot(rot(th)*[0 1]*a,'k--')
plot(rot(th)*(a*cos(phi)+sqrt(-1)*b*sin(phi)),'ko','markersize',10)
xlabel('Displacement East')
ylabel('Displacement North')

title('Sketch of Ellipse')
axis equal
axis([-3.4 3.4 -3.4 3.4])
vlines(0,':'),hlines(0,':')
ytick(1),xtick(1)

xi=(0:.1:th);
plot(1.25*rot(xi),'k');

xi=(th:.01:pi/1.71);
plot(1.75*rot(xi),'k');

text(1.5,2.3,'a')
text(-1,0.8,'b')
text(0.3,2,'$\phi$')
text(1.3,0.5,'$\theta$')

xtick([-3:3]),ytick([-3:3])
%fixlabels(-1)
fontsize jpofigure
set(gcf,'paperposition',[2 2 3.5 3.5])

if strcmpi(str,'print')
   print -deps ellipsesketch.eps
end
%!gv ellipsesketch.eps &  
%\********************************************************


%END of jlab_makefigs_bandwidth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[varargout]=jlab_makefigs_asilomar(str)
%/******************************************************************
load ebasnfloats
use ebasnfloats

len=cellength(lat);
index=find(len>200);
lato=30;
lono=-25;

vindex(id,num,lat,lon,p,t,index,1);


clear ga be fs p
for i=1:length(lat)
    if i==25||i==26
        %Special treatment for 36 (148000, now 25) and 37 (149000, now 26), 
        %       which are not getting ridges otherwise
        ga(i)=3;
        be(i)=8;
        fs{i}=morsespace(ga(i),be(i),{0.05,pi},2*pi/100,8);
    else
        ga(i)=3;
        be(i)=3;
        fs{i}=morsespace(ga(i),be(i),{0.05,2*pi/3},2*pi/100,8);
    end
    ppsi(i)=morseprops(ga(i),be(i));
end


%Compute wavelet transforms using generalized Morse wavelets
clear wx wy ir jr wxr wyr fxr fyr
for i=1:length(lat)
    i
    cx{i}=fillbad(latlon2xy(lat{i},lon{i},lato,lono));   
    [wx{i},wy{i}]=wavetrans(real(cx{i}),imag(cx{i}),{ga(i),be(i),fs{i},'bandpass'},'mirror');
    dt=num{i}(2)-num{i}(1);
    [wxr{i},wyr{i},ir{i},jr{i},fr{i}]=ridgewalk(dt,wx{i},wy{i},fs{i},ppsi(i),3); 
end


clear xhatr yhatr xres lathat lonhat latres lonres kap lam the phi
clear fxhat fyhat cxhat cyhat fbar 

for i=1:length(lat)
    %Small overlap of two weak ridges for i=9 and i=23 
    [xhatr{i},yhatr{i},fbar{i}]=...
         ridgemap(length(cx{i}),wxr{i},wyr{i},fr{i},ir{i},'collapse');
    xhatmag=sqrt(squared(xhatr{i})+squared(yhatr{i}));
    [kap{i},lam{i},the{i},phi{i}]=ellparams(xhatr{i},yhatr{i});  
    xhat{i}=real(xhatr{i})+sqrt(-1)*real(yhatr{i});
    xres{i}=cx{i}-vswap(xhat{i},nan,0);
    [latres{i},lonres{i}]=xy2latlon(xres{i},lato,lono);
end

clear cv cvres cvhat
for i=1:length(lat)
    cv{i}=latlon2uv(num{i},lat{i},lon{i});
    cv{i}=fillbad(cv{i});
    cvres{i}=latlon2uv(num{i},latres{i},lonres{i});
    cvres{i}=fillbad(cvres{i});
    cvhat{i}=cv{i}-cvres{i};
end
jj=find(cellength(ir)>0);
%\******************************************************************

%/******************************************************************
figure
ax=[-30.5 -19.5 18 36.5];

subplot(131),h=cellplot(lon,lat);linestyle -h h D
hold on,h=cellplot(lon(jj),lat(jj));linestyle -h h k
h1=plot(lon{24},lat{24});linestyle -h h1 5C
h1=plot(lon{24},lat{24});linestyle -h h1 1k


ylabel('Latitude')
xlabel('West Longitude'),latratio(30),axis(ax)
title('Float Trajectories')

subplot(132)
xlabel('West Longitude'),title('Signal Ellipses')

for i=1:length(lat)
    if anyany(~isnan(kap{i}))
        dt=round(frac(2*pi,fbar{i}));
        clear index
        index(1)=dt(find(~isnan(dt),1,'first'));

        while index(end)<length(kap{i})-dt(find(~isnan(dt),1,'last'))
            index(end+1)=index(end)+dt(index(end));
        end
        index=nonnan(index);
        index=index(~isnan(kap{i}(index)));
        if ~isempty(index)
            ar=[frac(360,2*pi*radearth) frac(360,2*pi*radearth*cosd(30))];
            h=ellipseplot(2*kap{i},lam{i},the{i},lonres{i}+sqrt(-1)*latres{i},ar,'index',index);hold on
        end
    end
end
latratio(30),axis(ax)
linestyle k D

subplot(133),h=cellplot(lonres,latres);linestyle -h h D
hold on,h=cellplot(lonres(jj),latres(jj));linestyle -h h k
xlabel('West Longitude'),latratio(30),axis(ax),title('Residuals')
h1=plot(lonres{24},latres{24});linestyle -h h1 5C
h1=plot(lonres{24},latres{24});linestyle -h h1 1k

for i=1:3
    subplot(1,3,i)
    xtick([-30:2.5:-20])
    xtl=get(gca,'xticklabel');
    set(gca,'xticklabel',xtl(:,2:end))
    %fixlabels([-1,0])
    boxon
end
letterlabels(4)
packfig(1,3,'columns')


orient landscape
fontsize 12 10 10 10
set(gcf,'paperposition',[1 1 11 7])

if strcmpi(str,'print')
   print -depsc ebasn-trajectories.eps
end
%\******************************************************************


%/*************************************************    
figure
ci=(0:5:65);
ii=24;
numo=datenum(1986,1,1)-1;
[h,hl]=wavespecplot(num{ii}-numo,cv{ii},2*pi./fs{ii},sqrt(abs(wx{ii}).^2+abs(wy{ii}).^2),1,ci);
linestyle -h hl k k--
axes(h(1)),ylim([-18 18]),ylabel('Current Speed (cm/s)'),title('Bivariate Ridge Method Example')
text(-90,15,'(a)')

axes(h(2)),caxis([0 40]),colormap gray,flipmap,ylim([3.6 60]),hold on
plot(num{ii}-numo,2*pi./fbar{ii},'w','linewidth',4)
plot(num{ii}-numo,2*pi./fbar{ii},'k','linewidth',2)

xlabel('Day of Year 1986'),ylabel('Period in Days')
set(gca,'ytick',2.^(2:.5:5.5))
set(gca,'yticklabel',[' 4';'  ';' 8';'  ';'16';'  ';'32';'  '])
inticks
text(-90,4.5,'(b)')


orient landscape
fontsize 12 10 10 10
set(gcf,'paperposition',[1 1 9 5])
if strcmpi(str,'print')
    print -deps ebasn-example.eps
end
%\*************************************************

%END of jlab_makefigs_asilomar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[varargout]=jlab_makefigs_morsies(str) 
%/************************************************************
%Linear chirp
N=2000;
t=(0:N-1)';t=t-mean(t);
x=cos(t.^2/5000).*exp(-(t./250).^2);
%plot(t,x)
 
om=vdiff(t.^2/5000,1);

fs=(1./(logspace(log10(10),log10(2000),50)'));
%fs=[-flipud(fs);fs];

psi1=morsewave(N,1,3,1,2*pi*fs,'energy');
psi2=morlwave(N,2*pi*1.5./2./pi,2*pi*fs,'energy');
wx1=wavetrans(x,psi1,'zeros');
wx2=wavetrans(x,psi2,'zeros');
wx1=wx1./maxmax(abs(wx1));
wx2=wx2./maxmax(abs(wx2));
figure,
subplot(3,1,1),plot(t,[x exp(-(t./250).^2)]),linestyle k 2k,xlim([-400 400]),hlines(0,'k:'),ylim([-1.2 1.2]),ytick(-1:.5:1),fixlabels([0 -1])
subplot(3,1,2),contourf(t,fs*100,(abs(wx1)'),(0:.1:1)),xlim([-400 400]),ylim([0 .1]/2*100), colormap gray,flipmap,nocontours
subplot(3,1,3),contourf(t,fs*100,(abs(wx2)'),(0:.1:1)),xlim([-400 400]),ylim([0 .1]/2*100), colormap gray,flipmap,nocontours

for i=2:3
    subplot(3,1,i)
    if i>1
        hold on,
        plot(t,100*om./2./pi,'w','linewidth',1)
        plot(t,100*om./2./pi,'k--','linewidth',1),
        plot(t,-100*om./2./pi,'w','linewidth',1)
        plot(t,-100*om./2./pi,'k--','linewidth',1),
        ytick(0:1:3),ylim([0 3.2])
        ylabel('Frequency $\omega/(2 \pi) \times 100$')
    end
    xtick(-400:200:400)
end
subplot(3,1,1),title('A Gaussian-Enveloped Chirp')

packfig(3,1,'rows')

fontsize 8 8 8 8
letterlabels(2)

orient tall
set(gcf,'paperposition',[1 1 2.75 6])
if strcmpi(str,'print')
    %
    print -deps morsies_linear_chirp.eps
end
%\************************************************************

%/************************************************************
figure

N=512*2*4;
fs=1/8/4/2;
t=1:N;t=t-mean(t);t=t.*fs;
clear psi psif 
P=zeros(4,4);
for i=1:4
    for j=0:3
        P(i,j+1)=sqrt(i*j);
        if j==0
           %This is due to the definition of the beta=0 case in morsefreq
           [psi(:,i,j+1) psif(:,i,j+1)]=morsewave(N,1,i,j,2*pi*fs*frac(morsefreq(i,0),morsefreq(i,1)), 'bandpass');
        else 
           [psi(:,i,j+1) psif(:,i,j+1)]=morsewave(N,1,i,j,2*pi*fs,'bandpass');
        end
        
    end
end
%P=P(:);
%P(1:4)=P(5:8);
P(:,1)=P(:,2);

for i=1:4
    for j=1:4
        subplot(4,5,i+(j-1)*5)
        x=[squeeze(real(psi(:,i,j))) squeeze(imag(psi(:,i,j))) squeeze(abs(psi(:,i,j)))]./max(abs(psi(:,i,j)));
        plot(t./(P(i,j)./pi),x),xlim([-2 2]),box on,ylim([-0.85 1.1])
        linestyle k k-- 1.5k
        hlines(0,'k:')
        set(gca,'xtick',[]), set(gca,'ytick',[])
        set(gca,'xcolor','w')
        set(gca,'ycolor','w')
        %xlim([-1.5 1.5]/2)
    end
end


for j=1:4
    subplot(4,5,5+(j-1)*5)
    plot(squeeze(psif(:,:,j)))
    xtick([0 65 210]),xlim([0 200]),ytick([-1 3]),%vlines(65,'k:')
    set(gca,'xticklabel',[]),set(gca,'yticklabel',[]),linestyle k k-- k-. k:
    boxoff
    set(gca,'ticklen',[0.025 0.025]*2)
end


subplot(4,5,1),title('Cauchy Family ($\gamma=1$)')
subplot(4,5,2),title('Gaussian Family ($\gamma=2$)')
subplot(4,5,3),title('Airy Family ($\gamma=3$)')
subplot(4,5,4),title('Hyper-Gaussian Family ($\gamma=4$)')
subplot(4,5,5),title('Frequency Domain')
subplot(4,5,1),text(-1.8,.9,'$\beta=0$')
subplot(4,5,6),text(-1.8,.9,'$\beta=1$')
subplot(4,5,11),text(-1.8,.9,'$\beta=2$')
subplot(4,5,16),text(-1.8,.9,'$\beta=3$')


packfig(4,5)
legend(gca,'1','2','3','4')
legend(gca,'Mod','Re','Im')

fontsize 10 10 10 10
orient landscape
set(gcf,'paperposition',[1 1 10 6])

if strcmpi(str,'print')
    print -deps morsie_families.eps
end
%\*********************************************************


%/********************************************************************
%The code for this figure is in MAKEFIGS_WIGDIST
h=gcf;
wigdist --f

figure(h.Number+1)
set(gcf,'paperposition',[0.25 1 7 5])
fontsize 12 10 10 10
if strcmpi(str,'print')
    print -deps morsie_morlet_wigdist.eps
end

figure(h.Number+2)
set(gcf,'paperposition',[0.25 1 7 5])
fontsize 12 10 10 10
if strcmpi(str,'print')
    print -deps morsie_morlet_wigdist_long.eps
end

figure(h.Number+3)
orient landscape
set(gcf,'paperposition',[1 1 10 6])
fontsize 12 10 8 10
if strcmpi(str,'print')
    print -deps morsie_wigdist_three.eps
end
%/*****************************************

%/*****************************************
p1=[ (1:.01:2) (2.1:.1:50) (55:5:100)];
alpha=(-1.5:.025:1.5)*2.3./1.5;

[p,alpha]=meshgrid(p1,alpha);
ga=p.*alpha+3;
be=squared(p)./ga;
index=find(ga<=0);
ga(index)=nan;
be(index)=nan;

%index=find(be>frac(ga-1,2));

[a,dt,dom]=morsebox(ga,be);
index=find(ga<1);
[fm,fe,fi,cf] = morsefreq(ga,be);

figure
subplot(2,3,1),h=contour(p./pi,alpha,ga,(3:10).^2,'k');hold on
h=contour(p./pi,alpha,ga,[2 3 4],'k','linewidth',2); 
h=contour(p./pi,alpha,ga,[2 3 4],'w','linewidth',1);xlog, 
[p1,a1]=morseprops(1,[(0:.1:100) (101:1:1000) (1000:10:10000)] );
plot(p1./pi,a1,'k','linewidth',2);
ylabel('Demodulate Skewness $\Im\{\alpha_{3;\beta,\gamma}\}$')
text(0.55,-2,'(a)  $\gamma$ parameter')

subplot(2,3,2),contour(p./pi,alpha,be,(1:10).^2,'k');xlog,hold on
h=contour(p./pi,alpha,be,[2 3 4],'k','linewidth',2);xlog,hold on
h=contour(p./pi,alpha,be,[2 3 4],'w','linewidth',1);xlog,hold on
contour(p./pi,alpha,be,[1 1],'k','linewidth',2);
title('Properties of Generalized Morse Wavelets')
text(0.55,-2,'(b)  $\beta$ Parameter')

subplot(2,3,3),contour(p./pi,alpha,a,(.52:.01:.59),'k');
xlog,hold on
contour(p./pi,alpha,a,[.51 0.51],'k','linewidth',2);
text(0.55,-2,'(c)  Heisenberg Area')

subplot(2,3,4),
x=(frac(fe-fm,fm));x(x>0.2)=0.2;
contour(p./pi,alpha,x,(0.025:.025:.2),'k'),hold on,xlog
x=-(frac(fe-fm,fm));x(x>0.2)=0.2;
contour(p./pi,alpha,x,(0.025:.025:.2),'k','linewidth',2),%
contour(p./pi,alpha,x,(0.025:.025:.2),'w:','linewidth',1),%
contour(p./pi,alpha,x,[0 0],'k','linewidth',2),%
ylabel('Demodulate Skewness $\Im\{\alpha_{3;\beta,\gamma}\}$')
xlabel('Duration $P_{\beta,\gamma}/\pi$')
%hold on,nocontours
text(0.55,-2,'(d)  Energy Freq. / Peak Freq. - 1 ')

subplot(2,3,5),
x=(frac(fi-fm,fm));x(x>0.2)=0.2;
contour(p./pi,alpha,x,(0.025:.025:.2),'k'),hold on,xlog
x=-(frac(fi-fm,fm));x(x>0.2)=0.2;
contour(p./pi,alpha,x,(0.025:.025:.2),'k','linewidth',2),%
contour(p./pi,alpha,x,(0.025:.025:.2),'w:','linewidth',1),%
contour(p./pi,alpha,x,[0 0],'k','linewidth',2),%
xlabel('Duration $P_{\beta,\gamma}/\pi$')
text(0.55,-2,'(e)  Inst. Freq. / Peak Freq. - 1')

subplot(2,3,6),
x=frac(1,2*pi)*(cf);x(x>0.2)=0.2;
contour(p./pi,alpha,x,(0.025:.025:.2),'k'),hold on,xlog
x=-frac(1,2*pi)*(cf);x(x>0.2)=0.2;
contour(p./pi,alpha,x,(0.025:.025:.2),'k','linewidth',2),%
contour(p./pi,alpha,x,(0.025:.025:.2),'w:','linewidth',1),%
contour(p./pi,alpha,x,[0 0],'k','linewidth',2),%
xlabel('Duration $P_{\beta,\gamma}/\pi$')
text(0.55,-2,'(f)  Frequency Curvature')

for i=1:6
    subplot(2,3,i)
    xlim([10^(-1/2) 10^1.5])
    xtick([1/2 1 2 4 8 16 32])
    %if ~verLessThan('matlab','8.4.0')
    %    %Fix weird labelling 
    %    h=get(gca,'xticklabel');
    %    h{2}='1';
    %    set(gca,'xticklabel',h)
    %end
    hlines(0,'k:')
    [p1,a1]=morseprops(1/100,[(0:.1:100) (101:1:1000) (1000:100:1000000)] );
    plot(p1./pi,a1,'k--','linewidth',1);
end

packfig(2,3)

orient landscape
set(gcf,'paperposition',[2  2 8 5])
fontsize 12 10 10 9
if strcmpi(str,'print')
    print -depsc morsie_parameters.eps
end
%\*****************************************


%END of jlab_makefigs_morsies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=jlab_makefigs_ridges(str)
%/********************************************************
%Load and prepare data

load npg2006
use npg2006

% This data was collected as part of the POMME experiment, see Le Cann,
% Assenbaum, Gascard, and Reverdin (2005), JGR.  

% num   -  Date in day number of 2001
% dt    -  Time step in days
% lat   -  Latitude
% lon   -  Longitude
% p     -  Pressure in decibar
% t     -  Potential temperature
% cx    -  Complex-valued position x+iy
% cv    -  Complex-valued displacement velocity u+iv

%Decide on frequencies
fs=morsespace(2,4,2*pi/10,2*pi/100,8);

cv=npg2006.cv;  %Annoyingly, I have to set cv by hand because of a Matlab bug with Simulink, 
                %which gets confused by 'use' on account of some other function called 'cv'

%Compute wavelet transforms using generalized Morse wavelets
P=sqrt(2*4);
[wx,wy]=wavetrans(real(cx),imag(cx),{2,4,fs,'bandpass'},'mirror');

%Convert to plus and minus transforms
[wp,wn]=vectmult(tmat,wx,wy);


%Uncomment to see plots of Cartesian and rotary transforms
%h=wavespecplot(num,cx,dt./fs,wx,wy,0.5);
%h=wavespecplot(num,cx,dt./fs,wp,wn,0.5);

%Form ridges of component time series
%Note... using hidden 'phase' option, read body of RIDGWEALK for details 
%In general, you will not need to use this option
[wrp,irp,jrp,frp]=ridgewalk(dt,wp,fs,P,7,'phase');   
[wrn,irn,jrn,frn]=ridgewalk(dt,wn,fs,P,7,'phase');   
[wrx,irx,jrx,frx]=ridgewalk(dt,wx,fs,P,7,'phase');   
[wry,iry,jry,fry]=ridgewalk(dt,wy,fs,P,7,'phase');    

%Map into time series locations
[wrx,frx]=ridgemap(length(cx),wrx,frx/2/pi,irx);
[wry,fry]=ridgemap(length(cx),wry,fry/2/pi,iry);
[wrp,frp]=ridgemap(length(cx),wrp,frp/2/pi,irp);
[wrn,frn]=ridgemap(length(cx),wrn,frn/2/pi,irn);

%Convert xy transforms to ellipse forms
[kappa,lambda,theta,phi]=ellparams(wrx,wry);

%Other ellipse properties from xy transforms
rm=ellrad(kappa,lambda,'geometric');
ra=ellrad(kappa,lambda,phi,'average');
ri=ellrad(kappa,lambda,phi,'instantaneous');
vm=ellvel(4*3600,kappa,lambda,theta,phi,1e5,'geometric');
vphi=ellvel(4*3600,kappa,lambda,theta,phi,1e5,'azimuthal');

%Other frequencies
[fphi,fth]=vdiff(2*pi*dt,phi,theta,1,'nans');
frp2=fphi+fth;frn2=fphi-fth;

%Elliptical signal
cxe=ellsig(kappa,lambda,theta,phi);
cxr=cx-cxe;

L=54;  %approximate region of edge-effects

%/********************************************************  
figure,
subplot(221),plot(num,[frx fry]),linestyle 2k k--
title('Diagnosed frequencies') 
axis([min(num) max(num) -.05 .37]),fixlabels([0 -2]),hlines(0,'k:'),vlines(num([L length(cx)-L]),'k:')
ylabel('Frequency (Cycles / day)'),
subplot(223),plot(num,frp),hold on,plot(num,frn),plot(num,frn2)
linestyle 2k k-- k-.
xlabel('Day of 2001'),ylabel('Frequency (Cycles / day)')
axis([min(num) max(num) -.05 .37]),fixlabels([0 -2]),hlines(0,'k:'),vlines(num([L length(cx)-L]),'k:')

eps=.01/2;
subplot(222),plot(num,frp+eps),hold on,plot(num,frp2-eps),linestyle 2k k--
title('Inferred frequencies') ,
axis([min(num) max(num) -.05 .37]),fixlabels([0 -2]),hlines(0,'k:'),noylabels,vlines(num([L length(cx)-L]),'k:')
%ylabel('\omega_+/2\pi (Cycles / day)'),
subplot(224),plot(num, [fphi fth]),linestyle 2k k--
%ylabel('\omega_\phi/2\pi and \omega_\theta/2\pi (Cycles / day)'),
axis([min(num) max(num) -.05 .37]),fixlabels([0 -2]),hlines(0,'k:'),noylabels,vlines(num([L length(cx)-L]),'k:')

letterlabels(2)
xlabel('Day of 2001'),packfig(2,2)
fontsize jpofigure
set(gcf,'paperposition',[1 1 7 4])

if strcmpi(str,'print')
   print -deps npg-2006-0054-f04.eps
end
%\********************************************************  

%/********************************************************  
figure
subplot(211)
plot(num,[ri rm rm ra abs(wn(:,27))/sqrt(2)]),linestyle  0.5k 3w 1.5k  0.5k-- k-.
ylabel('Radius (kilometers)'),
axis([min(num) max(num) 0 25]),vlines(num([L length(cx)-L]),'k:')
title('Radius and Temperature') ,
subplot(212)
plot(num,[t vfilt(t,12) vfilt(t,12)]),linestyle  0.5k 4w 1.5k  
xlabel('Day of 2001'),ylabel('Temperature ( $^\circ$ C)'),
axis([min(num) max(num) 12.1 12.86 ]),vlines(num([L length(cx)-L]),'k:')

letterlabels(1)
packfig(2,1,'rows')
fontsize jpofigure
set(gcf,'paperposition',[1 1 7 4])

if strcmpi(str,'print')
   print -deps npg-2006-0054-f05.eps
end
%\********************************************************  

%/\********************************************************  
figure

index=L:length(cx)-L;

r1=(1e-10:.1:25)';

subplot(121)
plot(ri(index),-vphi(index),'k'),hold on,
xlabel('Radius (kilometers)'),
ylabel('Azimuthal velocity (cm/s)'),
title('Instantaneous properties') ,
plot(r1,100./r1,'k:'),plot(r1,200./r1,'k:'),plot(r1,400./r1,'k:'),plot(r1,50./r1,'k:')
axis([0 22 0 25])
%plot(r1,-vq,'k--'),

subplot(122)
plot(ri(index),-vphi(index),'k'),hold on,linestyle 0.5D
plot(rm(index),-vm(index),'ko','markersize',2,'markerfacecolor','k'),
xlabel('Radius (kilometers)'),
title('Geometric mean properties') ,
plot(r1,100./r1,'k:'),plot(r1,200./r1,'k:'),plot(r1,400./r1,'k:'),plot(r1,50./r1,'k:')
axis([0 22 0 25])
%plot(r1,-vq,'k--'),


letterlabels(1)
packfig(1,2,'columns')

fontsize jpofigure
set(gcf,'paperposition',[1 1 7 3])

if strcmpi(str,'print')
   print -depsc npg-2006-0054-f06.eps
end
%\*****************************************************  


%/********************************************************
figure,
subplot(121)
h=plot(cx,'k');hold on
axis equal,axis([-90 80 -80 65]),
title('Eddy-trapped float')
xtick(-75:25:75),ytick(-75:25:75)
ylabel('Displacement North (km)')
xlabel('Displacement East (km)')
plot(cx(1),'k*','markersize',10)
subplot(122)
index=(6*4:6*4:length(kappa)-6*4);
ellipseplot(kappa,lambda,theta,cxr,'npoints',64,'index',index)
hold on,linestyle k,plot(cxr,'k:') 
axis equal,axis([-90 80 -80 65]),
xtick(-75:25:75),ytick(-75:25:75)
title('Ellipse extraction')
xlabel('Displacement East (km)')
letterlabels(1)
packfig(1,2,'columns')

fontsize jpofigure
set(gcf,'paperposition',[1 1 7 4])

if strcmpi(str,'print')
   print -deps npg-2006-0054-f01.eps
end
%\********************************************************  


%/********************************************************
figure
subplot(121),plot(num,real([cx cxe-100 cxr cxr]))
linestyle 0.5k 0.5k 3w 1.5k
hlines(-100,'k:')
title('Float displacement East')
ylabel('Kilometers')
xlabel('Day of 2001')
axis([min(num) max(num) -120 80]),vlines(num([L length(cx)-L]),'k:')
subplot(122),plot(num,imag([cx cxe-100*sqrt(-1) cxr cxr]))
linestyle 0.5k 0.5k 3w 1.5k
hlines(-100,'k:')
axis([min(num) max(num) -120 80]),vlines(num([L length(cx)-L]),'k:')
title('Float displacement North')
xlabel('Day of 2001')
letterlabels(1)
packfig(1,2,'columns')

fontsize jpofigure
set(gcf,'paperposition',[1 1 7 3])

if strcmpi(str,'print')
   print -deps npg-2006-0054-f02.eps
end
%\********************************************************  


%/********************************************************
%Ellipsesketch

a=3;
b=2;

phi=pi/4;
th=pi/6;

[k,l]=ab2kl(a,b);
figure
ellipseplot(k,l,th,'npoints',64,'phase',phi),hold on,linestyle k
plot(rot(th+pi/2)*[0 1]*b,'k--')
plot(rot(th)*[0 1]*a,'k--')
plot(rot(th)*(a*cos(phi)+sqrt(-1)*b*sin(phi)),'k*')

title('Sketch of ellipse')
axis equal
axis([-3 3 -3 3])
vlines(0,':'),hlines(0,':')
ytick(1),xtick(1)

xi=(0:.1:th);
plot(1.25*rot(xi),'k');

xi=(th:.01:pi/2.8);
plot(1.5*rot(xi),'k');

text(2,1,'a')
text(-1,1.1,'b')
text(1.2,1.2,'$\phi$')
text(1.4,0.35,'$\theta$')

%fixlabels(-1)
fontsize jpofigure
set(gcf,'paperposition',[2 2 3.5 3.5])

if strcmpi(str,'print')
   print -deps npg-2006-0054-f03.eps
end
%!gv ellipsesketch.eps &  
%\********************************************************


% Simple walk-through exampleSS
% load npg2006
% use npg2006
% 
% %Take the wavelet transform of eastward velocity from a Lagrangian float
% 
% %Decide on frequencies
% K=1;  %K is always 1; we rarely use the higher-order wavelets
% gamma=3;  %Gamma should always be 3, see Lilly and Olhede (2009)
% beta=3;   %Increase beta for more wiggles, decrease for less, but keep beta>1
% fhigh=2*pi/10;  %High frequency in radians per sample point
% flow=2*pi/1000;  %Low frequency in radians per sample point
% D=4;  %Density (or overlap) of wavelets in frequency; D=4 should be fine
% fs=morsespace(gamma,beta,fhigh,flow,D);
% 
% %Compute wavelet transforms using generalized Morse wavelets
% wx=wavetrans(real(cx),{K,gamma,beta,fs,'bandpass'},'mirror');
% h=wavespecplot(num,real(cx),2*pi./fs,wx);

%END of jlab_makefigs_ridges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



