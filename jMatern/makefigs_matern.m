function[varargout]=makefigs_matern(str)
%MAKEFIGS_MATERN  Make figures for Lilly et al. (2017).
% 
%   This function makes all figures for the paper 
%
%   Lilly, J. M., A. M. Sykulski, J. J. Early, and S. C. Olhede (2017).
%       Fractional Brownian motion, the Matern process, and stochastic 
%       modeling of turbulent dispersion.  Nonlinear Processes in
%       Geophysics, 24: 481--514.
%
%   To use it, you'll need to download the file materndata.zip from 
%   http://jmlilly.net/ftp/pub, unzip it, and put the files in it in a 
%   directory that is on your Matlab search path.  
%
%   Please note that some of these calculations are fairly intensive and
%   may take a little while. 
%
%   Usage: makefigs_matern
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2017--2023 J.M. Lilly --- type 'help jlab_license' for details

if nargin==0
    str='noprint';
end

% This portion creates a data file in Matlab format from the model output.
% You won't need to use it, as the Matlab format output is distributed 
% at http://jmlilly.net/ftp/pub.
% %/*************************************************************************
% %Create the qg6snapshot dataset
% filename='/Users/Shared/qg/Experiment_06/Experiment6TurbulenceIsotropic.nc';
% 
% x=FieldsFromTurbulenceFile(filename,[], 'x')'/1000;
% y=FieldsFromTurbulenceFile(filename,[], 'y')/1000;
% psi=FieldsFromTurbulenceFile(filename,1, 'psi');
% u=FieldsFromTurbulenceFile(filename,1, 'u')*100;
% v=FieldsFromTurbulenceFile(filename,1, 'v')*100;
% cv=u+sqrt(-1)*v;
% zeta=FieldsFromTurbulenceFile(filename,1, 'rv');
% N=FieldsFromTurbulenceFile(filename,1, 'strain_n');
% S=FieldsFromTurbulenceFile(filename,1, 'strain_s');
% 
% description=['Snapshot from forced-dissipative QG Experiment 6 by Jeffrey Early' setstr(10) setstr(10)...
%            'x:    Zonal position in km' setstr(10)...
%            'y:    Meridional position in km' setstr(10)...
%            'psi:  Streamfunction in m^2/s' setstr(10)...
%            'cv:   Complex velocity u+iv in cm/s' setstr(10)...
%            'zeta: Vorticity in 1/s'  setstr(10)...
%            'N:    Normal strain in 1/s' setstr(10)...
%            'S:    Shear strain in 1/s'];
% creator='J. M. Lilly -- makefigs_randomwalk script';
% timestamp=datestr(now);
% matsave qg6snapshot description creator timestamp x y psi cv zeta N S
% 
% %Verify that this is the same
% %xf=ncread(filename, 'x-position')'/1000;
% %yf=ncread(filename, 'y-position')'/1000;
% %\************************************************************************
% 
% %/*************************************************************************
% %Create the qg6trajectories dataset
% filename='/Volumes/Bifrost/early/qg/Experiment_06/Experiment6Trajectories50000.mat';
% cx=double((xpos+1i*ypos)/1000);  
% cv=vdiff(cx,1)*1000*100/(3600*6); 
% num=double(t)/86400;
% vindex(num,cx,cv,1:4:size(cv,1),1);
% dt=num(2)-num(1);   %Sample interval in days
% lat=24;             %Latitude
% L=2500;             %Domain size in km
% description=['Particle trajectories from forced-dissipative QG Experiment 6 by Jeffrey Early' setstr(10) setstr(10)...
%              'dt:   Sample interval in days' setstr(10)...
%              'lat:  Latitude of f-plane' setstr(10)...
%              'L:    Domain size in km' setstr(10)...
%              'num:  Time in days' setstr(10)...
%              'cx:   Complex position cx=x+iy in km'  setstr(10)...
%              'cv:   Complex velocity cv=u+iv in cm/s'];
% creator='J. M. Lilly -- makefigs_randomwalk script';
% timestamp=datestr(now);
% matsave qg6trajectories description creator timestamp dt lat L num cx cv
% %\*************************************************************************


% Because performing the fit to the Matern process is time-consuming, the
% fits are also distirbuted, in the file qgmodelfit.  However, you can run
% the fits yourself if you'd like to.
% %/*************************************************************************
% %Fit to background, 0 to 1.5 rad/day
% load qg6trajectories
% use qg6trajectories 
% vindex(num,cx,cv,1:1095,1);  %%Three years of daily data 1095=365*3
% 
% %Sort by mean spin
% spin=imag(vdiff(cv,1).*conj(cv));%./abs(cv);
% meanspin=mean(spin,1);
% [temp,sorter]=sort(log10(abs(meanspin)));
% vindex(meanspin,spin,cv,cx,sorter,2);
% 
% index=1:512;
% vindex(cx,cv,meanspin,index,2);
% psi=sleptap(size(cv,1),10,1);
% [f,spp,snn]=mspec(dt,cv,psi,'parallel');
% 
% make trajectories.model num cv cx f spp snn
% 
% %f(2)./mean(fcor(:)) = 0.044
% clear fit
% use trajectories.model
% tic;fit=maternfit(dt,cv,1,0,1.5,'both','taper',psi,'parallel');toc
% 
% %Matern fit
% use fit
% [tau,R]=materncov(dt,size(cv,1),sigma,alpha,lambda);
% [f,S]=blurspec(dt,R,'tapered',psi);  
% rng(0);
% cv=maternoise(dt,size(cv,1),sigma,alpha,lambda);
% cx=cumsum(cv,1).*dt;
% cx=cx+vrep(trajectories.model.cx(1,:)-cx(1,:),size(cx,1),1);
% [f,spp,snn]=mspec(dt,cv,psi,'parallel');
% make trajectories.matern num cv cx f spp snn S
% 
% matsave qgmodelfit fit trajectories
% %/************************************************************************
% 
% %/************************************************************************
% %Two supplemental experiments: 
% use qgmodelfit
% dt=1;
% %White noise
% use fit 
% rng(0);
% N=length(trajectories.model.num);
% kappa=maternprops(sigma,alpha,lambda);
% %Why do I need the 1/4 here?
% cv=vrep(sqrt(kappa/4)',N,1).*frac(1,sqrt(2)).*(randn(N,length(sigma))+1i*randn(N,length(sigma)));
% cx=cumsum(cv,1).*dt;
% %cx=cx+vrep(trajectories.model.cx(1,:)-cx(1,:),size(cx,1),1);
% psi=sleptap(size(cv,1),10,1);
% [f,spp,snn]=mspec(dt,cv,psi,'parallel');
% make trajectories.whitenoise num cv cx f spp snn 
% 
% %Make fBm that has the Matern slope
% dt=1;
% use fit 
% lambdao=2*pi/length(qgmodelfit.trajectories.model.num);
% sigmao=sigma.*frac(lambda,lambdao).^(alpha-1/2);
% [f,S]=maternspec(dt,size(cv,1),sigmao,alpha,lambdao); 
% cv=maternoise(dt,size(cv,1),sigmao,alpha,lambdao);
% cx=cumsum(cv,1).*dt;
% cx=cx+vrep(trajectories.model.cx(1,:)-cx(1,:),size(cx,1),1);
% psi=sleptap(size(cv,1),20,1);
% [f,spp,snn]=mspec(dt,cv,psi,'parallel');
% make trajectories.fBm num cv cx f spp snn S
% 
% matsave qgmodelfit fit trajectories
% %/************************************************************************
%End of data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%/**********************************************************************
%Heimdall: 1000 = 3x, 2000=7x, 4000=11x, 8000=45x

N=2000;  %In the text I use 8000, making shorter verison here for speed
%N=8000;

dt=1;
alpha=[1 1.5 2 3 4];
h=[.01 .02 .05 .2 1];
[alpha,h]=meshgrid(alpha,h);
h=h.*alpha;
h=h';alpha=alpha';

rng(1);  %set seed
tic;z=maternoise(dt,N,1,alpha,h);etime1=toc;
tic;z2=maternoise(dt,N,1,alpha,h,'fast');etime2=toc;

[psi,lambda]=sleptap(N,8);
[om,spp,snn]=mspec(dt,z,psi,lambda,'adaptive');    
[om,spp2,snn2]=mspec(dt,z2,psi,lambda,'adaptive');    

[f,sppo,snno]=maternspec(dt,N,1,alpha,h);

disp(['MATERNOISE fast algorithm was ' num2str(etime1./etime2) ' times faster than Cholesky algorithm.'])
for i=1:size(spp,2)
    spp(:,i)=spp(:,i).*sqrt(10).^(i-1);
    spp2(:,i)=spp2(:,i).*sqrt(10).^(i-1);
    sppo(:,i)=sppo(:,i).*sqrt(10).^(i-1);
end

figure
subplot(1,2,1),plot(om,spp,'k'),hold on,h=plot(om,sppo);
linestyle -h h D,xlog,ylog,title('Matern Process, Cholesky Algorithm')
axis tight,ylim([10^-5 10^12.5]),xlabel('Frequency (rad)'),text(10^-3,10^-4,'(a)')
ylabel('Power Spectral Density (offset)')
subplot(1,2,2),plot(om,spp2,'k'),hold on,h=plot(om,sppo);
linestyle -h h D,xlog,ylog,title('Matern Process, Fast Algorithm')
axis tight,ylim([10^-5 10^12.5]),xlabel('Frequency (rad)'),text(10^-3,10^-4,'(b)')
packfig(1,2,'columns')

cd /Users/lilly/Desktop/Dropbox/Projects/matern
set(gcf,'paperposition',[1 1 11 6])
fontsize 14 12 12 12
if strcmpi(str,'print')
    print -dpng matern_algorithms
    crop matern_algorithms.png
end
%/**********************************************************************


%/**********************************************************************
alpha=[1/2+1/100 5/8:1/8:2];
N=4096;
lambda=1/20;
[f,s]=maternspec(1,N,1,alpha,lambda);
[tau,R]=materncov(1,N,1,alpha,lambda);
tauG=tau/10;
G=maternimp(tauG,alpha,1/20);
 
figure
subplot(1,3,1),plot(f/lambda,s),xlog,ylog,axis tight,ylim([10^-3 10^2])
linestyle k E k E k E k E k E k E 3k
title('Matern Spectrum'),xlabel('Nondimensional Frequency $\omega/\lambda$')
ylabel('Power Spectral Density')
subplot(1,3,2),plot(tau*lambda,R),xlim([0 4]),ylim([0 1])
linestyle k E k E k E k E k E k E 3k
title('Matern Autocovariance Function'),xlabel('Nondimensional Time $\tau\lambda$')
ylabel('Autocovariance Amplitude')
subplot(1,3,3),plot(tauG*lambda,G),xlim([0 4]),ylim([0 0.3333])
linestyle k E k E k E k E k E k E 3k
title('Matern Green''s Function'),xlabel('Nondimensional Time $\tau\lambda$')
ylabel('Green''s Function Amplitude')
letterlabels(2)

set(gcf,'paperposition',[1 1 16 5])
fontsize 15 13 13 13
if strcmpi(str,'print')
    print -dpng matern_threepanel
    crop matern_threepanel.png
end
%\**********************************************************************

%/**********************************************************************
alpha=[0.501 5/8:1/8:3/2-1/8 1.499];

rng(1);
t=[1:4096];t=t-mean(t);t=1.1*t./max(t);
bool=abs(t)<1/4;
z=maternoise(1,length(t),1,alpha,0);
for i=1:size(z,2)
    z(:,i)=z(:,i)./vstd(z(:,i),1);
    z(:,i)=z(:,i)-vmean(z(bool,i),1);
end

alpha=fliplr(alpha);
z=fliplr(z);
for i=1:length(alpha)
    alphastr{i}=num2str(alpha(i),'%1.3f');
end


[x,y]=vzeros(5,length(alpha));
xo=(4);yo=(1/4).^(alpha-1/2);
x(1,:)=xo;y(1,:)=yo;
x(2,:)=-xo;y(2,:)=yo;
x(3,:)=-xo;y(3,:)=-yo;
x(4,:)=xo;y(4,:)=-yo;
x(5,:)=xo;y(5,:)=yo;

x=x/16;
xtilde=x*xo;
ytilde=y;
for i=1:size(x,2)
    ytilde(:,i)=y(:,i)./yo(i)+3*(i-1);
    y(:,i)=y(:,i)+3*(i-1);
end

ztilde=z;
for i=1:size(z,2)
    ztilde(:,i)=z(:,i)/yo(i);
end

figure
subplot(1,2,1),
plot(t(1:4:end),real(z(1:4:end,:))),yoffset 3, hold on, 
linestyle k E,ytick([-3:3:33])
h=plot(x,y);linestyle -h h 3w
h=plot(x,y);linestyle -h h G
title('Fractional Brownian Motion')
xlabel('Time'),ylabel('Real Part (offset)')
xlim([-1.1 1.1]),ylim([-2 27]),xtick([-1:.25:1]),text(-1.05,26.2,'(a)')
subplot(1,2,2),
plot(t*xo,real(ztilde)),yoffset 3, hold on,
linestyle k E,ytick([-33:3:3])
h=plot(xtilde,ytilde);linestyle -h h 3w
h=plot(xtilde,ytilde);linestyle -h h G
xlim([-1.1 1.1]),ylim([-2 27]),xtick([-1:.25:1]),text(-1.05,26.2,'(b)')
xlabel('Time'),
title('Fractional Brownian Motion, Similarity Rescaling')
%noxlabels,noylabels
for i=1:length(alpha),
    text(1.16,3*(i-1),['$\alpha$=' alphastr{i}],...
        'HorizontalAlignment','center','rotation',90);
end
h=packfig(1,2,'columns');
cd /Users/lilly/Desktop/Dropbox/Projects/matern
orient portrait
set(gcf,'paperposition',[1 1 11 8])
fontsize 14 12 12 10
if strcmpi(str,'print')
    print -dpng matern_brownian
    crop matern_brownian.png
end
%\**********************************************************************


%/**********************************************************************
alpha=[0.501 5/8:1/8:3/2-1/8 1.499];
figure
plot(z),xoffset 3,axis([-2 27 -3.1 3.1])
title('Fractional Brownian Motion, Plan View')
xlabel('Real Part (offset)'),ylabel('Imaginary Part')
set(gca,'DataAspectRatio',[1 1 1]),xtick([-3:3:27]),ytick([-3:1:3])
linestyle k E  
for i=1:length(alpha),
     text(3*(i-1),-2.7,['$\alpha=$' alphastr{i}],'HorizontalAlignment','center');
end
set(gca,'ticklen',[0.0100    0.0250]/2),set(gca,'tickdir','out')
set(gcf,'paperposition',[1 1 11 3])
fontsize 14 12 12 12
if strcmpi(str,'print')
    print -dpng matern_brownianplan
    crop matern_brownianplan.png
end
%\**********************************************************************

%/**********************************************************************
alpha=[1/2:1/8:2];
lambda=10.^[-3:-1];
rng(1);
t=[1:4096];t=t-mean(t);t=1.1*t./max(t);  
bool=abs(t)<1/4;
[lambdag,alphag]=meshgrid(lambda,alpha);
z=maternoise(1,length(t),1,alphag,lambdag);
for i=1:size(z,2)
    z(:,i)=z(:,i)./vstd(z(:,i),1);
    z(:,i)=z(:,i)-vmean(z(bool,i),1);
end
z=reshape(z,length(z),length(alpha),length(lambda));
for i=1:size(z,2)
    for j=1:size(z,3)
        z(:,i,j)=z(:,i,j)-1i*(j-1)*5-(i-1)*3;
    end
end
z=reshape(z,length(z),length(alpha)*length(lambda));


figure
[ax,p1,p2] = plotyy(real(z),imag(z),0,0);axis([-38 3 -13 3])
axes(ax(1))
title('Matern Process, Plan View')
set(gca,'DataAspectRatio',[1 1 1])
linestyle k E k E k E k E k E k E 3k
axis([-25.25 3 -2.5 3])
vlines(-25.25,'k:'),hlines(-2.5,'k:')
axis([-38 3 -13 3])
xlabel('(<== Increasing $\alpha$)         Real Part         (Decreasing $\alpha$ ==>)')
ylabel('Imaginary Part  ')
axes(ax(2))
set(gca,'ytick',[-3 -2 -1]),flipy,
ylim([-3.6 -0.4]),xlim([-3.6 -0.4]*41/16),set(gca,'DataAspectRatio',[1 1 1])
set(gca,'yticklabel',{'10^{-3}','10^{-2}','10^{-1}'})
ylabel('<== Increasing $\lambda$')
set(gca,'ycolor','k')

axes(ax(1))
cd /Users/lilly/Desktop/Dropbox/Projects/matern

set(gcf,'paperposition',[1 1 11 5])
fontsize 14 12 12 12
if strcmpi(str,'print')
    print -dpng matern_maternplan
    crop matern_maternplan.png
end
%\**********************************************************************

%/**********************************************************************
alpha=[1/2+1e-8:.01:3/2-1e-8]';
tau=[1:0.01:10]';
C=zeros(length(tau),length(alpha));
for i=1:size(C,2)
    C(:,i)=frac(1,2)*(abs(1+tau).^(2*alpha(i)-1)+abs(1-tau).^(2*alpha(i)-1)-2*abs(tau).^(2*alpha(i)-1));
end
figure,
jpcolor(tau,alpha,log10(abs(C'))),colormap gray,hold on, caxis([-2.5 0])
contour(tau,alpha,C',[.2:.2:2]/2,'linecolor','k')
[c,h]=contour(tau,alpha,C',[0 0],'linestyle','-','linewidth',2,'linecolor','w');
[c,h]=contour(tau,alpha,C',[-1:.02:0]/2,'linestyle','-','linecolor','w','linewidth',1);
hc=colorbar('EastOutside');
hc.Label.String='Log10 of Absolute Value of R_{zz}^{fGn}(\tau,1)/(A^2V_\alpha)';
title('Fractional Gaussian Noise Autocovariance')
ylabel('Slope Parameter $\alpha$')
xlabel('Nondimensional Time Offset $\tau/\Delta$')
text(5,0.55,'Anti-Persistence, $R_{zz}^{fGn}(\tau,1)< 0$','color','w')
text(5,1.05,'Persistence, $R_{zz}^{fGn}(\tau,1) > 0$','color','w')
ylim([1/2 3/2]),xtick([0:1:10])

set(gcf,'paperposition',[1 1 6 6])
fontsize 12 12 12 12
if strcmpi(str,'print')
    print -dpng matern_gaussiancovariance
    crop matern_gaussiancovariance.png
end
%\*************************************************************************

load qg6snapshot
load qgmodelfit
%load qg6trajectories
%/*************************************************************************
figure
subplot(1,2,1)
use qg6snapshot
jpcolor(x/1000,y/1000,abs(cv))
axis equal,axis([-1 1 -1 1]*1250*0.99/1000),hold on,colormap gray,caxis([0 18]),flipmap
title('Current Speed')

use qgmodelfit.trajectories.model
vindex(num,cx,cv,1:1095,1);  %Three years of daily data 1095=365*3
%Sort by mean spin
spin=imag(vdiff(cv,1).*conj(cv));%./abs(cv);
meanspin=mean(spin,1);
cx=trajwrap(cx,1250);
cx(abs(vdiff(cx,1))>100)=nan;
subplot(1,2,2),plot(cx/1000),axis equal,axis([-1 1 -1 1]*1250*0.99/1000),hold on
linestyle A B C D E F G H I J K 
plot(cx(1,:),'k.','markersize',6)
%xlabel('East-West Position (1000 km)')
title('Particle Trajectories')

h=packfig(1,2,'columns');
axes(h(2)),set(gca,'yaxislocation','right'),set(gca,'yticklabelmode','auto')

cd /Users/lilly/Desktop/Dropbox/Projects/matern
set(gcf,'paperposition',[1 1 11 6])
fontsize 14 12 12 12
if strcmpi(str,'print')
    print -dpng matern_snapshot
    crop matern_snapshot.png
end
%\************************************************************************


load qgmodelfit 
use qgmodelfit 
%/************************************************************************
figure
use trajectories.model
d=sqrt(vmean(abs(cx-vrep(cx(1,:),size(cx,1),1)).^2,2))/1000;
subplot(2,2,1),plot((cx-vrep(cx(1,:),size(cx,1),1))/1000),axis equal,
linestyle A B C D E F G H I J K 
axis([-1 1 -1 1]*1250/1000),hlines([-1250 1250]/1000,'2G'),vlines([-1250 1250]/1000,'2G')
axis([-1 1 -1 1]*2000*0.99/1000),plot(0,0,'k','markersize',8)
for i=1:6
    h=plot(d(floor(365/2)*i)*rot(0:0.01:2*pi+0.01));linestyle -h h 3w
    h=plot(d(floor(365/2)*i)*rot(0:0.01:2*pi+0.01));linestyle -h h 2k
end
%ylabel('North-South Displacement (1000 km)')
text(-1.73,2.2,'Dispersion for: (a) Geostrophic Turbulence, (b) Matern, (c) White Noise, (d) Power-Law Process','fontweight','bold')

use trajectories.matern
d=sqrt(vmean(abs(cx-vrep(cx(1,:),size(cx,1),1)).^2,2))/1000;
subplot(2,2,2),plot((cx-vrep(cx(1,:),size(cx,1),1))/1000),axis equal,
linestyle A B C D E F G H I J K 
axis([-1 1 -1 1]*1250/1000),hlines([-1250 1250]/1000,'2G'),vlines([-1250 1250]/1000,'2G')
axis([-1 1 -1 1]*2000*0.99/1000),plot(0,0,'k','markersize',8)
for i=1:6
    h=plot(d(floor(365/2)*i)*rot(0:0.01:2*pi+0.01));linestyle -h h 3w
    h=plot(d(floor(365/2)*i)*rot(0:0.01:2*pi+0.01));linestyle -h h 2k
end
%title('Integrated Matern Random Process')

use trajectories.whitenoise
d=sqrt(vmean(abs(cx-vrep(cx(1,:),size(cx,1),1)).^2,2))/1000;
subplot(2,2,3),plot((cx-vrep(cx(1,:),size(cx,1),1))/1000),axis equal,
linestyle A B C D E F G H I J K 
axis([-1 1 -1 1]*1250/1000),hlines([-1250 1250]/1000,'2G'),vlines([-1250 1250]/1000,'2G')
axis([-1 1 -1 1]*2000*0.99/1000),plot(0,0,'k','markersize',8)
for i=1:6
    h=plot(d(floor(365/2)*i)*rot(0:0.01:2*pi+0.01));linestyle -h h 3w
    h=plot(d(floor(365/2)*i)*rot(0:0.01:2*pi+0.01));linestyle -h h 2k
end
%xlabel('East-West Displacement (1000 km)')
%ylabel('North-South Displacement (1000 km)')
%title('Integrated White Noise Trajectories')

use trajectories.fBm
d=sqrt(vmean(abs(cx-vrep(cx(1,:),size(cx,1),1)).^2,2))/1000/1e6;
subplot(2,2,4),plot((cx-vrep(cx(1,:),size(cx,1),1))/1000/1e6),axis equal,hold on
linestyle A B C D E F G H I J K 
%axis([-1 1 -1 1]*1250/1000),hlines([-1250 1250]/1e9,'2G'),vlines([-1250 1250]/1e9,'2G')
axis([-1 1 -1 1]*2000*0.99/1000),plot(0,0,'k','markersize',8)
for i=1:6
    h=plot(d(floor(365/2)*i)*rot(0:0.01:2*pi+0.01));linestyle -h h 3w
    h=plot(d(floor(365/2)*i)*rot(0:0.01:2*pi+0.01));linestyle -h h 2k
end
set(gca,'YAxisLocation','right')
%xlabel('East-West Displacement (10^9 km)')
%title('Integrated Fractional Brownian Motion')
letterlabels(3)

h=packfig(2,2);
axes(h(2)),set(gca,'yaxislocation','right'),set(gca,'yticklabelmode','auto')
axes(h(4)),set(gca,'yaxislocation','right'),set(gca,'yticklabelmode','auto')
 
for i=1:4
    %axes(h(i)),xtick([-2:.5:2]),ytick([-2:.5:2])
    fixlabels(-1)
end

cd /Users/lilly/Desktop/Dropbox/Projects/matern
set(gcf,'paperposition',[1 1 8 7.6])
fontsize 10 10 10 10
if strcmpi(str,'print')
    print -dpng matern_dispersion
    crop matern_dispersion.png
end
%\************************************************************************


load qgmodelfit 
use qgmodelfit 
%/************************************************************************
figure
use trajectories.model
plot(f,vmean(spp,2));hold on 
use trajectories.matern
plot(f,vmean(spp,2));
use trajectories.whitenoise
plot(f,vmean(spp,2));
use trajectories.fBm
%plot(f,vmedian(s,2));
plot(f,vmean(spp,2));
xlog,ylog,ylim([10^-5.5 10^4.5]),xlim([min(f) max(f)])
plot(f,f.^(-8)/100)
%use trajectories.model
%[f,sppper,snnper]=mspec(1,cv,[]);
%plot(f,vmean(sppper,2));
ylabel('Power Spectral Density')
linestyle 2k k-- 2E 2C-- k
use trajectories.fBm
plot(f,0*f+maxmax(vmean(spp,2))/10^15,'k:')
ylabel('Power Spectral Density'),xlabel('Frequency (rad / day)')
%legend('Turbulence','Matern','White Noise','Power Law','Periodogram','location','best')
legend('Turbulence','Matern','White Noise','Power Law','Slope of -8','Noise Level','location','best')
title('Observed and Modeled Spectra')


cd /Users/lilly/Desktop/Dropbox/Projects/matern
set(gcf,'paperposition',[1 1 5 4.5])
fontsize 12 10 10 10
if strcmpi(str,'print')
    print -dpng matern_spectra
    crop matern_spectra.png
end
%\*************************************************************************

%% Compute some statistiscs of the fit if you're interested
% load qgmodelfit 
% disp('Alpha')
% mean(qgmodelfit.fit.alpha)
% %median(qgmodelfit.fit.alpha)
% std(qgmodelfit.fit.alpha)
% disp('Lambda')
% mean(qgmodelfit.fit.lambda)
% %median(qgmodelfit.fit.lambda)
% std(qgmodelfit.fit.lambda)
% disp('Sigma')
% mean(qgmodelfit.fit.sigma)
% %median(qgmodelfit.fit.sigma)
% std(qgmodelfit.fit.sigma)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Movies
if 0
%--------------------------------------------------------------------------
%Turbulence model movie
load qgmodelfit 
use qgmodelfit 
L=2500;             %Domain size in km
cx=qgmodelfit.trajectories.model.cx/1000;
cx=trajwrap(cx,L/1000,'nan');


figure
N=100;
cd /Users/lilly/Figures/maternmodelmovie

%Create the qg6snapshot dataset
%filename='/Users/Shared/qg/Experiment_06/Experiment6TurbulenceIsotropic.nc';
filename='/Users/lilly/Desktop/Experiment6TurbulenceIsotropic.nc';

x=FieldsFromTurbulenceFile(filename,[], 'x')'/1000/1000;
y=FieldsFromTurbulenceFile(filename,[], 'y')/1000/1000;

for i=1:5:1095
    clf
    u=FieldsFromTurbulenceFile(filename,i, 'u');
    v=FieldsFromTurbulenceFile(filename,i, 'v');
    %zeta=FieldsFromTurbulenceFile(filename,1, 'rv');
    %Ns=FieldsFromTurbulenceFile(filename,1, 'strain_n');
    %Ss=FieldsFromTurbulenceFile(filename,1, 'strain_s');
    %ow=(Ns.^2+Ss.^2)-zeta.^2;
    spd=sqrt(u.^2+v.^2);
    
    jpcolor(x,y,spd),colormap gray,caxis([-.2 .05]*1e-9),hold on,flipmap
    
    set(gca,'color',[248,248,255]/256)
    plot(cx(max(i-N+1,1):i,:)),hold on,plot(cx(max(i,2),:)+1i./10^9 ,'k.'),caxis([0 .2])
    %text(-2.4,-2.4,'Particle Dispersion in 2D Turbulence')
    
    text(0.6,-1.19,['Day ' int2str(i) '/' int2str(size(cx,1))]),

    set(gca,'xtick',[],'ytick',[]),
    axis equal,axis([-L L -L L]/2/1000)
    set(gca,'position',[0 0 1 1]);
    set(gcf,'papersize',[6 6]);
    set(gcf,'paperposition',[0 0 6 6]);
    eval(['print -dpng frame' int2str(i) '.png'])
    %eval(['crop frame' int2str(i) '.png'])
end

%--------------------------------------------------------------------------
%Dispersion movie
load qgmodelfit 
use qgmodelfit 
cx1=qgmodelfit.trajectories.model.cx/1000;
cx2=qgmodelfit.trajectories.matern.cx/1000;
cx1=cx1-vrep(cx1(1,:),size(cx1,1),1);
cx2=cx2-vrep(cx2(1,:),size(cx2,1),1);
std1=vstd(cx1,2);
std2=vstd(cx2,2);
L=2500;             %Domain size in km


%Web version
N=100;
cd /Users/lilly/Figures/materndispersionmovie

for i=1:5:100
    subplot(1,2,1)
    if i>1
        h=plot(rot([0:.01:2*pi+0.01])*std1(i));
        linestyle -h h 2G
    end
    hold on
    set(gca,'color',[248,248,255]/256)
    if i>1
        plot(cx1(max(i-N+1,1):i,:)),hold on,plot(cx1(i,:),'k.')
    else
        plot(0,0,'k.')
    end
    %text(-2.4,-2.4,'Particle Dispersion in 2D Turbulence')
    
    subplot(1,2,2)
    if i>1
        h=plot(rot([0:.01:2*pi+0.01])*std2(i));
        linestyle -h h 2G
    end
    hold on
    set(gca,'color',[248,248,255]/256)
    if i>1
        plot(cx2(max(i-N+1,1):i,:)),hold on,plot(cx2(i,:),'k.')
    else
        plot(0,0,'k.')
    end
    %text(-2.4,-2.4,'Particle Dispersion in a 3-Parameter Stochastic Model')
    text(1,-2.35,['Day ' int2str(i) '/' int2str(size(cx1,1))])

    h=packfig(1,2);
    for j=1:2
        axes(h(j)),set(gca,'xtick',[],'ytick',[]),box on
        %set(gca,'xcolor','none','ycolor','none')
        axis equal,axis([-L L -L L]/2/1000*2),
    end
    set(gcf,'paperposition',[1 1 12 7]);
    eval(['print -dpng frame' int2str(i) '.png'])
    eval(['crop frame' int2str(i) '.png'])
end
%--------------------------------------------------------------------------

N=100;
cd /Users/lilly/Figures/materndispersionmovie

for i=N+1:5:1095
    subplot(1,2,1)   
    h=plot(rot([0:.01:2*pi+0.01])*std1(i));hold on
    linestyle -h h 2G
    plot(cx1(i-N:i,:)),hold on,plot(cx1(i,:),'k.')
    title(['Dispersion of Simulation on Day ' int2str(i) '/' int2str(size(cx1,1))])
    xlabel('East / West Displacement (1000 km)')
    ylabel('North / South  Displacement (1000 km)')
    
    subplot(1,2,2)
    h=plot(rot([0:.01:2*pi+0.01])*std2(i));hold on
    linestyle -h h 2G
    plot(cx2(i-N:i,:)),hold on,plot(cx2(i,:),'k.')
    title(['Dispersion of Matern Process on Day '  int2str(i) '/' int2str(size(cx1,1))])
    xlabel('East / West Displacement (1000 km)')
    
    h=packfig(1,2);
    for j=1:2
        axes(h(j)),xtick([-2:1:2]),ytick([-2:1:2])
        axis equal,axis([-L L -L L]/2/1000*2)
    end
    set(gcf,'paperposition',[1 1 12 7]);
    eval(['print -dpng frame' int2str(i) '.png'])
end



load qgmodelfit 
use qgmodelfit 
cx1=trajwrap(qgmodelfit.trajectories.model.cx,L/2,nan)/1000;
cx2=trajwrap(qgmodelfit.trajectories.matern.cx,L/2,nan)/1000;

use qg6trajectories
vindex(num,cx,cv,1:1095,1);  %%Three years of daily data 1095=365*3
%Sort by mean spin
spin=imag(vdiff(cv,1).*conj(cv));%./abs(cv);
meanspin=mean(spin,1);%./mean(squared(cv),1);
[temp,sorter]=sort(log10(abs(meanspin)));
vindex(meanspin,spin,cv,cx,sorter,2);
cx=trajwrap(cx,L/2,nan)/1000;



N=100;
cd /Users/lilly/Figures/maternparticlemovie

for i=N+1:5:1095
    subplot(1,3,1)   
    plot(cx(i-N:i,:)),hold on,plot(cx(i,:),'k.')
    title(['All Simulation Particles on Day ' int2str(i) '/' int2str(size(cx1,1))])
    xlabel('East / West Displacement (1000 km)')
    ylabel('North / South  Displacement (1000 km)')
    
    subplot(1,3,2)   
    plot(cx(i-N:i,1:512)),hold on,plot(cx(i,1:512),'k.')
    title(['Eddy-Free Particles on Day ' int2str(i) '/' int2str(size(cx1,1))])
    xlabel('East / West Displacement (1000 km)')
    ylabel('North / South  Displacement (1000 km)')
    
    subplot(1,3,3)
    plot(cx2(i-N:i,:)),hold on,plot(cx2(i,:),'k.')
    title(['Matern Particles on Day '  int2str(i) '/' int2str(size(cx1,1))])
    xlabel('East / West Displacement (1000 km)')
    
    h=packfig(1,3);
    for j=1:3
        axes(h(j)),xtick([-2:1:2]/2),ytick([-2:1:2]/2)
        axis equal,axis([-L L -L L]/2/1000)
    end
    %figuresize(15,7,'inches')
    set(gcf,'paperposition',[1 1 15 7]);
    eval(['print -dpng frame' int2str(i) '.png'])
    %eval(['crop frame' int2str(i) '.png'])
end

end

