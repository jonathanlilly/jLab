function[varargout]=makefigs_bandwidth(str)
%MAKEFIGS_BANDWIDTH  Make figures for Lilly and Olhede (2010a).
% 
%   This function makes all figures for the paper 
%
%   Lilly, J. M., and S. C. Olhede (2010). Bivariate instantaneous 
%      frequency and bandwidth. IEEE Transactions on Signal Processing,
%      58 (2), 591--603.
%
%   Usage: makefigs_bandwidth
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2010--2023 J.M. Lilly --- type 'help jlab_license' for details

if nargin==0
    str='noprint';
end

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

