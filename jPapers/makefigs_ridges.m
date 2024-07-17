function[varargout]=makefigs_ridges(str)
%MAKEFIGS_RIDGES  Make figures for Lilly, Scott, and Olhede (2011).
% 
%   This function makes all figures for the paper 
%
%   Lilly, J. M., R. K. Scott, and S. C. Olhede (2011). Extracting waves 
%      and vortices from Lagrangian trajectories. Geophysical Research 
%      Letters, 38, L23605, 1--5.
%
%   Usage: makefigs_ridges
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2023 J.M. Lilly --- type 'help jlab_license' for details

if nargin==0
    str='noprint';
end

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



