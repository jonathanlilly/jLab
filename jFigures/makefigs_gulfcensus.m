function[varargout]=makefigs_gulfcensus(str)
%MAKEFIGS_GULFCENSUS  Makes all figures for Lilly and Perez-Brunius (2021b).
%
%   This function makes all figures for 
%
%       Lilly, J. M. and P. Perez-Brunius (2021). Extracting statistically
%           significant eddy signals from large Lagrangian datasets using
%           wavelet ridge analysis, with application to the Gulf of Mexico.
%           In press at Nonlinear Processes in Geohysics.
% 
%   To use this function, you'll need to download several datasets
%
%   gomed.nc                  https://zenodo.org/record/3978803
%   GulfDriftersOpen.nc       https://zenodo.org/record/3985916
%   GulfFlowOneTwelfth.nc     https://zenodo.org/record/3978793
%
%   You'll need to edit this file to reflect your search paths, and you'll 
%   need to have jLab installed.
%
%   The figures are not going to look the same as in the paper, because 
%   GulfDriftersOpen.nc is missing several proprietary datasets that are 
%   included in the figures in the paper. 
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2020--2021 J.M. Lilly --- type 'help jlab_license' for details

if nargin==0
  str='print';
end

if strcmp(str,'--f')
     makefigs_gulfcensus('noprint');return
end

try
    %str='noprint';
    %This is where your files are kept
    basedir='/Users/lilly/Desktop/Dropbox/NetCDF/';
    
    %This is where you want things to be printed to
    gulfdir='/Users/lilly/Desktop/Dropbox/Projects/gulfcensus/figures';
    %--------------------------------------------------------------------------
    %loading in the eddy and gridded data
    ncload([basedir 'GulfFlow_OneTwelfth.nc'])
    gulfflow_onetwelfth=GulfFlow_OneTwelfth;
    clear GulfFlow_OneTwelfth
    ncload([basedir 'gomed.nc']);
    ncload([basedir 'gomed_noise.nc']);
    load gulfnoisedrifters
    %--------------------------------------------------------------------------
    %loading in and preparing gulfdrifters
    %ncload([basedir 'GulfDriftersAll.nc']);
    %gulfdrifters=GulfDriftersAll;
    %if you're using GulfDriftersOpen
    ncload([basedir 'GulfDriftersOpen.nc']);
    gulfdrifters=GulfDriftersOpen;
    %if you're using GulfDriftersDWDE
    %ncload([basedir 'GulfDriftersDWDE.nc']);
    %gulfdrifters=GulfDriftersDWDE;
    %--------------------------------------------------------------------------
    load jtopo
catch
    disp('Looks like Matlab can''t file the required data files.  See MAKEFIGS_GULFCENSUS for details.')
end

alpha=1/80;  %stretch factor for quiver arrows
set(0,'DefaultFigureColormap',feval('lansey'));
%gulfdir='/Users/lilly/Desktop/Dropbox/web/jmlilly/talks/lilly-os20/figures';

tweakcolorbar=['pos=get(gca,''position'');'...
    'cpos=get(hc,''position'');'...
    'set(hc,''position'',[cpos(1)+0.025 cpos(2)+cpos(4)/2 cpos(3)-0.05 cpos(4)/2]);'...
    'set(gca,''position'',pos)'];

tweakcolorbar2=['pos=get(gca,''position'');'...
    'cpos=get(hc,''position'');'...
    'set(hc,''position'',[cpos(1)+0.025 cpos(2)-1.5*cpos(4) cpos(3)-0.05 cpos(4)/2]);'...
    'set(gca,''position'',pos)'];

tweakmap_fewercontours=['plot([-84.5+1i*21.5 -84.5+1i*22],''k--'',''linewidth'',0.5);hold on;'...
    'plot([-86.9+1i*21.5 -84.5+1i*21.5],''k--'',''linewidth'',0.5);' setstr(10) ...
    'jfig axis|[-98.35 -81.5 18 30.75] latratio|25 topoplot|gray '...
    'eval|topoplot([],[],-1/2,''3w'') eval|topoplot([],[],-1/2,''2D'') ticksout fontsize|[11 10 10 10] '];

tweakmap_fewercontours_thinner=['plot([-84.5+1i*21.5 -84.5+1i*22],''k--'',''linewidth'',0.5);hold on;'...
    'plot([-86.9+1i*21.5 -84.5+1i*21.5],''k--'',''linewidth'',0.5);' setstr(10) ...
    'jfig axis|[-98.35 -81.5 18 30.75] latratio|25 topoplot|gray '...
    'eval|topoplot([],[],-1/2,''2w'') eval|topoplot([],[],-1/2,''1.3D'') ticksout fontsize|[11 10 10 10] '];
%\*************************************************************************

%/*************************************************************************
%creation of mean flow contour
count=[]; %this is just to keep Matlab from complaining when called as a function
%Velocity vectors with polysmooth
use gulfflow_onetwelfth
[long,latg]=meshgrid(lon,lat);
spd=vmean(abs(u+1i*v),3)*100;%convert m/s to cm/s

%Speed contour
[ds,xs,ys,zs,ws]=spheresort(latg,long,spd,count,lat,lon,50); 
[zhat,beta,aux]=polysmooth(ds,xs,ys,[],zs,[],0,50);

load jtopo
jj=find(jtopo.lon>=lon(1)&jtopo.lon<=lon(end)+1/128);
ii=find(jtopo.lat>=lat(1)&jtopo.lat<=lat(end)+1/128);
maxmax(abs(lon-jtopo.lon(jj)'))  %ok, these are basically the same
maxmax(abs(lat-jtopo.lat(ii)))
depth=-jtopo.topo(ii,jj);
zhat(depth<0)=nan;

%set area outside of GoM to nans
[xg,yg]=meshgrid(lon,lat);
zhat(yg<21.5&xg>-86.9)=nan;


%finding 70 cm/s contour 
figure
jpcolor(lon,lat,zhat),hold on
c=contour(lon,lat,zhat,[70 70],'m');
close

clear xc yc
index=find(c(1,:)==70);
index(end+1)=size(c,2)+1;
for i=1:length(index)-1
    xc{i,1}=c(1,index(i)+1:index(i+1)-1)';
    yc{i,1}=c(2,index(i)+1:index(i+1)-1)';
end
[~,ii]=max(cellength(xc));
xc=xc{ii};
yc=yc{ii};
make gulfmeanflow xc yc
matsave gulfmeanflow
%\*************************************************************************

%/*************************************************************************
use gulfdrifters
use gulfmeanflow     
    
figure
use gulfdrifters
for j=1:length(lat)
    lat{j}(filled{j}==1)=inf;
end
cellplot(lon,lat),hold on,eval(tweakmap_fewercontours)
plot(cellfirst(lon),cellfirst(lat),'w.','markersize',10*1.75)
plot(cellfirst(lon),cellfirst(lat),'k.','markersize',8*1.75)
h=plot(xc,yc);linestyle -h h 4.5w
h=plot(xc,yc);linestyle -h h 3k    
fontsize 12 12 12 12 
if strcmp(str,'print')
%    jprint(gulfdir,'gulfcensus_spaghetti','-r300')
    jprint(gulfdir,'fig01','-r225')
end
%\*************************************************************************

%/*************************************************************************
%ellipseplots by rossby number with length colored
%this code is really ugly, sorry 

use gomed
bool=omega_ast_bar>-1/2&rho(:,5)<0.1;
gomed_noninertial_significant=structindex(gomed,bool);
use gomed_noninertial_significant

pindex=periodindex(3600,omega,1);
mom=R;
for i=1:length(time)
    mom{i}=omega_ast_bar(i)+0*R{i};
end
for i=1:length(time)
    [kappa{i},xi{i},theta{i},latres{i},lonres{i},lat{i},lon{i},omega{i},omega_ast{i},R{i},V{i},mom{i}]=...
        vindex(kappa{i},xi{i},theta{i},latres{i},lonres{i},lat{i},lon{i},omega{i},omega_ast{i},R{i},V{i},mom{i},pindex{i},1);
end
cell2col(kappa,xi,theta,latres,lonres,lat,lon,omega,omega_ast,R,V,omega_ast,mom);

[~,sorter]=sort(abs(omega_ast));
vindex(kappa,xi,theta,latres,lonres,lat,lon,omega,omega_ast,R,V,omega_ast,mom,sorter,1);
lambda=sign(xi).*sqrt(1-squared(xi));

use gulfmeanflow 

figure
subplot(1,2,1)
bool=mom<0;
h1=ellipseplot(kappa/111,lambda,theta,lonres+sqrt(-1)*latres,[1/cosd(25) 1],'index',bool);
linecolor(h1,abs(omega_ast(bool)),0.05,0.5,'lansey');hold on
set(h1,'linewidth',0.75)
eval(tweakmap_fewercontours_thinner),caxis([0.05 0.5])
text(-98,30.155,['(a) Anticyclones'],'color','w')
%--------------------------------------------------------------------------
subplot(1,2,2)
bool=mom>0;%&abs(lambda)<0.95;
h1=ellipseplot(kappa/111,lambda,theta,lonres+sqrt(-1)*latres,[1/cosd(25) 1],'index',bool);
set(h1,'linewidth',0.75)
linecolor(h1,abs(omega_ast(bool)),0.05,0.5,'lansey');hold on
eval(tweakmap_fewercontours_thinner),caxis([0.05 0.5])
text(-98,30.155,['(b) Cyclones'],'color','w')
%--------------------------------------------------------------------------
ha=packfig(1,2,'both');
for i=1:2
    axes(ha(i)),
    h=plot(xc,yc);linestyle -h h 3w
    h=plot(xc,yc);linestyle -h h 2k
    if i==2
        hc=colorbar('SouthOutside');eval(tweakcolorbar2)
        hc.Label.Interpreter='latex';
        hc.Label.String='Nondimensional Instantaneous Frequency Magnitude $\left|\omega_*(t)\right|\equiv\left|\omega(t)/f(t)\right|$';
    end
end
pos=get(hc,'position');
set(hc,'position',[pos(1)-.3 pos(2) pos(3)+.21 pos(4)])
orient tall
fontsize 11 11 11 11
set(gcf,'paperposition',[0.25 0.25 11 7.75])
if strcmp(str,'print')
%    jprint(gulfdir,'alleddieswithrossbycolor','-r300')
    jprint(gulfdir,'fig02','-r300')
end
%\*************************************************************************

%/*************************************************************************
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
vlines(0,'0.5:'),hlines(0,'0.5:')
ytick(1),xtick(1)

xi=(0:.1:th);
plot(1.25*rot(xi),'k');

xi=(th:.01:pi/1.71);
plot(1.75*rot(xi),'k');

text(1.5,2.3,'a')
text(-1,0.8,'b')
text(0.3,2,'$\varphi$')
text(1.3,0.5,'$\theta$')

xtick([-3:3]),ytick([-3:3])
%fixlabels(-1)
fontsize jpofigure
set(gcf,'paperposition',[2 2 3.5 3.5])

imag(log(a.*cos(phi)+1i*b.*sin(phi)))
sqrt(frac(a.^2+b.^2,2))
frac(2*a*b,a.^2+b.^2)

if strcmp(str,'print')
%    jprint(gulfdir,'ellipsesketch','jpeg -r500')
    jprint(gulfdir,'fig03','-r500')
end
%\*************************************************************************


%/*************************************************************************
%Synthetic example eddy 
dt=0.1;
T=1000;
t=[0:dt:T]'; 
phio=t/20.*[1+t/T];
thetao=-t/200+pi/2;
kappao=10+t./max(t).*50;
lambdao=-((4/3)*t/T-(1/3));
lambdao(1:floor(end/4))=-1e-10;
zo=ellsig(kappao,lambdao,thetao,phio);
xp=anatrans(real(zo-mean(zo)));
yp=anatrans(imag(zo-mean(zo)));
[kappahat,lambdahat,thetahat]=ellparams(xp,yp);
%--------------------------------------------------------------------------
figure,
subplot(1,2,1)
plot(zo,'linewidth',1),set(gca,'DataAspectRatio',[1 1 1]),axis([-85 85 -60 60]),hold on,plot(0,0,'k.')
xlabel('Displacement East (km)')
ylabel('Displacement North (km)')
subplot(1,2,2)
h=plot(zo);set(gca,'DataAspectRatio',[1 1 1]),hold on
linestyle -h h 0.5D
he=ellipseplot(kappahat,lambdahat,thetahat,'index',250:500:length(t)-250,'npoints',64);
set(he,'linewidth',1)
plot(zo(250:500:length(t)-250),'ko','markerfacecolor','w','markersize',3)
axis([-85 85 -60 60])
plot(0,0,'k.')
letterlabels(1)
xlabel('Displacement East (km)')
packfig(1,2,'columns')
fontsize 7 7 7 7 
if strcmp(str,'print')
%    jprint(gulfdir,'gulfcensus_synthetic_planview','-r300')
    jprint(gulfdir,'fig04','-r500')
end
%--------------------------------------------------------------------------
%Let t=days
%U=-1*3600*24/100/1000;  %5cm/s = 4.32 km/day  1 cm * (1 km/100/1000 cms)
rng(2);
U=-1/2*3600*24/100/1000;  %cms 2  kmd 
uo=maternoise(dt,length(zo),1/4,1,1/10,'fast');
[f,s]=maternspec(dt,length(zo),1/4,1,1/10);
zeps=-cumsum(uo*100*1000/3600/24)*dt;
z4=zeps+zo+U*t;
syntheticeddy=eddyridges(corfreq(25)/3600,t,z4,2,1/256,sqrt(6),0,0);
%--------------------------------------------------------------------------
figure,
jj=26;
subplot(2,1,1)
plot(z4,'linewidth',1),hold on,plot(zeps+U*t,'linewidth',1)
set(gca,'DataAspectRatio',[1 1 1])
use syntheticeddy
h=plot(zres{jj});linestyle -h h 0.5k
axis([-520 10 -59.999 70])
ylabel('Displacement North (km)')
subplot(2,1,2)
use syntheticeddy
h=plot(zres{jj});linestyle -h h 0.5k
lambda=sign(xi{jj}).*sqrt(1-squared(xi{jj}));
he=ellipseplot(kappa{jj},lambda,theta{jj},zres{jj},'index',250:500:length(kappa{jj})-250,'npoints',64);
hold on,set(he,'linewidth',1)
set(gca,'DataAspectRatio',[1 1 1])
axis([-520 10 -60 70])
xlabel('Displacement East (km)')
ylabel('Displacement North (km)')
letterlabels(1)
packfig(2,1,'rows')
fontsize 7 7 7 7 
pos=get(gcf,'paperposition');
set(gcf,'paperposition',[pos(1:3) 3.75])
if strcmp(str,'print')
%    jprint(gulfdir,'gulfcensus_synthetic_noisy','-r300')
    jprint(gulfdir,'fig05','-r500')
end
%--------------------------------------------------------------------------
figure,
subplot(1,2,1)
xp=anatrans(real(detrend(z4)));
yp=anatrans(imag(detrend(z4)));
h=plot(real(xp)+1i*real(yp));set(gca,'DataAspectRatio',[1 1 1]),hold on
linestyle -h h 0.5D
[kappahat4,lambdahat4,thetahat4]=ellparams(xp,yp);
he=ellipseplot(kappahat4,lambdahat4,thetahat4,'index',250:500:length(t)-250,'npoints',64);
set(he,'linewidth',1)
plot(zo(250:500:length(t)-250),'ko','markerfacecolor','w','markersize',3)
axis([-85 85 -60 60])
plot(0,0,'k.')
xlabel('Displacement East (km)')
ylabel('Displacement North (km)')
subplot(1,2,2)
use syntheticeddy
h=plot(zhat{jj});set(gca,'DataAspectRatio',[1 1 1]),hold on
lambda=sign(xi{jj}).*sqrt(1-squared(xi{jj}));
linestyle -h h 0.5D
he=ellipseplot(kappa{jj},lambda,theta{jj},'index',250:500:length(kappa{jj})-250,'npoints',64);
set(he,'linewidth',1)
plot(zo(250:500:length(t)-250),'ko','markerfacecolor','w','markersize',3)
axis([-85 85 -60 60])
plot(0,0,'k.')
letterlabels(1)
xlabel('Displacement East (km)')
packfig(1,2,'columns')
fontsize 7 7 7 7 
if strcmp(str,'print')
%    jprint(gulfdir,'gulfcensus_synthetic_planview_noisy','-r300')
    jprint(gulfdir,'fig06','-r500')
end
%--------------------------------------------------------------------------
%wavelet transform
figure
use syntheticeddy.params
wtot=sqrt(squared(wp)+squared(wn));
fcor=abs(corfreq(24))*24; 
cv=vdiff(z4,1)./(dt*3600*24)*100*1000;
cv(2500:2502)=nan;
cv=fillbad(cv);
h=wavespecplot(t,cv,2*pi./fs*dt,wtot);
cvhat=vdiff(syntheticeddy.zhat{26},1)./(dt*3600*24)*100*1000;
axes(h(1)), ylabel('Velocity (cm/s)'),hlines(0,'0.5k:')
h1=plot(t,real(cvhat));linestyle -h h1 D
h1=plot(t,imag(cvhat));linestyle -h h1 D
h1=plot(t,real(cv));linestyle -h h1 T
h1=plot(t,imag(cv));linestyle -h h1 U

axes(h(2))
h1=plot(t,2*pi./(vdiff(phio,1)./dt+sqrt(1-lambdao.^2).*vdiff(thetao,1)./dt),'w');
linestyle -h h1 D
use syntheticeddy
for i=1:length(ir)
    plot(t(ir{i}),2*pi./omega{i},'w')
end
plot(t(ir{26}),2*pi./omega{26},'k','linewidth',2),ylim([5 200])
ylabel('Period (Days)')
xlabel('Time (Days)')
axes(h(1)),text(20,-13.5,'(a)')
axes(h(2)),text(20,160,'(b)','color','w')
set(gcf,'paperposition',[1 1 4 4]),fontsize 7.5 7.5 7.5 7.5
if strcmp(str,'print')
%    jprint(gulfdir,'synthetictransform','-r300')
    jprint(gulfdir,'fig08','-r300')
end
%\*************************************************************************

%/*************************************************************************
%wavelet plot
figure
t=1:1001;t=t-mean(t);
omegao=1/5;
[psi,psif]=morsewave(length(t),3,2,omegao,'bandpass');%using the 3,2 wavelet
subplot(1,2,1)
psi=psi./maxmax(abs(psi));
[h1,h2]=uvplot(t*omegao,psi);xlim([-8 8 ]),hold on,ylim([-1 1]*1.05)
set(h1,'linewidth',1),set(h2,'linewidth',1)
h=plot(t*omegao,abs(psi)); linestyle -h h 0.5D
h=plot(t*omegao,-abs(psi)); linestyle -h h 0.5D
xlabel('Nondimensional Time ($\breve t = \omega_{\beta,\gamma}t$)')
hlines(0,'0.25k:')
vlines(0,'0.25k:')
vlines(sqrt(6),'0.5D--')
vlines(-sqrt(6),'0.5D--')
text(15*omegao,0.8,'$\breve t=\pm P_{\beta,\gamma}$','interpreter','latex')
subplot(1,2,2)
om=2*pi*[0:length(psif)-1]/1001;
plot(om/omegao,psif),xlim(2*pi*[0 80]/1001/omegao),ylim([0 2.05]),hold on
psigauss=2*exp(-1/2*6*(om./omegao-1).^2);
h=plot(om./omegao,psigauss);linestyle -h h 0.5D
h=plot(om./omegao,psif);linestyle -h h T
fm=morsefreq(3,2);
vlines(1,'0.52k:'),vlines(1+1/sqrt(6),'0.5D--'),vlines(1-1/sqrt(6),'0.5D--')
plot(1+1/sqrt(6),2*exp(-1/2),'ko','markerfacecolor','w','markersize',3);
plot(1-1/sqrt(6),2*exp(-1/2),'ko','markerfacecolor','w','markersize',3);
xlabel('Nondimensional Frequency ($\breve \omega=\omega/\omega_{\beta,\gamma}$)')
text(0.3/omegao,1.6,'$\breve \omega=1\pm 1/P_{\beta,\gamma}$','interpreter','latex')
letterlabels(1)

set(gcf,'paperposition',[1 1 8 3]),fontsize 8 8 8 8
if strcmp(str,'print')
%   jprint(gulfdir,'gulfcensus_wavelet','-r300')
    jprint(gulfdir,'fig07','-r500')
end
%\*************************************************************************

%/*************************************************************************
%Example of ridge analysis
jj=1315;

use gulfdrifters
onegulfeddy=eddyridges(time{jj},lat{jj},lon{jj},2,1/64,sqrt(6),1,0);
timeo=datenum(1+floor(yearfrac(time{jj}(1))),1,1);%   2013

use onegulfeddy
time=num;
[mult,ztot]=ridgemult(ir,kr,zhat);

latresfull=latres;
lonresfull=lonres;
for i=1:length(ir)
    [latresfull{i},lonresfull{i}]=...
        sig2latlon(real(ztot{i}),imag(ztot{i}),lat{i},lon{i});
end

figure,
use gulfdrifters
subplot(1,3,1),eval(tweakmap_fewercontours)
h=plot(lon{jj},lat{jj},'linewidth',1.5);linestyle -h h T1.5
axis([-96.6 -93.0000 18.5 21.25])
text(-93.25,21.1,'(a)')
xlabel('Longitude'),ylabel('Latitude')
latratio(19.875)

subplot(1,3,2)
use gulfdrifters
h0=plot(time{jj}-timeo,lat{jj},'linewidth',1.5);hold on
use onegulfeddy
time=num;
ii=1;h1=plot(time{ii}-timeo,lat{ii}-latres{ii}+20);linestyle(h1,'1.5K')
ii=2;h2=plot(time{ii}-timeo,lat{ii}-latres{ii}+20);linestyle(h2,'1.5X')
ii=3;h3=plot(time{ii}-timeo,lat{ii}-latres{ii}+20);linestyle(h3,'1.5K')
ii=4;h4=plot(time{ii}-timeo,lat{ii}-latres{ii}+20);linestyle(h4,'1.5U')
ii=5:12;h5=cellplot(celladd(time(ii),-timeo),celladd(celladd(lat(ii),cellmult(-1,latres(ii))),20),'1.5');
linestyle(h5,'1.5V')
ylim([18.5 21.25]),xlim([-45 118]),xtick([-40:20:150])
h6=cellplot(celladd(time,-timeo),celladd(latresfull,0),'0.5W');
legend([h0 h1 h4 h2 h5(1) h6(1)],'Original','Low','Middle','High','Inertial','Residual','Orientation','horiz','location','south')
text(112.9,21.1,'(b)')
xlabel('Day of Year 2013')
h=packfig(1,3,'columns');
delete(h(3))
axes(h(2))
pos=get(gca,'position');
set(gca,'position',[pos(1) pos(2) pos(3)*2 pos(4)])

set(gcf,'paperposition',[1 1 12 3]),fontsize 9 9 9 9 
if strcmp(str,'print')
%    jprint(gulfdir,'exampleeddyplot','-r300')
    jprint(gulfdir,'fig09','-r500')
end
%--------------------------------------------------------------------------
figure
use gulfdrifters
cv=cellpair(u,v);
use onegulfeddy.params
wtot=sqrt(squared(wp)+squared(wn));
fcor=abs(corfreq(lat{jj}))*24; 
bool=abs(wp)>abs(wn);
[~,~,it,jt,ft,et]=ridgewalk(1/24,wp,wn,fs,sqrt(6),1,0,[2*fcor,1/64*fcor]);
[~,~,ip,jp,fp,ep]=ridgewalk(1/24,wp,wn,fs,sqrt(6),1,0,[2*fcor,1/64*fcor],'mask',bool);
[~,~,in,jn,fn,en]=ridgewalk(1/24,wp,wn,fs,sqrt(6),1,0,[2*fcor,1/64*fcor],'mask',~bool);
%wp(~bool)=nan;wn(bool)=nan;

timet=it;timet(~isnan(it))=time{jj}(it(~isnan(it)))-timeo;
timep=ip;timep(~isnan(ip))=time{jj}(ip(~isnan(ip)))-timeo;
timen=in;timen(~isnan(in))=time{jj}(in(~isnan(in)))-timeo;
Tp=2*pi./fp;Tn=2*pi./fn;Tt=2*pi./ft;
col2cell(timet,Tt);
col2cell(timep,Tp);
col2cell(timen,Tn);
h=wavespecplot(time{jj}-timeo,cv{jj}*100,2*pi./fs/24,wtot,wp,wn);
axes(h(1)),ylim([-99 99])
for i=2:4
    use gulfdrifters,axes(h(i))
    contour(time{jj}-timeo,2*pi./fs/24,abs(wp./wn)',[1 1],'w','linewidth',1/2);
    plot(time{jj}-timeo,2*pi./corfreq(lat{jj})/24,...
        'color',[1 1 1]/2)
    ylabel('Period (Days)')
end
axes(h(2)),use onegulfeddy
hl=cellplot(timep,Tp);linestyle(hl([1 3]),'1.75K'),linestyle(hl(4),'1.75U'),linestyle(hl(2),'1.75X')
hl=cellplot(timen,Tn,'1.75V');
axes(h(3)),
hl=cellplot(timep,Tp);linestyle(hl([1 3]),'1.75K'),linestyle(hl(4),'1.75U'),linestyle(hl(2),'1.75X')
axes(h(4))
hl=cellplot(timen,Tn,'1.75V');
axes(h(1)),ytick([-100:25:100])
for i=1:4,axes(h(i)),xtick([-40:10:130]),end
xlabel('Day of Year 2013') 
axes(h(1)),text(-41,-80,'(a) Velocity')
axes(h(2)),text(-41,60,'(b) Total')
axes(h(3)),text(-41,60,'(c) Cyclonic')
axes(h(4)),text(-41,60,'(d) Anticyclonic')

set(gcf,'paperposition',[1 1 9 8]),fontsize 10 10 10 10
if strcmp(str,'print')
%    jprint(gulfdir,'exampletransform','-r300')
    jprint(gulfdir,'fig10','-r300')
end
%\*************************************************************************    

%/*************************************************************************
use gulfnoisedrifters{1}
use gulfmeanflow     

figure
N=length(gulfdrifters.time);
eval(tweakmap_fewercontours),cellplot(lon(1:N),lat(1:N)),
plot(cellfirst(lon(1:N)),cellfirst(lat(1:N)),'w.','markersize',10*1.75)
plot(cellfirst(lon(1:N)),cellfirst(lat(1:N)),'k.','markersize',8*1.75)
topoplot([],[],-1/2,'3w') 
topoplot([],[],-1/2,'2D')
topoplot(axis,[],0,'3w')
h=plot(xc,yc);linestyle -h h 4.5w
h=plot(xc,yc);linestyle -h h 3k    
fontsize 12 12 12 12 
if strcmp(str,'print')
    jprint(gulfdir,'gulfcensus_noisespaghetti','-r300')
    jprint(gulfdir,'fig11','-r180')
end
%\*************************************************************************

%/*************************************************************************
%Scatter plot of all radius vs. velocity and rossby number
figure
use gomed
bool=omega_ast_bar>-1/2&rho(:,5)<0.1;
gomed_noninertial_significant=structindex(gomed,bool);
for i=1:3
    if i==1
        use gomed
    elseif i==2
        use gomed_noise
    elseif  i==3
        use gomed_noninertial_significant
    end
    %--------------------------------------------------------------------------
    %bool=cellfirst(len)>3.5&cellmean(Ro)>-0.5; 
    if i==2 
        ii=min(find(diff(gomed_noise.segment_id)<0));
        vindex(omega_ast,omega,xi,kappa,R_bar,V_bar,L,xi_bar,latres,1:ii,1);
        length(omega) %2025
    end  
    [L,sorter]=sort(L);
    vindex(R_bar,V_bar,omega_ast_bar,xi_bar,sorter,1);
    %length(find(omega_ast_bar>-1/2))
     
    subplot(3,3,i)
    %scatter(kappabar,V_bar,5*L_bar,L_bar,'filled'),hold on
    scatter(R_bar,V_bar*100,5*L,L,'filled'),hold on
    caxis([2*sqrt(3)/pi 10]),colormap lansey
    axis([0 149.99 -90 90])%,  axis([0 149.99 -115 115]),
    boxon,hlines(0,'0.5k:')
    ylabel('Average Velocity $\overline{V}$ (cm/s)','interpreter','latex')
    %title('Mean Eddy Properties in Radius / Velocity Space')
    xtick([0:25:200]),ytick([-100:20:100])
     
    %V/(fR)=fmax, V= fmax * f * R
    fmax=2;fmin=1/64;
    fact=(corfreq(24)/3600)*100*1000; 
    h(1)=plot([0+sqrt(-1)*0;200+200*sqrt(-1)*fmax*fact],'k','linewidth',2);
    plot([0+sqrt(-1)*0;200-200*sqrt(-1)*fmax*fact],'k','linewidth',2)
    h(2)=plot([0+sqrt(-1)*0;200-200*sqrt(-1)*1*fact],'-','color',[1 1 1]/2,'linewidth',2);
    h(3)=plot([0+sqrt(-1)*0;200+200*sqrt(-1)*0.1666*fact],'k--','linewidth',1);
    plot([0+sqrt(-1)*0;200-200*sqrt(-1)*0.1666*fact],'k--','linewidth',1)
    h(4)=plot([0+sqrt(-1)*0;200+200*sqrt(-1)*fmin*fact],'k','linewidth',1);
    h(5)=plot([0+sqrt(-1)*0;200-200*sqrt(-1)*0.5*fact],'color',[1 1 1]*0.7)
    plot([0+sqrt(-1)*0;200-200*sqrt(-1)*fmin*fact],'k','linewidth',1)
    %legend(h,' V/(fR) = \pm 2',' V/(fR) = - 1',' V/(fR) = \pm 1/6',' V/(fR) = \pm 1/64')
    if i==1
        text(137,-82,'(a)')
    elseif i==2
        text(137,-82,'(b)')
    elseif i==3
        text(137,-82,'(c)'),vlines([10 50],'0.5D:'),
    end
    %--------------------------------------------------------------------------
    subplot(3,3,i+3)
    scatter(R_bar,omega_ast_bar,5*L,L,'filled')
    caxis([2*sqrt(3)/pi 10]),colormap lansey
    axis([0 149.99 -1.19999 1.1999]),boxon,hlines(0,'0.5k:')
    %title('Mean Eddy Properties in Radius / Ro Space')
    xtick([0:25:200]),ytick([-1.2:.2:1.2])
    h(2)=hlines(1/6,'k--');hlines(-1/6,'k--'),h(1)=hlines(2,'2k');
    h(3)=hlines(1/64,'k');hlines(-1/64,'k'),h(4)=hlines(-1,'2F');hlines(-0.5,'C')
    xlabel('Average Radius $\overline{R}$ (km)','interpreter','latex')
    %xlabel('Root-Mean-Square Radius \kappa_\star (km)','interpreter','tex')
      ylabel('Average Nondimensional Frequency $\overline{\omega_*}$','interpreter','latex')
    if i==1
        text(137,-1.1,'(d)')
    elseif i==2
        text(137,-1.1,'(e)')
    elseif i==3
        text(137,-1.1,'(f)'),vlines([10 50],'0.5D:'),
    end
    if i==2
        %legend(h([1 2 3 4]),'Ro = \pm 2', 'Ro = \pm 1/6', 'Ro = \pm 1/64','Ro = -1','interpreter','latex');
        legend(h([1 4 5 2 3]),'$\overline{\omega_*} = \pm 2$', '$\overline{\omega_*} = -1$','$\overline{\omega_*} = -1/2$','$\overline{\omega_*} = \pm 1/6$', '$\overline{\omega_*} = \pm 1/64$','interpreter','latex');
    end
    if i==2
        hc=colorbar('South');
        hc.Label.String='Eddy Ridge Duration (Periods)';
        hc.Ticks=[1:10];
    end
    subplot(3,3,i+6)
    scatter(R_bar,xi_bar,5*L,L,'filled')
    caxis([2*sqrt(3)/pi 10]),colormap lansey
    axis([0 149.99 -1 1]),boxon,hlines(0,'0.5k:')
    xtick([0:25:200]),ytick([-1.2:.2:1.2])
    xlabel('Average Radius $\overline{R}$ (km)','interpreter','latex')
    ylabel('Average Circularity $\overline{\xi}$','interpreter','latex')
    if i==1
        text(137,-0.92,'(g)')
    elseif i==2
        text(137,-0.92,'(h)')
    elseif i==3
        text(137,-0.92,'(i)'),vlines([10 50],'0.5D:'),
    end
end
ha=packfig(3,3,'both');
axes(ha(5)),eval(tweakcolorbar)
pos=get(hc,'position');
set(hc,'position',[pos(1) pos(2)-.009 pos(3:4)])

set(gcf,'paperposition',[1 1 11 11]),fontsize 11 11 11 11
if strcmp(str,'print')
%    jprint(gulfdir,'rossbyscatterplot','-r300')
    jprint(gulfdir,'fig12','-r400')
end
%\*************************************************************************

%/*************************************************************************
%scatter plot of rotary coefficient

figure,    
use gomed
[L,sorter]=sort(L);
vindex(R_bar,V_bar,omega_ast_bar,xi_bar,sorter,1);
    
subplot(1,3,1)
scatter(abs(xi_bar),omega_ast_bar,5*L,L,'filled'),hold on
axis([0 1 -1.19999 1.1999]),boxon,hlines(0,'0.5k:')
ytick([-1.2:.2:1.2])
hlines(1/6,'k--'),hlines(-1/6,'k--')
hlines(1/64,'k'),hlines(-1/64,'k'),hlines(-1,'2F'),hlines(-0.5,'C')
xlabel('Magnitude of Average Circularity $|\overline\xi|$','interpreter','latex')
ylabel('Average Nondimensional Frequency $\overline{\omega_*}$','interpreter','latex')
caxis([2*sqrt(3)/pi 10])
text(0.05,1.1,'(a)')
xtick([0:0.1:1])

use gomed_noise%_noninertial_significant
ii=min(find(diff(gomed_noise.segment_id)<0));
vindex(L,R_bar,V_bar,omega_ast_bar,xi_bar,1:ii,1);
[L,sorter]=sort(L);
vindex(R_bar,V_bar,omega_ast_bar,xi_bar,sorter,1);
    
subplot(1,3,2)
scatter(abs(xi_bar),omega_ast_bar,5*L,L,'filled'),hold on
axis([0.001 1 -1.19999 1.1999]),boxon,hlines(0,'0.5k:')
ytick([-1.2:.2:1.2])
hlines(1/6,'k--'),hlines(-1/6,'k--')
hlines(1/64,'k'),hlines(-1/64,'k'),hlines(-1,'2F'),hlines(-0.5,'C')
xlabel('Magnitude of Average Circularity $|\overline{\xi}|$','interpreter','latex')
caxis([2*sqrt(3)/pi 10])
text(0.05,1.1,'(b)')
xtick([0:0.1:1])

subplot(1,3,3)
use gomed
bool=omega_ast_bar>-1/2&rho(:,5)<0.1;
gomed_noninertial_significant=structindex(gomed,bool);
use gomed_noninertial_significant
[L,sorter]=sort(L);
vindex(R_bar,V_bar,omega_ast_bar,xi_bar,sorter,1);

scatter(abs(xi_bar),omega_ast_bar,5*L,L,'filled'),hold on
axis([0.001 1 -1.19999 1.1999]),boxon,hlines(0,'0.5k:')
ytick([-1.2:.2:1.2])
hlines(1/6,'k--'),hlines(-1/6,'k--')
hlines(1/64,'k'),hlines(-1/64,'k'),hlines(-1,'2F'),hlines(-0.5,'C')
xlabel('Magnitude of Average Circularity $|\overline\xi|$','interpreter','latex')
caxis([2*sqrt(3)/pi 10])
text(0.05,1.1,'(c)')
xtick([0:0.1:1])

hc=colorbar('South');
hc.Label.String='Eddy Ridge Duration (Periods)';
hc.Ticks=[1:10];

ha=packfig(1,3,'columns');
axes(ha(3)),eval(tweakcolorbar)
pos=get(hc,'position');
set(hc,'position',[pos(1) pos(2)+.06 pos(3:4)])

set(gcf,'paperposition',[1 1 11 5]),fontsize 10 10 10 10
if strcmp(str,'print')
%    jprint(gulfdir,'polarizationscatterplot','-r300')
    jprint(gulfdir,'fig13','-r500')
end
%\*************************************************************************

%/*************************************************************************
bool=gomed.omega_ast_bar>-1/2;
bool_noise=gomed_noise.omega_ast_bar>-1/2;
r=0.1;

% %checking out the dependence of detected events on significance parameter
% for i=1:8
%     %length(find(gomed.rho(bool,i)<r))    
%     floor(sum(gomed.L(gomed.rho(bool,i)<r)))
% end
% %suggests using Lxi^4 ...maximized the number of oscillations

% %computing some numbers
% length(find(gomed.omega_ast_bar>-1/2))    %2520
% length(find(gomed.omega_ast_bar>-1/2&gomed.rho(:,5)<0.1))    %1033
% length(find(gomed.omega_ast_bar>-1/2&gomed.rho(:,5)<0.1&gomed.L>4))    %279
% length(find(gomed_noise.omega_ast_bar>-1/2&gomed_noise.rho(:,5)<0.1&gomed_noise.L>4))  

dx=1/4;
xmid=[-2:0.01:2];
Lbins=[0:0.01:20];
xibins=[0:0.001:1];
N=10;

x=gomed.omega_ast_bar;
xnoise=gomed_noise.omega_ast_bar;
y=(gomed.L).*abs(gomed.xi_bar).^4;
ynoise=(gomed_noise.L).*abs(gomed_noise.xi_bar).^4;
[rho_L_xi4v2,ymid,S,Snoise,Srat]=eddylevels(dx,xmid,Lbins,x,y,xnoise,ynoise,N,'symmetric','sort');
%aresame(gomed.rho(:,5),rho_L_xi4v2)  %true

figure,
MS=4;
N=length(find(isfinite(cell2col(GulfDriftersAll.time))));
subplot(1,3,1),jpcolor(xmid,ymid,log10(S/N)),hold on
subplot(1,3,2),jpcolor(xmid,ymid,log10(Snoise/N)),hold on
subplot(1,3,3),jpcolor(xmid,ymid,log10(Srat)),hold on
subplot(1,3,1),use gomed,plot(omega_ast_bar,L.*squared(xi_bar),'k.','markersize',MS)
text(-1.9,5.7,'(a)')
contour(xmid,ymid,log10(S/N),[-6 -5],'w','linewidth',1.5)
ii=min(find(diff(gomed_noise.segment_id)<0));
subplot(1,3,2),use gomed_noise,plot(omega_ast_bar(1:ii),L(1:ii).*squared(xi_bar(1:ii)),'k.','markersize',MS)
text(-1.9,5.7,'(b)')
contour(xmid,ymid,log10(Snoise/N),[-6 -5],'w','linewidth',1.5)
subplot(1,3,3),use gomed
text(-1.9,5.7,'(c)')
bool=rho(:,5)<0.01;
plot(omega_ast_bar(bool),L(bool).*(xi_bar(bool)).^4,'k.','markersize',MS)
bool=rho(:,5)>0.01&rho(:,5)<0.1;
plot(omega_ast_bar(bool),L(bool).*(xi_bar(bool)).^4,'.','color',[1 1 1]/2,'markersize',MS)
bool=rho(:,5)>0.1;
plot(omega_ast_bar(bool),L(bool).*(xi_bar(bool)).^4,'w.','markersize',MS)
contour(xmid,ymid,log10(Srat),[-1 -1],'k','linewidth',3)
contour(xmid,ymid,log10(Srat),[-2 -2],'w','linewidth',3)
contour(xmid,ymid,log10(Srat),[-1 -1],'w','linewidth',2)
contour(xmid,ymid,log10(Srat),[-2 -2],'k','linewidth',2)

for i=1:3
    subplot(1,3,i)
    axis([-1.99 1.99 0 6])
    xlabel('Average Nondimensional Frequency $\overline{\omega_*}$','interpreter','latex')
    ylabel('Significance Parameter $L\overline{\xi}^4$','interpreter','latex')
    xtick([-2:1/2:2]),ytick([1:6])
    if i~=3,caxis([-7 -4]),else caxis([-3 0]),end
    set(gca,'DataAspectRatio',[1 2.25*(6/8) 1])
    vlines(-0.5,'1.5w')
    vlines(-0.5,'1C')
end
ha=packfig(1,3,'columns');
axes(ha(1)),hc=colorbar('SouthOutside');
hc.Label.Interpreter='latex';
hc.Label.String='Log10 Observed Survival Function $\mathcal{S}_{L\xi^4}$';
eval(tweakcolorbar)
axes(ha(2)),hc=colorbar('SouthOutside');
hc.Label.Interpreter='latex';
hc.Label.String='Log10 Noise Survival Function $\mathcal{S}^\varepsilon_{L\xi^4}$';
eval(tweakcolorbar)
axes(ha(3)),hc=colorbar('SouthOutside');
hc.Label.Interpreter='latex';
hc.Label.String='Log10 Confidence Level $\rho_{L\xi^4}\equiv \mathcal{S}_{L\xi^4}^\varepsilon/ \mathcal{S}_{L\xi^4}$';
eval(tweakcolorbar)
%--------------------------------------------------------------------------
set(gcf,'paperposition',[1 1 11 7.3]),fontsize 10 10 10 10
if strcmp(str,'print')
%    jprint(gulfdir,'significanceplot','-r300')
    jprint(gulfdir,'fig14','-r300')
end
%\*************************************************************************

% check for tides
% use gomed
% omega_ast_bar=cellmean(omega_ast,'weight',kappa);
% R_bar=sqrt(cellmean(cellmult(R,R)));
% bool=omega_ast_bar>-1/2&rho_bar<0.1&R_bar<5;
% gomed_noninertial_small=structindex(gomed,bool); 
% use gomed_noninertial_small
% figure,cellplot(time,omega),
% hlines(tidefreq('k1')/3600)

