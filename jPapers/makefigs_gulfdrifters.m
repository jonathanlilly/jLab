function[varargout]=makefigs_gulfdrifters(str)
%MAKEFIGS_GULFDRIFTERS  Make figures for Lilly and Perez-Brunius (2021a).
%
%   This function makes all figures for the paper
%
%   Lilly, J. M. and P. Perez-Brunius (2021a). A gridded surface current
%       product for the Gulf of Mexico from consolidated  drifter
%       measurements.  Earth System Science Data, 13: 645–669.
% 
%   To use this function, you'll need to download several datasets
%
%   GulfDriftersOpen.nc       https://zenodo.org/record/3985917
%   GulfFlowOneQuarter.nc     https://zenodo.org/record/3978794
%   GulfFlowOneTwelfth.nc     https://zenodo.org/record/3978794
%   drifter_monthlymeans.mat  https://www.aoml.noaa.gov/phod/gdp/mean_velocity.php
%
%   You'll need to edit this file to reflect your search paths, and you'll 
%   need to have jLab installed.
%
%   This file is provided for completeness.  Some of the figures are not 
%   going to look the same as in the paper, because GulfDriftersOpen.nc is
%   missing several proprietary datasets that are included in the figures 
%   in the paper.  Also, one of the figures uses CMEMS altimetry that is 
%   not redistributed here.
%
%   Usage: makefigs_gulfdrifters
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2020--2023 J.M. Lilly --- type 'help jlab_license' for details

if nargin==0
    str='noprint';
end

try
    %str='noprint';
    %This is where your files are kept
    basedir='/Users/lilly/Desktop/Dropbox/NetCDF/';
    
    %This is where you want things to be printed to
    gulfdir='/Users/lilly/Desktop/Dropbox/Projects/gulfdrifters/figures';
    %--------------------------------------------------------------------------
    %loading in the gridded data
    ncload([basedir 'GulfFlow_OneQuarter.nc'])
    ncload([basedir 'GulfFlow_OneTwelfth.nc'])
    gulfflow_onequarter=GulfFlow_OneQuarter;
    gulfflow_onetwelfth=GulfFlow_OneTwelfth;
    clear GulfFlow_OneQuarter GulfFlow_OneTwelfth
    %--------------------------------------------------------------------------
    %loading in and preparing gulfdrifters
    %ncload([basedir 'GulfDriftersAll.nc']);gulfdrifters=GulfDriftersAll;
    %if you're using GulfDriftersOpen
    ncload([basedir 'GulfDriftersOpen.nc']);gulfdrifters=GulfDriftersOpen;
    %if you're using GulfDriftersDWDE
    %ncload([basedir 'GulfDriftersDWDE.nc']);gulfdrifters=GulfDriftersDWDE;
    %--------------------------------------------------------------------------
    load jtopo
catch
    disp('Looks like Matlab can''t file the required data files.  See MAKEFIGS_GULFDRIFTERS for details. ')
end

alpha=1/80;  %stretch factor for quiver arrows
set(0,'DefaultFigureColormap',lansey);
tweakcolorbar=['pos=get(gca,''position'');'...
    'cpos=get(hc,''position'');'...
    'set(hc,''position'',[cpos(1)+0.025 cpos(2)+cpos(4)/2 cpos(3)-0.05 cpos(4)/2]);'...
    'set(gca,''position'',pos)'];

tweakmap_wide=[  'plot([-84.5+1i*21.5 -84.5+1i*22],''k--'',''linewidth'',1);hold on;'...
    'plot([-86.9+1i*21.5 -84.5+1i*21.5],''k--'',''linewidth'',1);'...
    'plot([-81.5+1i*23.04 -81.5+1i*25.82],''w:'');' setstr(10) ...
    'jfig axis|[-98.35 -80.5 18 30.75] latratio|25 topoplot|gray '...
    'eval|topoplot([],[],-5:1:-1,''E'') eval|topoplot([],[],-1/2,''k'') eval|topoplot([],[],-0.005,''E'') ' ...
    'ticksout portrait fontsize|[11 10 10 10] '];

tweakmap=['plot([-84.5+1i*21.5 -84.5+1i*22],''k--'',''linewidth'',1);hold on;'...
    'plot([-86.9+1i*21.5 -84.5+1i*21.5],''k--'',''linewidth'',1);' setstr(10) ...
    'jfig axis|[-98.35 -81.5 18 30.75] latratio|25 topoplot|gray '...
    'eval|topoplot([],[],-5:1:-1,''E'') eval|topoplot([],[],-1/2,''k'') eval|topoplot([],[],-0.005,''E'') ' ...
    'ticksout portrait fontsize|[11 10 10 10] '];

tweakmap_thin=['plot([-84.5+1i*21.5 -84.5+1i*22],''k--'',''linewidth'',0.5);hold on;'...
    'plot([-86.9+1i*21.5 -84.5+1i*21.5],''k--'',''linewidth'',0.5);' setstr(10) ...
    'jfig axis|[-98.35 -81.5 18 30.75] latratio|25 topoplot|gray '...
    'eval|topoplot([],[],-5:1:-1,''0.5E'') eval|topoplot([],[],-1/2,''0.5k'') eval|topoplot([],[],-0.005,''0.5E'') ' ...
    'ticksout portrait fontsize|[11 10 10 10] '];
%\*************************************************************************

if 0
%/*************************************************************************
%Information for the table ...

% 1 latex:6 6
% 2 sculp1:1.5 1.5
% 3 sculp2:1.5 1.5
% 4 gdp:6 6
% 5 hargos:1 1.6136
% 6 aoml:1 1
% 7 sgom:1 1.0083
% 8 ngom:1 3.101
% 9 ocg:1 1.9121
% 10 glad:0.25 0.25
% 11 hercules:0.083333 0.33189
% 12 hgps:1 1.1738
% 13 laser:0.25 0.33095
% 14 dwde:1 1.0007
% 15 splash:0.083333 0.31442

ii=0;
clear driftertype tracking delta
ii=ii+1;driftertype{ii}='WOCE & 7.5$^*$ m';tracking{ii}='Argos';delta{ii}='6.0';%LATEX1
ii=ii+1;driftertype{ii}='CODE & 1 m';tracking{ii}='Argos';delta{ii}='1.5';%SCULP1
ii=ii+1;driftertype{ii}='CODE & 1 m';tracking{ii}='Argos';delta{ii}='1.5';%SCULP2
ii=ii+1;driftertype{ii}='SVP & 15 m';tracking{ii}='Argos';delta{ii}='6.0';%GDP
ii=ii+1;driftertype{ii}='SVP & 15 m';tracking{ii}='Argos';delta{ii}='1.0';%HARGOS
ii=ii+1;driftertype{ii}='CODE & 1 m';tracking{ii}='Argos';delta{ii}='Irreg.';%AOML
ii=ii+1;driftertype{ii}='FHD & 45 m';tracking{ii}='GPS';delta{ii}='1.0';%SGOM
ii=ii+1;driftertype{ii}='FHD & 45 m';tracking{ii}='GPS';delta{ii}='1.0';%NGOM
ii=ii+1;driftertype{ii}='CODE & 1 m';tracking{ii}='Argos';delta{ii}='0.5/1.0';%OCG
ii=ii+1;driftertype{ii}='CODE & 1 m';tracking{ii}='GPS';delta{ii}='0.25';%GLAD
ii=ii+1;driftertype{ii}='Tube & 1 m';tracking{ii}='GPS';delta{ii}='5 min';%Hercules
ii=ii+1;driftertype{ii}='SVP & 15 m';tracking{ii}='GPS';delta{ii}='1.0';%HGPS
ii=ii+1;driftertype{ii}='CARTHE & 1 m';tracking{ii}='GPS';delta{ii}='0.25';%LASER
ii=ii+1;driftertype{ii}='Various & 1 m';tracking{ii}='GPS';delta{ii}='1.5';%DWDE
ii=ii+1;driftertype{ii}='CARTHE & 1 m';tracking{ii}='GPS';delta{ii}='5 min';%SPLASH
ii

for i=1:15
    %i
    use gulfdrifters
    if length(find(source==i))>0
        vindex(time,filled,lat,source==i,1);
        cell2col(time,filled);
        Ltotal=length(find(isfinite(filled)));
        Lfilled=length(find(isfinite(filled)&filled==1));
        
        fillpct=num2str(100*Lfilled./Ltotal,5);
        if length(fillpct)>5, fillpct=fillpct(1:5);end
        
        disp([exp_names(i,:) ' & ' driftertype{i} ' & ' tracking{i} ' &  ' delta{i} '   & '...
            int2str(length(lat)) ' & ' ...
            int2str(Ltotal) ' & ' ...
            fillpct ' & ' ...
            datestr(minmin(time),20) ' & ' ...
            datestr(maxmax(time),20) ' & ' ...
            int2str(mean(cellength(lat)/24)) ' $\pm$ ' int2str(std(cellength(lat)/24,1))  ' & ' ...
            int2str(max(cellength(lat)/24)) ' \\'])
        if i==15, disp('\hline'),end
    end
end

use gulfdrifters
cell2col(time,filled);
Ltotal=length(find(isfinite(filled)));
Lfilled=length(find(isfinite(filled)&filled==1));
fillpct=num2str(100*Lfilled./Ltotal,5);
if length(fillpct)>5, fillpct=fillpct(1:5);end

disp(['All & - & - & -  & 1.0 & '...
    int2str(length(lat)) ' & ' ...
    int2str(Ltotal) ' & ' ...
    fillpct ' & ' ...
    datestr(minmin(time),20) ' & ' ...
    datestr(maxmax(time),20) ' & ' ...
    int2str(mean(cellength(lat)/24)) ' $\pm$ ' int2str(std(cellength(lat)/24,1))  ' & ' ...
    int2str(max(cellength(lat)/24)) ' \\'])

use gulfdrifters
bool=(source~=7)&(source~=8)&(source~=14);
exp_names([7 8 14],:)
vindex(id,lat,time,filled,bool,1);
cell2col(lat,time,filled);
col2cell(lat);
Ltotal=length(find(isfinite(filled)));
Lfilled=length(find(isfinite(filled)&filled==1));
fillpct=num2str(100*Lfilled./Ltotal,5);
if length(fillpct)>5, fillpct=fillpct(1:5);end

disp(['Open & - & - & -  & 1.0 & '...
    int2str(length(id)) ' & ' ...
    int2str(Ltotal) ' & ' ...
    fillpct ' & ' ...
    datestr(minmin(time),20) ' & ' ...
    datestr(maxmax(time),20) ' & ' ...
    int2str(mean(cellength(lat)/24)) ' $\pm$ ' int2str(std(cellength(lat)/24,1))  ' & ' ...
    int2str(max(cellength(lat)/24)) ' \\'])
%\*************************************************************************
end

%/*************************************************************************
%Bathymetry 
figure
axis([-98.35 -80.5 18 30.75])
topoplot([],[-7:1/8:0]), eval(tweakmap_wide),caxis([-4 0])
hc=colorbar('SouthOutside');hc.Label.String='Bottom Depth (km)';eval(tweakcolorbar)
hcl=get(hc,'ticklabels');
for i=1:length(hcl)
    if strcmp(hcl{i}(1),'-')
        hcl{i}=hcl{i}(2:end);
    end
end
set(hc,'ticklabels',hcl)
caxis([-4 0])
fontsize 10 10 10 10        
if strcmp(str,'print')
    jprint(gulfdir,'fig1','-r300')
end
%\*************************************************************************

%/*************************************************************************
%Individual trajectory plots
figure
for i=1:6
    ii=i+0
    subplot(3,2,i)
    use gulfdrifters
    vindex(id,time,lon,lat,filled,source==ii,1);
    for j=1:length(lat)
        lat{j}(filled{j}==1)=inf;
    end
    eval(tweakmap),hold on,cellplot(lon,lat)
    plot(cellfirst(lon),cellfirst(lat),'w.','markersize',10)
    plot(cellfirst(lon),cellfirst(lat),'k.','markersize',8)
    text(-98,30.15,['(' setstr(96+ii) ') ' exp_names(ii,:)])
end
packfig(3,2,'both')
orient tall
fontsize 8 8 8 8
set(gcf,'paperposition',[0.25 0.25 8 9.5])
if strcmp(str,'print')
    jprint(gulfdir,'fig4-1','-r300')
end

figure
for i=1:6
    ii=i+6
    subplot(3,2,i)
    use gulfdrifters
    if length(find(source==ii))>0
        vindex(id,time,lon,lat,filled,source==ii,1);
        for j=1:length(lat)
            lat{j}(filled{j}==1)=inf;
        end
        eval(tweakmap),hold on,cellplot(lon,lat)
        plot(cellfirst(lon),cellfirst(lat),'w.','markersize',10)
        plot(cellfirst(lon),cellfirst(lat),'k.','markersize',8)
    else
        eval(tweakmap)
    end
    text(-98,30.15,['(' setstr(96+ii) ') ' exp_names(ii,:)])
end
packfig(3,2,'both')
orient tall
fontsize 8 8 8 8
set(gcf,'paperposition',[0.25 0.25 8 9.5])
if strcmp(str,'print')
    jprint(gulfdir,'fig4-2','-r300')
end

figure
for i=1:4
    ii=i+12
    subplot(3,2,i)
    use gulfdrifters
    if length(find(source==ii))>0||i==4
        if i<=3
            vindex(id,time,lon,lat,filled,source==ii,1);
        end
        for j=1:length(lat)
            lat{j}(filled{j}==1)=inf;
        end
        eval(tweakmap),hold on,cellplot(lon,lat)
        plot(cellfirst(lon),cellfirst(lat),'w.','markersize',10)
        plot(cellfirst(lon),cellfirst(lat),'k.','markersize',8)
    else
        eval(tweakmap)
    end
    if i<=3
        text(-98,30.15,['(' setstr(96+ii) ') ' exp_names(ii,:)])
    else
        text(-98,30.15,['(' setstr(96+ii) ') All ' ])
    end
end
h=packfig(3,2,'both');
delete(h(5:6))
axes(h(3)),set(gca,'xticklabelmode','auto')
axes(h(4)),set(gca,'xticklabelmode','auto')
orient tall
fontsize 8 8 8 8
set(gcf,'paperposition',[0.25 0.25 8 9.5])
if strcmp(str,'print')
    jprint(gulfdir,'fig4-3','-r300')
end
%\*************************************************************************

%/*************************************************************************
%Velocity vectors with polymap
use gulfflow_onetwelfth
[long,latg]=meshgrid(lon,lat);
cv=vmean(u+1i*v,3)*100;

%load jtopo
jj=find(jtopo.lon>=lon(1)&jtopo.lon<=lon(end)+1/128);
ii=find(jtopo.lat>=lat(1)&jtopo.lat<=lat(end)+1/128);
maxmax(abs(lon-jtopo.lon(jj)'))  %ok, these are basically the same
maxmax(abs(lat-jtopo.lat(ii)))
depth=-jtopo.topo(ii,jj);

total_count=vsum(total_count,3);
index=total_count>0;
%former version
%[ds,xs,ys,zs,ws]=spheresort(latg(index),long(index),cv(index),total_count(index),lat,lon,50); 
%[zhat,beta,aux]=polysmooth(ds,xs,ys,[],zs,ws,0,50);
tic;[zhat,beta,aux]=polymap(long(index),latg(index),[],total_count(index),...
    cv(index),lon,lat,{0,40,2,1},'parallel','sphere');toc
%this will differ very slightly from the published version because the 
%Gaussian kernel is no longer supported

zhat(depth<0)=nan;
%zhat(count<=1)=nan;
%111/12*[1 2 3 4 5]
zhat0=zhat(:,:,1);zhat0(aux(:,:,1)<=4)=nan;

figure
subplot(1,2,1)
%jpcolor(lon,lat,abs(zhat(:,:,1))),hold on,eval(tweakmap_thin)
jpcolor(lon,lat,abs(zhat0)),hold on,eval(tweakmap_thin)
caxis([0 70])
ar=get(gca,'dataaspectratio');
ii=1:3:size(cv,1);jj=1:3:size(cv,2);
long1=vcolon(long(ii,jj));
latg1=vcolon(latg(ii,jj));
zhat1=vcolon(zhat0(ii,jj));
%quiver(long(ii,jj),latg(ii,jj),alpha*real(zhat(ii,jj))*ar(1),alpha*imag(zhat(ii,jj)),0,'k');hold on
quiver([long1(:);-90.2],[latg1(:);20.5],alpha*[real(zhat1);30]*ar(1),alpha*[imag(zhat1);0],0,'k');hold on
text(-98,30.155,'(a)')
text(-90.2,19.95,'30 cm/s')
%--------------------------------------------------------------------------
%Streamline plot
[xg,yg]=meshgrid(lon,lat);
[xo,yo]=latlon2xy(yg,xg,24,-90); 
xy= stream2(xo,yo,real(zhat0),imag(zhat0),xo,yo);
xy=xy(find(cellength(xy)>1));
lats=xy;lons=xy;
for i=1:length(xy)
    [lats{i},lons{i}]=xy2latlon(xy{i}(:,1),xy{i}(:,2),24,-90);
end

subplot(1,2,2)
lons1=lons(1:5:end);
lats1=lats(1:5:end);
%h1=cellplot(lons,lats);eval(tweakmap_thin)
h1=cellplot(lons1,lats1);hold on,eval(tweakmap_thin)
c1=h1;
for i=1:length(h1)
    c1(i)=lons1{i}(1);
end
linecolor(h1,c1,-95,-82,'lansey')
text(-98,30.155,'(b)')
%--------------------------------------------------------------------------
h=packfig(1,2);
axes(h(1))
hc=colorbar('SouthOutside');hc.Label.String='Speed of Mean Flow (cm/s)';eval(tweakcolorbar)
axes(h(2))
hc=colorbar('SouthOutside');hc.Label.String='Initial Longitude for Streamlines';eval(tweakcolorbar)
caxis([0 1]),set(hc,'ticks',[0:2/13:12/13])
%Have to carefully set ticklabels manually 
%1  2  3  4  5  6  7  8  9  10 11 12 13 14
%95 94 93 92 91 90 89 88 87 86 85 84 83 82
%0    2/13  4/13  6/13  8/13 10/13 12/13
set(hc,'TickLabels',{'95','93','91','89','87','85','83'})
orient tall
fontsize 10 10 10 10
set(gcf,'paperposition',[0.25 0.25 11 7.75])
if strcmp(str,'print')
    jprint(gulfdir,'fig3','-r300')
end
%\*************************************************************************

%/*************************************************************************
%Velocity vectors from various sources
foundaviso=true;
try  %Sorry, I'm not redistributing the aviso files
    ncload([basedir 'gom_aviso_madt_means'],'lat','lon','ubar','vbar','uvbar','zetabar','spdbar');
catch
    foundaviso=false;
end
figure
%--------------------------------------------------------------------------
subplot(3,2,1)
if foundaviso
    use gom_aviso_madt_means
    [long,latg]=meshgrid(lon,lat);
    cv=ubar+1i*vbar;
    jpcolor(lon,lat,abs(cv)),hold on,eval(tweakmap_thin),caxis([0 70])
    ar=get(gca,'dataaspectratio');
    quiver([long(:);-90.2],[latg(:);20.5],alpha*[real(cv(:));30]*ar(1),alpha*[imag(cv(:));0],0,'k');hold on
else
    eval(tweakmap_thin)
end
text(-98,30.155,'(a) CMEMS','interpreter','latex')
text(-90.2,19.95,'30 cm/s')
%--------------------------------------------------------------------------
subplot(3,2,2)
load drifter_monthlymeans %Velocity vectors from Laurindo
%https://www.aoml.noaa.gov/phod/gdp/mean_velocity.php
ii=Lat>17.5&Lat<31.25;
jj=Lon>-99&Lon<-81;
lon=Lon(jj);lat=Lat(ii);
u=U(ii,jj)*100;v=V(ii,jj)*100;

[long,latg]=meshgrid(lon,lat);
cv=vmean(u+1i*v,3);
jpcolor(lon,lat,abs(cv)),hold on,eval(tweakmap_thin),caxis([0 70])
ar=get(gca,'dataaspectratio');
h=quiver(long,latg,alpha*real(cv)*ar(1),alpha*imag(cv),0,'k');hold on
text(-98,30.155,'(b) NSVC','interpreter','latex')
%--------------------------------------------------------------------------
subplot(3,2,3)
%from Gulfflow, recreating straight average
use gulfflow_onequarter
[long,latg]=meshgrid(lon,lat);
cv=vsum(total_count.*(u+1i*v),3)./vsum(total_count,3)*100;
depth=-interp2(jtopo.lon,jtopo.lat,jtopo.topo,long,latg);
cv(vsum(total_count,3)<=6)=nan;
cv(depth<0)=nan;

jpcolor(lon,lat,abs(cv)),hold on,eval(tweakmap_thin),caxis([0 70])
ar=get(gca,'dataaspectratio');
quiver(long,latg,alpha*real(cv)*ar(1),alpha*imag(cv),0,'k');hold on
text(-98,30.155,'(c) Bin Average','interpreter','latex')
%--------------------------------------------------------------------------
subplot(3,2,4)
%from Gulfflow
use gulfflow_onequarter
[long,latg]=meshgrid(lon,lat);
cv=vmean(u+1i*v,3)*100;
cv=vmean(cv,3);
depth=interp2(jtopo.lon,jtopo.lat,-jtopo.topo,long,latg);
cv(vsum(total_count,3)<=6)=nan;
cv(depth<0)=nan;

jpcolor(lon,lat,abs(cv)),hold on,eval(tweakmap_thin),caxis([0 70])
ar=get(gca,'dataaspectratio');
quiver(long,latg,alpha*real(cv)*ar(1),alpha*imag(cv),0,'k');hold on
text(-98,30.155,'(d) GulfFlow 1/4$^\circ$','interpreter','latex')
%--------------------------------------------------------------------------

h=packfig(3,2,'both');
delete(h(5:6))
orient tall
fontsize 8 8 8 8
set(gcf,'paperposition',[0.25 0.25 8 9.5])
axes(h(3)),set(gca,'xticklabelmode','auto')
axes(h(4)),set(gca,'xticklabelmode','auto')

axes(h(3))
hc=colorbar('South');hc.Label.String='Speed of Mean Flow (cm/s)';
pos=hc.Position;
hc.Position=[0.32 0.34 0.4 pos(4)/2];
set(hc,'AxisLocation','out')

if strcmp(str,'print')
    jprint(gulfdir,'fig2','-r300')
end
%\*************************************************************************

%/************************************************************************
%data distributions for Gulfflow
%--------------------------------------------------------------------------
figure
subplot(2,1,1)
use gulfflow_onequarter
N=size(total_count,3);
total_count=vsum(total_count>0,3);
total_count(total_count==0)=nan;
%N=length(find(total_count(:)>0));
jpcolor(lon,lat,100*total_count./N),hold on,eval(tweakmap),caxis([0 20])
contour(lon,lat,100*total_count./N,[10 10],'color','k','linewidth',2)
contour(lon,lat,100*total_count./N,[10 10],'color','w','linewidth',1.5)
text(-82.5,30.155,'(a)')
%total possible monthly samples are 
%3D fraction sampled ...  
hc1=colorbar('SouthOutside');hc1.Label.String='Percent of Monthly Time Slices Sampled';
ha=gca;
%--------------------------------------------------------------------------
subplot(2,1,2)
use gulfflow_onequarter
%How many data points are observed, ever, at each latitude
observedever=vsum(total_count>0,3)>0;
N=squeeze(vsum(observedever,1));
%what I want to know is how many January observations are hypothetically 
%possible?  The answer is N*27 
total_count=squeeze(sum(total_count>0,1));  %number of observed latitudes at each longitude and time
%figure,jpcolor(lon,yearfrac(num),total_count')
[yf,mf]=yearfrac(time);
mf=reshape(vshift(mf,-15,1),24,28);
total_count=vsum(reshape(total_count,74,24,28),3)';
m=[1:25]';total_count=[total_count;total_count(1,:)];
total_count=vshift(total_count,-15,1);
%jpcolor(lon,m,log10(total_count./vrep(N*27,25,1))),hold on,caxis([log10(30/max(N)/27) log10(20/100)])
%contour(lon,m,total_count./vrep(N*27,25,1)*100,[8 8],'color','k','linewidth',1.5)
%contour(lon,m,total_count./vrep(N*27,25,1)*100,[8 8],'color','w','linewidth',1)
jpcolor(lon,m,total_count./vrep(N*27,25,1)*100),hold on,caxis([0 20])
%contour(lon,m,total_count./vrep(N*27,25,1)*100,[10 10],'color','k','linewidth',1.5)
%contour(lon,m,total_count./vrep(N*27,25,1)*100,[10 10],'color','w','linewidth',1)
%figure,jpcolor(lon,m,log10(total_count./vrep(N*27,25,1)))
ytick([2:2:24]),ylim([1 25]),xlim([-98.35 -81.5])
text(-82.5,24,'(b)')
set(gca,'yticklabels',['JFMAMJJASOND']');
hc2=colorbar('SouthOutside');hc2.Label.String='Percent of Longitude/Month Bins Sampled';

h=gca;
pos=get(h,'Position');
pos(2)=0.35;
set(h,'position',pos)

%despite valiant efforts, I am unable to thin the colorbars
% posa=get(ha,'position');
% cpos1=get(hc1,'position');
% set(hc1,'position',[cpos1(1) cpos1(2)+cpos1(4)/2 cpos1(3) cpos1(4)/2]);
% set(ha,'position',posa)
% cpos2=get(hc2,'position');
% set(hc2,'position',[cpos2(1) cpos2(2)+cpos2(4)/2 cpos2(3) cpos2(4)/2]);

fontsize 8 8 8 8
set(gcf,'paperposition',[0.25 0.25 4 10])
if strcmp(str,'print')
    jprint(gulfdir,'fig8','-r300')
end
%\************************************************************************

%/************************************************************************
%Figure of drifter temporal ranges
use gulfdrifters
%sty=['TUVWXYZEbgrcmyk'];
%sty=['DTUVWXYZbgrcmyG'];
%sty=['GTUVWXYZbgrcmyC'];
sty=['IbgrcmyCTUVWXYZ'];
%sty=['DTUVWXYZymGbgrk'];

map=zeros(15,3);
figure
for i=1:length(sty)
    h=plot(1:10);
    linestyle(h,sty(i));
    map(i,:)=get(h,'color');
end
map(3,:)=map(3,:).^1.5;
map(5,:)=map(5,:).^(1/8);
map(6,:)=map(6,:).^(1/2);
map(7,:)=map(7,:).^(1/2);
close
colormap(map)

clear longnames
for i=1:size(exp_names,1)
    longnames{i,1}=['(' setstr(96+i) ') ' exp_names(i,:)];
end
    
figure
subplot(1,2,1)
clear h1
for i=1:size(exp_names,1)
    use gulfdrifters
    if length(find(source==i))>0
        vindex(time,source==i,1);
        [xx,sorter]=sort(cellfirst(time));
        time=time(sorter);
        trajnum=time;
        for j=1:length(time)
            trajnum{j}=j+0*trajnum{j};
        end
        h=cellplot(yearfrac(time),trajnum);hold on
        set(h,'color',map(i,:),'linewidth',1.5);
        h1(i)=h(1);
    end
end
ylabel('Trajectory Segments Sorted by Start Date'),xlabel('Time')
axis tight, boxon,vlines([1992:2020],'B:'),outticks,ylim([0 500])
%legend(h1(h1~=0),longnames,'location','North')

cd(gulfdir)
orient portrait
fontsize 12 12 12 12
text(2018.5,480,'(a)')

subplot(1,2,2)
clear h
for i=1:size(exp_names,1)
    use gulfdrifters
    if length(find(source==i))>0
        
        vindex(time,source==i,1);
        sorted=sort(cellength(time),'descend');
        h(i)=plot(sorted/24,1:length(time));
        set(h(i),'color',map(i,:),'linewidth',1.5);hold on
    end
end
ylog,xlabel('Segment Duration (Days)'),ylabel('Number of Segments')
axis tight,outticks
%legend(h,legendnames,'location','North')
legend(h,longnames,'location','East')


text(640,750,'(b)')
h=gca;
pos=get(h,'position');
set(h,'position',[0.525 pos(2:end)])

cd(gulfdir)
orient tall
fontsize 10 10 10 10
set(gcf,'paperposition',[0.25 0.25 11 4])

if strcmp(str,'print')
    jprint(gulfdir,'fig5','-r300')
end
%\************************************************************************

%/*************************************************************************
%Data distributions
use gulfdrifters
source=time;
for i=1:length(time)
    source{i}=gulfdrifters.source(i)+0*source{i};
end

cell2col(time,lon,lat,source,filled);
vindex(time,lon,lat,source,isfinite(filled)&filled==0,1);

[mat,xmid,ymid,numz]=twodstats(lon,lat,yearfrac(time),[-99:1/4:-78],[18:1/4:31]);
%Modal source #
mat=nan*mat;
pctmode=nan*mat;
for i=1:size(mat,1)
    i
    for j=1:size(mat,2)
         bool=abs(lat-ymid(i))<1/8&abs(lon-xmid(j))<1/8;
         mat(i,j)=mode(source(bool));
         %percentage associated with modal source
         pctmode(i,j)=length(find(source(bool)==mat(i,j)))./length(find(source(bool)));
    end
end
figure
%--------------------------------------------------------------------------
h(1)=subplot(1,2,1);
jpcolor(xmid,ymid,log10(numz)),hold on,eval(tweakmap)
caxis([2 4.1])
contour(xmid,ymid,log10(numz),[3 3],'color','k','linewidth',1.5)
contour(xmid,ymid,log10(numz),[3 3],'color','w','linewidth',1)
text(-98,30.155,'(a)')
colormap lansey
%--------------------------------------------------------------------------
h(2)=subplot(1,2,2);
jpcolor(xmid,ymid,mat),hold on,eval(tweakmap)
text(-98,30.155,'(b)')
%contour(xmid,ymid,mat,[6 6],'w')
caxis([0.5 15.5]),colormap(gca,map)
%contour(xmid,ymid,log10(numz),[3 3],'color','k','linewidth',1.5)
%contour(xmid,ymid,log10(numz),[3 3],'color','w','linewidth',1)
%--------------------------------------------------------------------------
h=packfig(1,2);
axes(h(1))
hc=colorbar('SouthOutside');hc.Label.String='Log10 Number of Observations';eval(tweakcolorbar)
hc.Ticks=[2:.25:4.25];
axes(h(2))
hc=colorbar('SouthOutside');hc.Label.String='Most Common Data Source';eval(tweakcolorbar)
hc.Ticks=[1:15]';
%shortnames={'L','S1','S2','G','HA','A','SG','NG','OCG','GL','Her','HG','LSR','D','SPL'}
shortnames=['abcdefghijklmno']';
hc.TickLabels=shortnames;% set(hc,'direction','reverse')
%hc.TickLabelsMode='auto';
%pos1=get(h(1),'position');
%pos=get(h(2),'position');
%posc=get(hc,'position');
%set(h(2),'position',[0.57 pos(2:4)])
%set(h(2),'yticklabelsmode','auto')
orient tall
fontsize 10 10 10 10
%set(gcf,'paperposition',[0.25 0.25 5.75 8])
set(gcf,'paperposition',[0.25 0.25 11 7.75])
if strcmp(str,'print')
    jprint(gulfdir,'fig6','-r300')
end
%return colormap to size 64
set(0,'DefaultFigureColormap',lansey(64));
%\*************************************************************************

%/************************************************************************
%subgridscale variance
use gulfflow_onequarter
cv=(u+1i*v)*100;
sig=sqrt(vstd(real(cv),3).^2+vstd(imag(cv),3).^2);

figure,
subplot(1,2,1),jpcolor(lon,lat,100*sqrt(vmean(epsuu+epsvv,3))),hold on,eval(tweakmap_thin),caxis([10 40])
subplot(1,2,2),jpcolor(lon,lat,sig),hold on,eval(tweakmap_thin),caxis([10 90])
h=packfig(1,2);
axes(h(1))
hc=colorbar('SouthOutside');hc.Label.String='Local Standard Deviation (cm/s)';eval(tweakcolorbar)
text(-98,30.155,'(a)')
axes(h(2))
hc=colorbar('SouthOutside');hc.Label.String='Bulk Standard Deviation (cm/s)';eval(tweakcolorbar)
text(-98,30.155,'(b)')
orient tall
fontsize 10 10 10 10
set(gcf,'paperposition',[0.25 0.25 11 7.75])
if strcmp(str,'print')
    jprint(gulfdir,'fig7','-r300')
end
%\************************************************************************

%--------------------------------------------------------------------------
%additional figures, not used in the paper


if 0
%this one is not shown in the paper but was a step in my processing
%/*************************************************************************
%Figure for mitigating SCULP issues
load gulfdrifters_secondary
use gulfdrifters_secondary
cell2col(cv,source,id,lat,lon,filled);
spd=abs(cv);
%--------------------------------------------------------------------------
%try the median-based robustification of Cleveland
bool=(filled==0);
[mspd,xmid,ymid]=twodstats(lon(bool),lat(bool),abs(cv(bool)),[-99:1/4:-78],[18:1/4:31]);

spdi=interp2(xmid,ymid,mspd,lon,lat,'linear');
spd0=interp2(xmid,ymid,mspd,lon,lat,'nearest');
spdi(~isfinite(spdi))=spd0(~isfinite(spdi));

res=abs(spd-spdi);
[s,xmid,ymid]=twodmed(lon,lat,res,[-99:1/4:-78],[18:1/4:31]);
%[s,xmid,ymid]=twodmed(lon,lat,res,[-99:1/2:-78],[18:1/2:31]);

si=interp2(xmid,ymid,s,lon,lat,'linear');
s0=interp2(xmid,ymid,s,lon,lat,'nearest');
si(~isfinite(si))=s0(~isfinite(si));

eps=abs(frac(res,si));
%delta=squared(1-squared(frac(res,6*si)));
%delta(frac(res,6*si)>1)=0;  %set extreme outliers to zero
%make gulfdrifters delta
%-------------------------------------------------------------------------
epscutoff=7.5;
%figure creation
for k=1:2
    %index into anomalous SCULP data
    %bool=delta==0&(source==9|source==10);
    if k==1
        bool=eps>epscutoff&(source==2|source==3);
    else
        bool=eps>5&eps<epscutoff&(source==2|source==3);
    end
    clear lat1 lon1 spd1
    [N,a,b]=blocknum(bool);
    n=0;
    for i=1:length(a)
        if bool(a(i))
            n=n+1;
            lon1{n,1}=lon(a(i):b(i));
            lat1{n,1}=lat(a(i):b(i));
            spd1{n,1}=spd(a(i):b(i));
        end
    end
    %index into anomalous non-SCULP data
    %bool=delta==0&~(source==9|source==10);
    if k==1
        bool=eps>epscutoff&~(source==2|source==3);
    else
        bool=eps>5&eps<epscutoff&~(source==2|source==3);
    end
    clear lat2 lon2 spd2
    [N,a,b]=blocknum(bool);
    n=0;
    for i=1:length(a)
        if bool(a(i))
            n=n+1;
            lon2{n,1}=lon(a(i):b(i));
            lat2{n,1}=lat(a(i):b(i));
            spd2{n,1}=spd(a(i):b(i));
        end
    end
    
    figure
    subplot(2,2,1),cellplot(lon1,lat1)
    subplot(2,2,2),cellplot(lon2,lat2)
    subplot(2,2,3),
    cell2col(lon1,lat1,spd1);
    [~,index]=sort(spd1);
    scatter(lon1(index),lat1(index),5+0*lat1,spd1(index),'filled');
    subplot(2,2,4),
    cell2col(lon2,lat2,spd2);
    [~,index]=sort(spd2);
    scatter(lon2(index),lat2(index),5+0*lat2,spd2(index),'filled');
    
    for i=1:4
        subplot(2,2,i),
        
        jfig axis|[-98.35 -81.5 18 30.75] latratio|25 topoplot|gray ...
            eval|topoplot([],[],-5:1:-1,'E') eval|topoplot([],[],-1/2,'k') ...
            ticksout portrait fontsize|[11 10 10 10]
        caxis([50 250]),ylim([26 30.75])
        text(-97.8,30,['(' setstr(96+i),') '],'color','w')
    end
    h=packfig(2,2,'both');
    
    axes(h(3))
    hc=colorbar('South');hc.Label.String='Speed (cm/s)';
    %pos=hc.Position;
    hc.Ticks=[50:50:250];
    eval(tweakcolorbar)
    pos=hc.Position;
    set(hc,'position',[pos(1)+0.05/2 0.13 pos(3)-0.1/2 0.02])   
    %hc.Position=[0.32 -0.1 0.4 0.04];
    set(hc,'AxisLocation','in')
    set(hc,'color','m')

    set(gcf,'paperposition', [1 1 10 3.1 ])
    if strcmp(str,'print')
        if k==1
            jprint(gulfdir,'gulfflow_sculpcleanup_epsabovecutoff','-r300')
        else
            jprint(gulfdir,'gulfflow_sculpcleanup_eps5tocutoff','-r300')
        end
    end
    %jprint(gulfdir,'gulfflow_sculpcleanup_after')
    %jprint(gulfdir,'gulfflow_sculpcleanup_interpolated')
end
%\*************************************************************************
end

if 0
%/*************************************************************************
%mean circulation and streamlines without ngom
%use GulfDriftersAll
use GulfDriftersOpen

lono=[-99-1/24:1/12:-80.5];
lato=[18-1/24:1/12:31];
vindex(time,lat,lon,u,v,filled,source,source~=8,1);
gulfflow_nongom=griddrifters(time,lat,lon,u,v,filled,source,lato,lono,1992,7,672);
use gulfflow_nongom
total_count=count;

%The rest is just copied from above
%Velocity vectors with polymap
[long,latg]=meshgrid(lon,lat);
cv=vmean(u+1i*v,3)*100;

%load jtopo
jj=find(jtopo.lon>=lon(1)&jtopo.lon<=lon(end)+1/128);
ii=find(jtopo.lat>=lat(1)&jtopo.lat<=lat(end)+1/128);
maxmax(abs(lon-jtopo.lon(jj)'))  %ok, these are basically the same
maxmax(abs(lat-jtopo.lat(ii)))
depth=-jtopo.topo(ii,jj);

total_count=vsum(total_count,3);
index=total_count>0;
%former version:
%[ds,xs,ys,zs,ws]=spheresort(latg(index),long(index),cv(index),total_count(index),lat,lon,50); 
%[zhat,beta,aux]=polysmooth(ds,xs,ys,[],zs,ws,0,50);
tic;[zhat,beta,aux]=polymap(long(index),latg(index),[],total_count(index),...
    cv(index),lon,lat,{0,40,2,1},'parallel','sphere');toc

zhat(depth<0)=nan;
%zhat(count<=1)=nan;
%111/12*[1 2 3 4 5]
zhat0=zhat(:,:,1);zhat0(aux(:,:,1)<=4)=nan;

figure
subplot(1,2,1)
%jpcolor(lon,lat,abs(zhat(:,:,1))),hold on,eval(tweakmap_thin)
jpcolor(lon,lat,abs(zhat0)),hold on,eval(tweakmap_thin)
caxis([0 70])
ar=get(gca,'dataaspectratio');
ii=1:3:size(cv,1);jj=1:3:size(cv,2);
long1=vcolon(long(ii,jj));
latg1=vcolon(latg(ii,jj));
zhat1=vcolon(zhat0(ii,jj));
%quiver(long(ii,jj),latg(ii,jj),alpha*real(zhat(ii,jj))*ar(1),alpha*imag(zhat(ii,jj)),0,'k');hold on
quiver([long1(:);-90.2],[latg1(:);20.5],alpha*[real(zhat1);30]*ar(1),alpha*[imag(zhat1);0],0,'k');hold on
text(-98,30.155,'(a)')
text(-90.2,19.95,'30 cm/s')
%--------------------------------------------------------------------------
%Streamline plot
[xg,yg]=meshgrid(lon,lat);
[xo,yo]=latlon2xy(yg,xg,24,-90); 
xy= stream2(xo,yo,real(zhat0),imag(zhat0),xo,yo);
xy=xy(find(cellength(xy)>1));
lats=xy;lons=xy;
for i=1:length(xy)
    [lats{i},lons{i}]=xy2latlon(xy{i}(:,1),xy{i}(:,2),24,-90);
end

subplot(1,2,2)
lons1=lons(1:5:end);
lats1=lats(1:5:end);
%h1=cellplot(lons,lats);eval(tweakmap_thin)
h1=cellplot(lons1,lats1);hold on,eval(tweakmap_thin)
c1=h1;
for i=1:length(h1)
    c1(i)=lons1{i}(1);
end
linecolor(h1,c1,-95,-82,'lansey')
text(-98,30.155,'(b)')
%--------------------------------------------------------------------------
h=packfig(1,2);
axes(h(1))
hc=colorbar('SouthOutside');hc.Label.String='Speed of Mean Flow (cm/s)';eval(tweakcolorbar)
axes(h(2))
hc=colorbar('SouthOutside');hc.Label.String='Initial Longitude for Streamlines';eval(tweakcolorbar)
caxis([0 1]),set(hc,'ticks',[0:2/13:12/13])
%Have to carefully set ticklabels manually 
%1  2  3  4  5  6  7  8  9  10 11 12 13 14
%95 94 93 92 91 90 89 88 87 86 85 84 83 82
%0    2/13  4/13  6/13  8/13 10/13 12/13
set(hc,'TickLabels',{'95','93','91','89','87','85','83'})
orient tall
fontsize 10 10 10 10
set(gcf,'paperposition',[0.25 0.25 11 7.75])
if strcmp(str,'print')
    jprint(gulfdir,'gulfdrifters_circulation_nongom','-r300')
end
%\*************************************************************************
end

if 0
%/*************************************************************************
%Figure 2d, but for drogued/undrogued/unknown/all data

lono=[-99:1/4:-80.5];
lato=[18:1/4:31];

use GulfDriftersAll
cell2col(time,lat,lon,u,v,filled,drogue);
vindex(time,lat,lon,u,v,filled,drogue==1|isnan(drogue),1);
col2cell(time,lat,lon,u,v,filled,drogue);
gulfflow_drogued=griddrifters(time,lat,lon,u,v,filled,source,lato,lono,1992,7,672);
use GulfDriftersAll
cell2col(time,lat,lon,u,v,filled,drogue);
vindex(time,lat,lon,u,v,filled,drogue==0|isnan(drogue),1);
col2cell(time,lat,lon,u,v,filled,drogue);
gulfflow_undrogued=griddrifters(time,lat,lon,u,v,filled,source,lato,lono,1992,7,672);
use GulfDriftersAll
cell2col(time,lat,lon,u,v,filled,drogue);
vindex(time,lat,lon,u,v,filled,drogue==inf|isnan(drogue),1);
col2cell(time,lat,lon,u,v,filled,drogue);
gulfflow_unknown=griddrifters(time,lat,lon,u,v,filled,source,lato,lono,1992,7,672);
%--------------------------------------------------------------------------
figure
subplot(3,2,1)
use gulfflow_drogued
[long,latg]=meshgrid(lon,lat);
cv=vmean(u+1i*v,3)*100;
cv=vmean(cv,3);
depth=interp2(jtopo.lon,jtopo.lat,-jtopo.topo,long,latg);
%cv(vsum(total_count,3)<=6)=nan;
cv(depth<0)=nan;

jpcolor(lon,lat,abs(cv)),hold on,eval(tweakmap_thin),caxis([0 70])
ar=get(gca,'dataaspectratio');
quiver(long,latg,alpha*real(cv)*ar(1),alpha*imag(cv),0,'k');hold on
text(-98,30.155,'(a) GulfFlow 1/4$^\circ$ drogued','interpreter','latex')
%--------------------------------------------------------------------------
subplot(3,2,2)
use gulfflow_undrogued
[long,latg]=meshgrid(lon,lat);
cv=vmean(u+1i*v,3)*100;
cv=vmean(cv,3);
depth=interp2(jtopo.lon,jtopo.lat,-jtopo.topo,long,latg);
%cv(vsum(total_count,3)<=6)=nan;
cv(depth<0)=nan;

jpcolor(lon,lat,abs(cv)),hold on,eval(tweakmap_thin),caxis([0 70])
ar=get(gca,'dataaspectratio');
quiver(long,latg,alpha*real(cv)*ar(1),alpha*imag(cv),0,'k');hold on
text(-98,30.155,'(b) GulfFlow 1/4$^\circ$ undrogued','interpreter','latex')
%--------------------------------------------------------------------------
subplot(3,2,3)
use gulfflow_unknown
[long,latg]=meshgrid(lon,lat);
cv=vmean(u+1i*v,3)*100;
cv=vmean(cv,3);
depth=interp2(jtopo.lon,jtopo.lat,-jtopo.topo,long,latg);
%cv(vsum(total_count,3)<=6)=nan;
cv(depth<0)=nan;

jpcolor(lon,lat,abs(cv)),hold on,eval(tweakmap_thin),caxis([0 70])
ar=get(gca,'dataaspectratio');
quiver(long,latg,alpha*real(cv)*ar(1),alpha*imag(cv),0,'k');hold on
text(-98,30.155,'(d) GulfFlow 1/4$^\circ$ unknown','interpreter','latex')
%--------------------------------------------------------------------------
subplot(3,2,4)
use gulfflow_onequarter
[long,latg]=meshgrid(lon,lat);
cv=vmean(u+1i*v,3)*100;
cv=vmean(cv,3);
depth=interp2(jtopo.lon,jtopo.lat,-jtopo.topo,long,latg);
%cv(vsum(total_count,3)<=6)=nan;
cv(depth<0)=nan;

jpcolor(lon,lat,abs(cv)),hold on,eval(tweakmap_thin),caxis([0 70])
ar=get(gca,'dataaspectratio');
quiver(long,latg,alpha*real(cv)*ar(1),alpha*imag(cv),0,'k');hold on
text(-98,30.155,'(d) GulfFlow 1/4$^\circ$ all','interpreter','latex')
%--------------------------------------------------------------------------

h=packfig(3,2,'both');
delete(h(5:6))
orient tall
fontsize 8 8 8 8
set(gcf,'paperposition',[0.25 0.25 8 9.5])
axes(h(3)),set(gca,'xticklabelmode','auto')
axes(h(4)),set(gca,'xticklabelmode','auto')

axes(h(3))
hc=colorbar('South');hc.Label.String='Speed of Mean Flow (cm/s)';
pos=hc.Position;
hc.Position=[0.32 0.34 0.4 pos(4)/2];
set(hc,'AxisLocation','out')

if strcmp(str,'print')
    jprint(gulfdir,'gulfdrifters_meanvelocities_droguestatus','-r300')
end
%\*************************************************************************
end

if 0
%/*************************************************************************
%Figure 2d, but for difference drogue depths

lono=[-99:1/4:-80.5];
lato=[18:1/4:31];

use GulfDriftersAll
vindex(time,lat,lon,u,v,filled,(type~=1&type~=3&type~=4),1);
gulfflow_shallow=griddrifters(time,lat,lon,u,v,filled,source,lato,lono,1992,7,672);
use GulfDriftersAll
vindex(time,lat,lon,u,v,filled, type==3|type==1,1);
gulfflow_medium=griddrifters(time,lat,lon,u,v,filled,source,lato,lono,1992,7,672);
use GulfDriftersAll
vindex(time,lat,lon,u,v,filled,type==4,1);
gulfflow_deep=griddrifters(time,lat,lon,u,v,filled,source,lato,lono,1992,7,672);
%--------------------------------------------------------------------------
figure
subplot(3,2,1)
use gulfflow_shallow
[long,latg]=meshgrid(lon,lat);
cv=vmean(u+1i*v,3)*100;
cv=vmean(cv,3);
depth=interp2(jtopo.lon,jtopo.lat,-jtopo.topo,long,latg);
%cv(vsum(total_count,3)<=6)=nan;
cv(depth<0)=nan;

jpcolor(lon,lat,abs(cv)),hold on,eval(tweakmap_thin),caxis([0 70])
ar=get(gca,'dataaspectratio');
quiver(long,latg,alpha*real(cv)*ar(1),alpha*imag(cv),0,'k');hold on
text(-98,30.155,'(a) GulfFlow 1/4$^\circ$ shallow','interpreter','latex')
%--------------------------------------------------------------------------
subplot(3,2,2)
use gulfflow_medium
[long,latg]=meshgrid(lon,lat);
cv=vmean(u+1i*v,3)*100;
cv=vmean(cv,3);
depth=interp2(jtopo.lon,jtopo.lat,-jtopo.topo,long,latg);
%cv(vsum(total_count,3)<=6)=nan;
cv(depth<0)=nan;

jpcolor(lon,lat,abs(cv)),hold on,eval(tweakmap_thin),caxis([0 70])
ar=get(gca,'dataaspectratio');
quiver(long,latg,alpha*real(cv)*ar(1),alpha*imag(cv),0,'k');hold on
text(-98,30.155,'(b) GulfFlow 1/4$^\circ$ medium','interpreter','latex')
%--------------------------------------------------------------------------
subplot(3,2,3)
use gulfflow_deep
[long,latg]=meshgrid(lon,lat);
cv=vmean(u+1i*v,3)*100;
cv=vmean(cv,3);
depth=interp2(jtopo.lon,jtopo.lat,-jtopo.topo,long,latg);
%cv(vsum(total_count,3)<=6)=nan;
cv(depth<0)=nan;

jpcolor(lon,lat,abs(cv)),hold on,eval(tweakmap_thin),caxis([0 70])
ar=get(gca,'dataaspectratio');
quiver(long,latg,alpha*real(cv)*ar(1),alpha*imag(cv),0,'k');hold on
text(-98,30.155,'(d) GulfFlow 1/4$^\circ$ deep','interpreter','latex')
%--------------------------------------------------------------------------
subplot(3,2,4)
use gulfflow_onequarter
[long,latg]=meshgrid(lon,lat);
cv=vmean(u+1i*v,3)*100;
cv=vmean(cv,3);
depth=interp2(jtopo.lon,jtopo.lat,-jtopo.topo,long,latg);
%cv(vsum(total_count,3)<=6)=nan;
cv(depth<0)=nan;

jpcolor(lon,lat,abs(cv)),hold on,eval(tweakmap_thin),caxis([0 70])
ar=get(gca,'dataaspectratio');
quiver(long,latg,alpha*real(cv)*ar(1),alpha*imag(cv),0,'k');hold on
text(-98,30.155,'(d) GulfFlow 1/4$^\circ$ all','interpreter','latex')
%--------------------------------------------------------------------------

h=packfig(3,2,'both');
delete(h(5:6))
orient tall
fontsize 8 8 8 8
set(gcf,'paperposition',[0.25 0.25 8 9.5])
axes(h(3)),set(gca,'xticklabelmode','auto')
axes(h(4)),set(gca,'xticklabelmode','auto')

axes(h(3))
hc=colorbar('South');hc.Label.String='Speed of Mean Flow (cm/s)';
pos=hc.Position;
hc.Position=[0.32 0.34 0.4 pos(4)/2];
set(hc,'AxisLocation','out')

if strcmp(str,'print')
    jprint(gulfdir,'gulfdrifters_meanvelocities_droguedepths','-r300')
end
%\*************************************************************************
end

if 0
%/*************************************************************************
%What is the relative count of drogued/undrogued/unknown drifters?
figure
count=[]%this is just to keep Matlab from complaining
use gulfflow_onequarter
num_drogued=vsum(vsum(count(:,:,:,1:15),4),3);
num_undrogued=vsum(vsum(count(:,:,:,16:30),4),3);
num_unknown=vsum(vsum(count(:,:,:,31:45),4),3);
%--------------------------------------------------------------------------
subplot(3,2,1);
jpcolor(lon,lat,log10(num_drogued)),hold on,eval(tweakmap)
caxis([2 4.1])
contour(lon,lat,log10(num_drogued),[3 3],'color','k','linewidth',1.5)
contour(lon,lat,log10(num_drogued),[3 3],'color','w','linewidth',1)
text(-98,30.155,'(a)')
colormap lansey
subplot(3,2,2);
jpcolor(lon,lat,log10(num_undrogued)),hold on,eval(tweakmap)
caxis([2 4.1])
contour(lon,lat,log10(num_undrogued),[3 3],'color','k','linewidth',1.5)
contour(lon,lat,log10(num_undrogued),[3 3],'color','w','linewidth',1)
text(-98,30.155,'(b)')
colormap lansey
subplot(3,2,3);
jpcolor(lon,lat,log10(num_unknown)),hold on,eval(tweakmap)
caxis([2 4.1])
contour(lon,lat,log10(num_unknown),[3 3],'color','k','linewidth',1.5)
contour(lon,lat,log10(num_unknown),[3 3],'color','w','linewidth',1)
text(-98,30.155,'(c)')
colormap lansey
subplot(3,2,4);
jpcolor(lon,lat,log10(vsum(total_count,3))),hold on,eval(tweakmap)
caxis([2 4.1])
contour(lon,lat,log10(vsum(total_count,3)),[3 3],'color','k','linewidth',1.5)
contour(lon,lat,log10(vsum(total_count,3)),[3 3],'color','w','linewidth',1)
text(-98,30.155,'(d)')
colormap lansey

h=packfig(3,2,'both');
delete(h(5:6))
orient tall
fontsize 8 8 8 8
set(gcf,'paperposition',[0.25 0.25 8 9.5])
axes(h(3)),set(gca,'xticklabelmode','auto')
axes(h(4)),set(gca,'xticklabelmode','auto')

axes(h(3))
hc=colorbar('South');hc.Label.String='Log10 Number of Observations';
pos=hc.Position;
hc.Position=[0.32 0.34 0.4 pos(4)/2];
set(hc,'AxisLocation','out')

if strcmp(str,'print')
    jprint(gulfdir,'gulfdrifters_distribution_droguestatus','-r300')
end

%total percentage of these three classes?
use gulfflow_onequarter
N=sum(total_count(:));
pct_drogued=100*vsum(vsum(vsum(vsum(count(:,:,:,1:15),4),3),2),1)/N
pct_undrogued=100*vsum(vsum(vsum(vsum(count(:,:,:,16:30),4),3),2),1)/N
pct_unknown=100*vsum(vsum(vsum(vsum(count(:,:,:,31:45),4),3),2),1)/N
%\*************************************************************************
end

if 0
%/*************************************************************************
%What is the relative count of shallow/intermediate/deep drifters?
use gulfflow_shallow
num_shallow=vsum(count,3);
use gulfflow_medium
num_medium=vsum(count,3);
use gulfflow_deep
num_deep=vsum(count,3);
%--------------------------------------------------------------------------
subplot(3,2,1);
jpcolor(lon,lat,log10(num_shallow)),hold on,eval(tweakmap)
caxis([2 4.1])
contour(lon,lat,log10(num_shallow),[3 3],'color','k','linewidth',1.5)
contour(lon,lat,log10(num_shallow),[3 3],'color','w','linewidth',1)
text(-98,30.155,'(a)')
colormap lansey
subplot(3,2,2);
jpcolor(lon,lat,log10(num_medium)),hold on,eval(tweakmap)
caxis([2 4.1])
contour(lon,lat,log10(num_medium),[3 3],'color','k','linewidth',1.5)
contour(lon,lat,log10(num_medium),[3 3],'color','w','linewidth',1)
text(-98,30.155,'(b)')
colormap lansey
subplot(3,2,3);
jpcolor(lon,lat,log10(num_deep)),hold on,eval(tweakmap)
caxis([2 4.1])
contour(lon,lat,log10(num_deep),[3 3],'color','k','linewidth',1.5)
contour(lon,lat,log10(num_deep),[3 3],'color','w','linewidth',1)
text(-98,30.155,'(c)')
colormap lansey
subplot(3,2,4);
use gulfflow_onequarter
jpcolor(lon,lat,log10(vsum(total_count,3))),hold on,eval(tweakmap)
caxis([2 4.1])
contour(lon,lat,log10(vsum(total_count,3)),[3 3],'color','k','linewidth',1.5)
contour(lon,lat,log10(vsum(total_count,3)),[3 3],'color','w','linewidth',1)
text(-98,30.155,'(d)')
colormap lansey

h=packfig(3,2,'both');
delete(h(5:6))
orient tall
fontsize 8 8 8 8
set(gcf,'paperposition',[0.25 0.25 8 9.5])
axes(h(3)),set(gca,'xticklabelmode','auto')
axes(h(4)),set(gca,'xticklabelmode','auto')

axes(h(3))
hc=colorbar('South');hc.Label.String='Log10 Number of Observations';
pos=hc.Position;
hc.Position=[0.32 0.34 0.4 pos(4)/2];
set(hc,'AxisLocation','out')

if strcmp(str,'print')
    jprint(gulfdir,'gulfdrifters_distribution_droguedepths','-r300')
end
%\*************************************************************************
end

if 0
%/************************************************************************
%seasonal cycle of subgrid-scale variance
use gulfflow_onequarter
[yf,mf]=yearfrac(time);
clear bool
%bool{1}=(floor(mf)==12)|(floor(mf)==1)|(floor(mf)==2);
%bool{2}=(floor(mf)==3)|(floor(mf)==4)|(floor(mf)==5);
%bool{3}=(floor(mf)==6)|(floor(mf)==7)|(floor(mf)==8);
%bool{4}=(floor(mf)==9)|(floor(mf)==10)|(floor(mf)==11);
bool{1}=(floor(mf)==1)|(floor(mf)==2)|(floor(mf)==3);
bool{2}=(floor(mf)==4)|(floor(mf)==5)|(floor(mf)==6);
bool{3}=(floor(mf)==7)|(floor(mf)==8)|(floor(mf)==9);
bool{4}=(floor(mf)==10)|(floor(mf)==11)|(floor(mf)==12);
for i=1:4
    epsac(:,:,i)=vmean(epsuu(:,:,bool{i})+epsvv(:,:,bool{i}),3);
end
figure
for i=1:4
    subplot(3,2,i),jpcolor(lon,lat,100*sqrt(epsac(:,:,i))),hold on,eval(tweakmap_thin),caxis([10 40])
    text(-98,30.155,['(' setstr(96+i) ')'])
end

h=packfig(3,2,'both');
delete(h(5:6))
orient tall
fontsize 8 8 8 8
set(gcf,'paperposition',[0.25 0.25 8 9.5])
axes(h(3)),set(gca,'xticklabelmode','auto')
axes(h(4)),set(gca,'xticklabelmode','auto')

axes(h(3))
hc=colorbar('South');hc.Label.String='Local Standard Deviation (cm/s)';
pos=hc.Position;
hc.Position=[0.32 0.34 0.4 pos(4)/2];
set(hc,'AxisLocation','out')

if strcmp(str,'print')
    jprint(gulfdir,'gulfdrifters_gridvariance_annualcycle','-r300')
end
%\************************************************************************
end

if 0
%/*************************************************************************
%computing error for table 2 in Lilly and Perez-Brunius (2021a)
%these sources files aren't redistributed, but are available upon request
load gulfflow_cmems
load gulfflow_hycom
load gulfflow_nemo
load gulfflow_roms
clear cvo cv1 cv2 latcell loncell cvcell
%--------------------------------------------------------------------------
%truth in cm/s
cvcell{1}=100*vmean(ncread([readdir 'gom_aviso_madt.nc'],'u'),3)+...
    1i*100*vmean(ncread([readdir 'gom_aviso_madt.nc'],'v'),3);
latcell{1}=ncread([readdir 'gom_aviso_madt.nc'],'lat');
loncell{1}=ncread([readdir 'gom_aviso_madt.nc'],'lon');
cvcell{2}=ncread([readdir 'hycom_surface.nc'],'ubar')+...
    1i*ncread([readdir 'hycom_surface.nc'],'vbar');
latcell{2}=ncread([readdir 'hycom_surface.nc'],'lat');
loncell{2}=ncread([readdir 'hycom_surface.nc'],'lon');
cvcell{3}=ncread([readdir 'nemo_surface.nc'],'ubar')+...
    1i*ncread([readdir 'nemo_surface.nc'],'vbar');
latcell{3}=ncread([readdir 'nemo_surface.nc'],'lat');
loncell{3}=ncread([readdir 'nemo_surface.nc'],'lon');
cvcell{4}=ncread([readdir 'roms_surface.nc'],'ubar')+...
    1i*ncread([readdir 'roms_surface.nc'],'vbar');
latcell{4}=ncread([readdir 'roms_surface.nc'],'lat');
loncell{4}=ncread([readdir 'roms_surface.nc'],'lon');
use gulfflow
for i=1:4
    [long,latg]=meshgrid(loncell{i},latcell{i});
    [longo,latgo]=meshgrid(lon,lat);
    cvo(:,:,i)=interp2(long,latg,cvcell{i},longo,latgo);
end

%--------------------------------------------------------------------------
%single averaging
use gulfflow_cmems
cv=100*(u+1i*v);cv1(:,:,1)=vsum(count(:,:,1:2:end).*cv(:,:,1:2:end),3)./vsum(count(:,:,1:2:end),3);
use gulfflow_hycom
cv=100*(u+1i*v);cv1(:,:,2)=vsum(count(:,:,1:2:end).*cv(:,:,1:2:end),3)./vsum(count(:,:,1:2:end),3);
use gulfflow_nemo
cv=100*(u+1i*v);cv1(:,:,3)=vsum(count(:,:,1:2:end).*cv(:,:,1:2:end),3)./vsum(count(:,:,1:2:end),3);
use gulfflow_roms
cv=100*(u+1i*v);cv1(:,:,4)=vsum(count(:,:,1:2:end).*cv(:,:,1:2:end),3)./vsum(count(:,:,1:2:end),3);
%--------------------------------------------------------------------------
%double averaging
use gulfflow_cmems
cv=100*(u(:,:,1:2:end)+1i*v(:,:,1:2:end));cv2(:,:,1)=vmean(cv,3);
use gulfflow_hycom
cv=100*(u(:,:,1:2:end)+1i*v(:,:,1:2:end));cv2(:,:,2)=vmean(cv,3);
use gulfflow_nemo
cv=100*(u(:,:,1:2:end)+1i*v(:,:,1:2:end));cv2(:,:,3)=vmean(cv,3);
use gulfflow_roms
cv=100*(u(:,:,1:2:end)+1i*v(:,:,1:2:end));cv2(:,:,4)=vmean(cv,3);
%--------------------------------------------------------------------------

for i=1:4
    spd1(i)=sqrt(vmean(vcolon(squared(cvo(:,:,i))),1));
    err1(i)=sqrt(vmean(vcolon(squared(cvo(:,:,i)-cv1(:,:,i))),1));
    err2(i)=sqrt(vmean(vcolon(squared(cvo(:,:,i)-cv2(:,:,i))),1));
end

%numbers for table 2
[spd1' err1' err2' (err1-err2)'./err1'*100 100*err2'./spd1']
%\*************************************************************************


