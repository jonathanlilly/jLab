function[]=makefigs_polysmooth2
%MAKEFIGS_POLYSMOOTH2  Makes a sample figure for POLYSMOOTH using TPJAOS.MAT.
    
%This may take a while...
if isempty(which('tpjaos.mat'))
    disp('Sorry, POLYSMOOTH can''t find TPJAOS.MAT.')
    return
end

load tpjaos,use tpjaos

%Standard deviation
z=vstd(tpjaos.ssh,3);

%Map standard deviation into sorted grid format using grid.index

load tpjaos
use tpjaos
lono=(-180:1:180)';
lato=(-66:1:66)';
tic;[ds,xs,ys,index]=spheresort(lat,lon,1:length(lat(:)),lato,lono,200,'parallel');toc
%tic;[ds,xs,ys,index]=spheresort(lat,lon,1:length(lat(:)),lato,lono,200);toc

make grid lato lono ds xs ys index


use grid
zs=nan*zeros(size(xs));
zs(~isnan(index))=z(index(~isnan(index)));

%Zeroth-order fit with 100 points included
tic;[zhat0,weight0,beta0,b0]=polysmooth(ds,xs,ys,zs,100,0,'sphere','variable');toc

%First-order fit with 70 points included
[zhat1,weight1,beta1,b1]=polysmooth(ds,xs,ys,zs,70,1,'sphere','variable');

%When the first-order fit fails (e.g. near coasts), use the zeroth-order
zhat1(isnan(zhat1))=zhat0(isnan(zhat1));

zhat2=zhat1;
zhat2(isnan(zhat1))=zhat0(isnan(zhat1));

figure
contourf(lono,lato,zhat1,[0:1/2:50]),nocontours,caxis([2.5 45]),colormap lansey
latratio(30),ylim([-66 66])
title('Standard Deviation of SSH from TPJAOS.MAT, mapped using POLYSMOOTH')
hc=colorbar('EastOutside');hc.Label.String='SSH Standard Deviation (cm)';
topoplot continents
set(gcf,'paperposition',[1 1 12 8])

use tpjaos
lon(abs(vdiff(lon,1))>90)=nan;
hold on,plot(lon,lat,'w:','linewidth',0.1)

%To print

if 0
    currentdir=pwd;
    cd([whichdir('jlab_license') '/figures'])
    print -dpng polysmooth
    crop polysmooth.png
    cd(currentdir)
end
