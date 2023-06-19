function[]=makefigs_polymap2
%MAKEFIGS_POLYMAP2  Makes a sample figure for POLYMAP using TPJAOS.MAT.
    
%This may take a while...
if isempty(which('tpjaos.mat'))
    disp('Sorry, POLYMAP can''t find TPJAOS.MAT.')
    return
end

load tpjaos,use tpjaos

%Standard deviation
z=vstd(tpjaos.ssh,3);

%Map standard deviation into sorted grid format using grid.index

lono=(-180:1:180)';
lato=(-66:1:66)';
tic;[ds,xs,ys,zs]=pm_sort(lon,lat,z,lono,lato,200,'sphere','parallel');toc

%70 points requires at most about 180 km bandwidth

%Zeroth-order fit with 70 points included
[zs,amat,xmat,wmat,Hmat,C]=pm_weight(ds,xs,ys,zs,[],[],{0,200,2,1},'population',70);
zhat0=pm_apply(zs,amat,xmat);

%First-order fit with 70 points included
[zs,amat,xmat,wmat,Hmat,C]=pm_weight(ds,xs,ys,zs,[],[],{1,200,2,1},'population',70);
zhat1=pm_apply(zs,amat,xmat);

%When the first-order fit fails (e.g. near coasts), use the zeroth-order
zhat1(isnan(zhat1))=zhat0(isnan(zhat1));

figure
contourf(lono,lato,zhat1,[0:1/2:50]),nocontours,caxis([2.5 45]),colormap lansey
latratio(30),ylim([-66 66])
title('Standard Deviation of SSH from TPJAOS.MAT, mapped using POLYMAP')
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
    print -dpng polymap
    crop polymap.png
    cd(currentdir)
end
