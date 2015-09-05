function[]=makefigs_polysmooth2
%MAKEFIGS_POLYSMOOTH2  Makes a sample figure for POLYSMOOTH using ALONGTRACK.MAT.
    
%This may take a while...
if isempty(which('alongtrack.mat'))
    disp('Sorry, POLYSMOOTH can''t find ALONGTRACK.MAT.')
    return
end

load alongtrack,use alongtrack

%Standard deviation
z=vstd(alongtrack.ssh,3);

%Map standard deviation into sorted grid format using grid.index
use grid
zs=nan*zeros(size(xs));
zs(~isnan(index))=z(index(~isnan(index)));

%Zeroth-order fit with 100 points included
[zhat0,weight,beta,b]=polysmooth(ds,xs,ys,zs,100,0,'sphere','variable');

%First-order fit with 70 points included
[zhat1,weight,beta,b]=polysmooth(ds,xs,ys,zs,70,1,'sphere','variable');

%When the first-order fit fails (e.g. near coasts), use the zeroth-order
zhat1(isnan(zhat1))=zhat0(isnan(zhat1));

figure
contourf(lono,lato,zhat1,[0:1/2:50]),nocontours,caxis([2.5 45]),colormap lansey
latratio(30),ylim([-66 66])
title('Standard Deviation of SSH from ALONGTRACK.MAT, mapped using POLYSMOOTH')
hc=colorbar('EastOutside');hc.Label.String='SSH Standard Deviation (cm)';
topoplot continents
set(gcf,'paperposition',[1 1 12 8])

use alongtrack
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
