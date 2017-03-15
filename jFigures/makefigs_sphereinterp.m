function[]=makefigs_sphereinterp
%MAKEFIGS_SPHEREINTERP  Makes two sample figure for SPHEREINTERP.

load goldsnapshot, use goldsnapshot

%You can run this commented-out part yourself but it takes a long time to
%generate the initial fields, so these are included in file GOLDSNAPSHOT
%tic;[dx,dy,index,bool,C]=sphereinterp(geolat,geolon,lat,lon,'parallel','periodic');toc;

tic;[ssh,ssh1]=sphereinterp(dx,dy,index,bool,ssho);toc

figure
subplot(1,2,1),jpcolor(lono,lato,ssho),caxis([-250 250])
title('Snapshot of SSH on GOLD''s tripolar grid')
subplot(1,2,2),jpcolor(lon,lat,ssh),caxis([-250 250])
title('Snapshot of SSH on regular lat/lon grid')
fontsize 14 12 12 12
 
figure
[fx,fy]=spheregrad(lat,lon,ssh);
jpcolor(lon,lat,sqrt(fx.^2+fy.^2))
caxis([0 1]/1000)
topoplot continents
latratio(40),ylim([-78 90]),axis off
fontsize 14 12 12 12

%To print
if 0
    currentdir=pwd;
    cd([whichdir('jlab_license') '/figures'])
    print -dpng sphereinterp1.png
    crop('sphereinterp1.png')
  
    close
    set(gcf,'paperposition',[1 1 15 5])
    print -dpng sphereinterp2.png
    crop('sphereinterp2.png')
    cd(currentdir) 
end


