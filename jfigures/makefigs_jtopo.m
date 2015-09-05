function[]=makefigs_jtopo
%MAKEFIGS_JTOPO  Makes a sample figure for ABOUT_JTOPO.

figure
load jtopo,use jtopo
jpcolor(lon,lat,topo),latratio(30),colormap lansey,caxis([-6 3.5])
xlabel('Longitude'),ylabel('Latitude')
title('One-sixth degree global topography from JTOPO.MAT')
set(gcf,'paperposition',[1 1 12 8])

%To print
if 0
    currentdir=pwd;
    cd([whichdir('jlab_license') '/figures'])
    print -dpng jtopo
    crop jtopo.png
    cd(currentdir)
end
