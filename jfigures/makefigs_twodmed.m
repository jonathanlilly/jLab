function[]=makefigs_twodmed
%MAKEFIGS_TWODMED  Makes a sample figure for TWODMED.

if isempty(which('drifters.mat'))
    disp('Sorry, TWODMED can''t find DRIFTERS.MAT.')
    return
end

%This make take a few minutes...

load drifters,use drifters
%Decimate to speed things up
vindex(lat,lon,cv,1:10:length(lat),1);

[mat,xmid,ymid]=twodmed(lon,lat,cellabs(cv),-180.5:180.5,-89.5:89.5);
figure,jpcolor(xmid,ymid,mat)
axis([-180 180 -70 90]),latratio(30),topoplot
xlabel('Longitude'),ylabel('Latitude')
title('Median Speed from the Global Surface Drifter Dataset')
caxis([8 60]),h=colorbar('EastOutside');colormap lansey
h.Label.String='Median Speed (cm/s)';

%To print
if 0
    currentdir=pwd;
    cd([whichdir('jlab_license') '/figures'])
    print -dpng twodmed
    crop twodmed.png
    cd(currentdir)
end

