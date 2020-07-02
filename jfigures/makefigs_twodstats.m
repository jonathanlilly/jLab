function[]=makefigs_twodstats
%MAKEFIGS_TWODSTATS  Makes a sample figure for TWODSTATS.

if isempty(which('drifters.mat'))
    disp('Sorry, TWODSTATS can''t find DRIFTERS.MAT.')
    return
end

%This make take a few minutes...

load drifters,use drifters
cv=cellpair(u,v);
vindex(lat,lon,cv,1:10:length(lat),1);

tic;[mat,xmid,ymid]=twodstats(lon,lat,cellabs(cv),-180.5:180.5,-89.5:89.5);etime1=toc;

figure,jpcolor(xmid,ymid,mat)
axis([-180 180 -70 90]),latratio(30),topoplot
xlabel('Longitude'),ylabel('Latitude')
title('Mean Speed from the Global Surface Drifter Dataset')
caxis([8 60]),h=colorbar('EastOutside');colormap lansey
h.Label.String='Mean Speed (cm/s)';

%To print
if 0
    currentdir=pwd;
    cd([whichdir('jlab_license') '/figures'])
    print -dpng twodstats
    crop twodstats.png
    cd(currentdir)
end