function[]=makefigs_twodhist
%MAKEFIGS_TWODHIST  Makes a sample figure for TWODHIST.

if isempty(which('drifters.mat'))
    disp('Sorry, TWODHIST can''t find DRIFTERS.MAT.')
    return
end

%This make take a few minutes...

load drifters,use drifters
%Decimate to speed things up
vindex(lat,lon,1:10:length(lat),1);

lon=cell2col(lon);lat=cell2col(lat);
tic;[mat,xmid,ymid]=twodhist(lon,lat,-180.5:180.5,-89.5:89.5);etime1=toc;
tic;[mat,xmid,ymid]=twodhist(lon,lat,-180.5:180.5,-89.5:89.5,'jlab');etime1=toc;

figure,jpcolor(xmid,ymid,log10(mat))
axis([-180 180 -70 90]),latratio(30),topoplot,colormap lansey
xlabel('Longitude'),ylabel('Latitude')
title('Histogram of Data from the Global Surface Drifter Dataset')
%caxis([1.5 3.5]),h=colorbar('EastOutside');
caxis([1 3]),h=colorbar('EastOutside');
h.Label.String='Log10 Number of Observations';

%To print
if 0
    currentdir=pwd;
    cd([whichdir('jlab_license') '/figures'])
    print -dpng twodhist
    crop twodhist.png
    cd(currentdir)
end
