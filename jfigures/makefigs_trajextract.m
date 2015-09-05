function[]=makefigs_trajextract
%MAKEFIGS_TRAJEXTRACT  Makes a sample figure for TRAJEXTRACT.

load ebasnfloats
use ebasnfloats
region=[-30 -21 24 35];
[lat1,lon1]=trajextract(region,lat,lon);

figure
cellplot(lon,lat,'D'),latratio,regionplot(region),hold on
cellplot(lon1,lat1),axis tight,boxon
title('Floats from Eastern Basin experiment')
xlabel('Longitude'),ylabel('Latitude')
