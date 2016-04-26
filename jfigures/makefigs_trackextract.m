%MAKEFIGS_TRACKEXTRACT  Makes a sample figure for TRACKEXTRACT.

load tpjaos
use tpjaos

region=[-66 -41 52 67];
[lat,lon,ssh,mss,atd]=trackextract(lat,lon,ssh,mss,atd,region);

figure
axis(region)
topoplot continents
plot(lon,lat,'linewidth',2)
latratio(58.5)

