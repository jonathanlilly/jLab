function[]=makefigs_regionplot
%MAKEFIGS_REGIONPLOT  Makes some sample figures for REGIONPLOT.

load ebasnfloats
use ebasnfloats
region=[-30 -21 24 35];
figure,cellplot(lon,lat),latratio,regionplot(region),axis tight
title('Floats from Eastern Basin experiment')
xlabel('Longitude'),ylabel('Latitude'),boxon
ax=axis;

if exist('m_map')==7
    figure,
    m_proj('albers equal-area conic','lon',ax(1:2),'lat',ax(3:4));
    cellplot(lon,lat,'m_map')
    m_grid('linestyle','none')
    regionplot(region,'m_map')
    title('Floats from Eastern Basin experiment using M-MAP')
else
    disp('REGIONPLOT skipping M_MAP example figure since M_MAP not detected.')
end