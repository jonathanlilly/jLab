function[]=makefigs_patchcontourf
%MAKEFIGS_PATCHCONTOURF  Makes a sample figure for PATCHCONTOURF.

load jtopo
use jtopo

figure
h=patchcontourf(lon,lat,topo,0,'k');latratio(30)
axis([-180 180 -75 75]),boxon,axis tight
title('Continents filled with PATCHCONTOURF')

%Make a second figure if m_map is installed
if exist('m_map')==7
    figure
    m_proj('miller','lat',75);
    m_grid('linestyle','none')
    patchcontourf(lon,lat,topo,0,'k','m_map');
    title('Continents filled with PATCHCONTOURF, M-MAP option')
end