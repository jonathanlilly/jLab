function[]=makefigs_topoplot
%MAKEFIGS_TOPOPLOT  Makes a sample figure for TOPOPLOT.

figure,topoplot,latratio(30)
title('Continents (black) and Continental Shelves (gray), from TOPOPLOT')

%To print
if 0
    currentdir=pwd;
    cd([whichdir('jlab_license') '/figures'])
    print -dpng topoplot
    crop topoplot.png
    cd(currentdir)
end

load ebasnfloats
use ebasnfloats
region=[-40.2800  -17.2400   19.1820   41.1340];

figure,topoplot(region,[-6:1/4:0],-2,'2w'),hold on
cellplot(lon,lat),latratio,axis tight,
title('Eastern Basin floats, with bathymetry')
xlabel('Longitude'),ylabel('Latitude'),boxon
h=colorbar;
if verLessThan('matlab','8.4.0')
    axes(h),hlines(-2,'3w'),ylabel('Bathymetry (km)')
else
    h.Label.String='Bathymetry (km)';
end
colormap(gca,'gray');flipmap
fontsize 14 12 12 12

if exist('m_map')==7
    figure,
    m_proj('albers equal-area conic','lon',region(1:2),'lat',region(3:4));
    cellplot(lon,lat,'m_map')
    m_grid('linestyle','none'),hold on
    topoplot(region,[-6:1/4:0],-2,'m_map','2w'),
    title('Eastern Basin floats using M-MAP, with bathymetry')
    %h=colorbar;axes(h),hlines(-2,'3w'),ylabel('Bathymetry (km)')
    colormap(gca,'gray');flipmap
    fontsize 14 12 12 12
else
    disp('REGIONPLOT skipping M_MAP example figure since M_MAP not detected.')
end 

