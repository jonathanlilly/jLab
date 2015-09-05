function[]=makefigs_closedcurves
%MAKEFIGS_CLOSEDCURVES  Makes a sample figure for CLOSEDCURVES.
 
load qgsnapshot, use qgsnapshot
dx=qgsnapshot.x(2)-qgsnapshot.x(1);
[cv,zeta,N,S,P]=psi2fields(dx,qgsnapshot.psi);
P=frac(P,std(P(:)));

[xc,yc]=closedcurves(x,y,P,-2);
[xp,yp]=closedcurves(x,y,P,-2,'periodic',100,100);
[xpax,ypax,fp]=periodize(100,100,x,y,P);

figure,jpcolor(xpax,ypax,fp),axis equal,axis tight
colormap gray, flipmap, caxis([-5 5]),
vlines(xpax([100 end-100]),'w'),hlines(ypax([100 end-100]),'w')
xtick([-5:1:5]*1000),ytick([-5:1:5]*1000),noxlabels,noylabels
cellplot(xp,yp,'2b'),hold on,cellplot(xc,yc,'r'),fontsize 11
title('Periodized (blue) non-periodized (red) curves for a QG model')
h=colorbar;
h.Label.String='Okubo-Weiss Parameter normalized by its own standard deviation';

%To print
if 0
    currentdir=pwd;
    cd([whichdir('jlab_license') '/figures'])
    print -dpng closedcurves
    crop closedcurves.png
    cd(currentdir)
end