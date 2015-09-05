function[]=makefigs_curvemoments
%MAKEFIGS_CURVEMOMENTS  Makes a sample figure for CURVEMOMENTS.
 
load qgsnapshot

[cv,zeta,N,S,P]=psi2fields(qgsnapshot.psi);
P=frac(P,std(P(:)));

[xc,yc]=closedcurves(qgsnapshot.x,qgsnapshot.y,P,-2);
[xo,yo,L,R,D,a,b,theta]=curvemoments(xc,yc);

figure,jpcolor(qgsnapshot.x,qgsnapshot.y,P),axis equal, axis tight,
hold on,colormap gray,flipmap,cellplot(xc,yc,'2b'),
[k,l]=ab2kl(a,b);ellipseplot(k,l,theta,xo+sqrt(-1)*yo,'r')
title('Curves of Okubo-Weiss (blue), moment ellipses (red), and centroids (black)')
caxis([-5 5]),xtick([-5:1:5]*1000),ytick([-5:1:5]*1000),plot(xo,yo,'k.')
noxlabels,noylabels,fontsize 11
h=colorbar;
h.Label.String='Okubo-Weiss Parameter normalized by its own standard deviation';

%To print
if 0
    currentdir=pwd;
    cd([whichdir('jlab_license') '/figures'])
    print -dpng curvemoments
    crop curvemoments.png
    cd(currentdir)
end