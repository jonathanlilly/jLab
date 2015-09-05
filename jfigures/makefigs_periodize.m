function[]=makefigs_periodize
%MAKEFIGS_PERIODIZE  Makes a sample figure for PERIODIZE.

load qgsnapshot
[xp,yp,fp]=periodize(100,200,qgsnapshot.x,qgsnapshot.y,qgsnapshot.psi);

figure
pcolor(xp,yp,fp),set(gca,'dataaspectratio',[1 1 1]),axis tight, shading interp
vlines(xp([100 end-100]),'w')
hlines(yp([200 end-200]),'w')
xtick([-5:1:5]*1000),ytick([-5:1:5]*1000)
title('Periodization of 1024x1024 QG turbulence with N=100 and M=200.')

