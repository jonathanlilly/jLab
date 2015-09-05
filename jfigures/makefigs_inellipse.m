function[]=makefigs_inellipse
%MAKEFIGS_INELLIPSE  Makes a sample figure for INELLIPSE.

x=[-100:2:100];
y=[-100:2:100];
[xg,yg]=meshgrid(x,y);
z=xg+sqrt(-1)*yg;
z=conj(z(:))';

figure
bool=inellipse(z,30,1/2,pi/3,20+sqrt(-1)*30);
plot(z(bool),'r.'),hold on,plot(z(~bool),'b.')
ellipseplot(30,1/2,pi/3,20+sqrt(-1)*30,'5k');
title('Illustration of INELLIPSE')
