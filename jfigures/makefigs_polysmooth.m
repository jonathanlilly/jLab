function[]=makefigs_polysmooth
%MAKEFIGS_POLYSMOOTH  Makes some sample figures for POLYSMOOTH.

%Example of polysmooth using Matlab's "Peaks" function,
%comparing constant, linear, and quadratic fits.
  
%Use peaks for testing... random sampling
[x,y,z]=peaks;
index=randperm(length(z(:)));
index=index(1:600);
[xdata,ydata,zdata]=vindex(x(:),y(:),z(:),index,1);

xo=(-3:.125:3);
yo=(-3:.125:3);

B=1;

[ds,xs,ys,zs]=twodsort(xdata,ydata,zdata,xo,yo,B);    
z0=polysmooth(ds,xs,ys,zs,B,0,'epan');
z1=polysmooth(ds,xs,ys,zs,B,1,'epan');
z2=polysmooth(ds,xs,ys,zs,B,2,'epan');


figure
subplot(2,2,1),contourf(xo,yo,z,(-10:10)),title('Gaussian Topography')
subplot(2,2,2),contourf(xo,yo,z0,(-10:10)),title('Constant Fit')
subplot(2,2,3),contourf(xo,yo,z1,(-10:10)),title('Linear Fit')
subplot(2,2,4),contourf(xo,yo,z2,(-10:10)),title('Quadratic Fit')

for i=1:4
    subplot(2,2,i),hold on,caxis([-7 7]),nocontours,plot(xdata,ydata,'k+')
    plot(B*rot(0:.05:2*pi+.1),'k'),axis square,axis equal,axis tight
    xtick(-3:1:3),ytick(-3:1:3)
    colorbar('vertical')
end

figure
subplot(2,2,1),contourf(xo,yo,z,(-10:10)),caxis([-7 7]),title('Gaussian Topography')
subplot(2,2,2),contourf(xo,yo,z0-z,(-10:10)/3),caxis([-3 3]),title('Constant Deviation')
subplot(2,2,3),contourf(xo,yo,z1-z,(-10:10)/3),caxis([-3 3]),title('Linear Deviation')
subplot(2,2,4),contourf(xo,yo,z2-z,(-10:10)/3),caxis([-3 3]),title('Quadratic Deviation')

for i=1:4
    subplot(2,2,i),hold on,nocontours,plot(xdata,ydata,'k+')
    plot(B*rot(0:.05:2*pi+.1),'k'),axis square,axis equal,axis tight
    xtick(-3:1:3),ytick(-3:1:3)
    colorbar('vertical')
end


figure,
subplot(1,3,1),plot(z,z0,'.'),title('Constant Fit')
subplot(1,3,2),plot(z,z1,'.'),title('Linear Fit')
subplot(1,3,3),plot(z,z2,'.'),title('Quadratic Fit')

for i=1:3
   subplot(1,3,i), axis square,axis equal,axis([-7 7 -7 7])  
   xtick(-6:2:6),ytick(-6:2:6),hlines(0,'k:'),vlines(0,'k:')
   plot([-1-sqrt(-1);1+sqrt(-1)]*7,'k')
   xlabel('Original Surface'),ylabel('Estimated Surface')
end
