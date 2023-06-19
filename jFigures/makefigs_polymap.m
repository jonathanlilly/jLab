function[]=makefigs_polymap
%MAKEFIGS_POLYMAP  Makes some sample figures for POLYMAP.

%Example of polymap using Matlab's "Peaks" function,
%comparing constant, linear, and quadratic fits.
  
%Use peaks for testing... random sampling
[x,y,z]=peaks;
rng(1);
index=randperm(length(z(:)));
index=index(1:600);
[xdata,ydata,zdata]=vindex(x(:),y(:),z(:),index,1);

xo=(-3:.125:3);
yo=(-3:.125:3);

H=0.5;

%[xg,yg]=meshgrid(xo,yo);bool=(xg.^2+yg.^2)>1;
zhat=zeros(length(yo),length(xo),3);
for i=1:3
    zhat(:,:,i)=polymap(xdata,ydata,[],[],zdata,xo,yo,{i-1,H,2,1});
%    zhat(:,:,i)=polymap(xdata,ydata,[],[],zdata,xo,yo,{i-1,H,2,1},'mask',bool);
end

%vswap(zhat,inf,nan);
figure
subplot(2,2,1),contourf(xo,yo,z,(-10:10)),title('Gaussian Topography')
subplot(2,2,2),contourf(xo,yo,zhat(:,:,1),(-10:10)),title('Constant Fit')
subplot(2,2,3),contourf(xo,yo,zhat(:,:,2),(-10:10)),title('Linear Fit')
subplot(2,2,4),contourf(xo,yo,zhat(:,:,3),(-10:10)),title('Quadratic Fit')

for i=1:4
    subplot(2,2,i),hold on,caxis([-7 7]),nocontours,plot(xdata,ydata,'k+')
    plot(H*rot(0:.05:2*pi+.1),'k'),axis square,axis equal,axis tight
    xtick(-3:1:3),ytick(-3:1:3)
    colorbar('vertical')
end

figure
subplot(2,2,1),contourf(xo,yo,z,(-10:10)),caxis([-7 7]),title('Gaussian Topography')
subplot(2,2,2),contourf(xo,yo,zhat(:,:,1)-z,(-10:10)/3),caxis([-3 3]),title('Constant Deviation')
subplot(2,2,3),contourf(xo,yo,zhat(:,:,2)-z,(-10:10)/3),caxis([-3 3]),title('Linear Deviation')
subplot(2,2,4),contourf(xo,yo,zhat(:,:,3)-z,(-10:10)/3),caxis([-3 3]),title('Quadratic Deviation')

for i=1:4
    subplot(2,2,i),hold on,caxis([-7 7]),nocontours,plot(xdata,ydata,'k+')
    plot(H*rot(0:.05:2*pi+.1),'k'),axis square,axis equal,axis tight
    xtick(-3:1:3),ytick(-3:1:3)
    colorbar('vertical')
    if i>1,caxis([-1 1]),end
end

figure,
subplot(1,3,1),plot(z,zhat(:,:,1),'.'),title('Constant Fit')
subplot(1,3,2),plot(z,zhat(:,:,2),'.'),title('Linear Fit')
subplot(1,3,3),plot(z,zhat(:,:,3),'.'),title('Quadratic Fit')

for i=1:3
   subplot(1,3,i), axis square,axis equal,axis([-7 7 -7 7])  
   xtick(-6:2:6),ytick(-6:2:6),hlines(0,'k:'),vlines(0,'k:')
   plot([-1-sqrt(-1);1+sqrt(-1)]*7,'k')
   xlabel('Original Surface'),ylabel('Estimated Surface')
end
