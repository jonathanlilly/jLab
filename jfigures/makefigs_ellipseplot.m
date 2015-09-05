function[]=makefigs_ellipseplot
%MAKEFIGS_ELLIPSEPLOT  Makes a sample figure for ELLIPSEPLOT.

a=ones(10,1);
b=(1:10)'./10;
[k,l]=ab2kl(a,b);
th=linspace(0,pi,10)';
x=linspace(0,1,10)'*20;

figure,
ellipseplot(k,l,th,x,'phase',0*th,'npoints',16,'g')
title('Counterclockwise rotating ellipse becoming circle')
set(gca,'dataaspectratio',[1 1 1])

h=ellipseplot(k,l,th,x,'npoints',16,'g','ellipses');
lato=5:10:85;
lono=0;

phi=[0:1:1000]'/10;
kappa=100+0*phi;
lambda=phi./maxmax(phi);
theta=(pi/2).*phi./maxmax(phi);
z=ellsig(kappa,lambda,theta,phi);

figure
for i=1:9
    subplot(3,3,i)
    [lat,lon]=xy2latlon(z,lato(i),lono);
    plot(lon,lat,'r'),latratio(lato(i)),hold on,
    ellipseplot(kappa,lambda,theta,lono+sqrt(-1)*lato(i)+0*kappa,latratio(lato(i)),'skip',100,'scale',1);
    axis tight
    if i==2
        title('Checking plotted ellipse size and shape at various latitudes')
    end
    xlabel('Longitude'),ylabel('Latitude')
end


