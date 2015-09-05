function[]=makefigs_ellband
%MAKEFIGS_ELLBAND  Makes a sample figure for ELLBAND.

t=(0:1:925)';
cxe=zeros(length(t),3);

kappa=3*exp(2*0.393*(t/1000-1));
lambda=0.4+0*t;
phi=(t/1000*5)*2*pi;
theta=pi/4+0*t;

cxe(:,1)=ellsig(kappa,lambda,theta,phi);
om=vdiff(phi,1);  %Since theta is constant

kappa=2.5+0*t;
lambda=zeros(size(t));
for i=2:length(t)
    lambda(i)=real(lambda(i-1)+2*sqrt(1-lambda(i-1).^2)*0.025.*om(i));
end
lambda(1)=nan;
lambda(lambda>1)=1;  

cxe(:,2)=ellsig(kappa,lambda,theta,phi);

[kappa,lambda]=ab2kl(3+zeros(size(t)),2+zeros(size(t)));

%theta=-phi/16.3;  Retrograde precession
theta=phi/14.45;
cxe(:,3)=ellsig(kappa,lambda,theta,phi);

titlestr{1}='Increasing Magnitude';
titlestr{2}='Increasing Eccentricity';
titlestr{3}='Precession';

figure
for i=1:3
    subplot(1,3,i)
    plot(cxe(:,i)),hold on
    plot(cxe(1:200,i)),linestyle 1k 2k
    plot(cxe(2,i),'ko','markersize',10), plot(cxe(end,i),'kx','markersize',10)
    title(titlestr{i})
    xlabel('Displacement East'),axis equal,axis square
    axis([-3.2 3.2 -3.2 3.2])
    ylabel('Displacement North')
    xtick(-3:3)   
    ytick(-3:3)
end
packfig(1,3,'columns')
letterlabels(1)