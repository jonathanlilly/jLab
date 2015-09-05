function[]=makefigs_ecconv
%MAKEFIGS_ECCONV  Makes a sample figure for ECCONV.

figure
nu=pi/4*(0:.01:1)';
ecc=ecconv(nu,'nu2ecc');
lin=ecconv(nu,'nu2lin');
zeta=sqrt(1-lin.^2);

plot(nu/pi,[zeta,lin,ecc]);
linestyle 3k k k-- 
ytick(0:.2:1)
fixlabels(-1)
xtick([0:1/16:1/4])
set(gca,'xticklabel',['0 pi   ';'1/16 pi';'1/8 pi ';'3/16 pi';'1/4 pi ']);
axis tight
set(gca,'dataaspectratio',[1/5 1 1])
legend('Circularity \zeta','Linearity \lambda','Eccentricity e')

hold on
ii=(6:10:length(lin))';
ellipseplot(1/25+0*ii,lin(ii),pi/2+0*ii,nu(ii)./pi+sqrt(-1)*0.08,[1/5 1],'k')
ellipseplot(1/25+0*ii,lin(ii),0*ii,nu(ii)./pi+sqrt(-1)*0.08,[1/5 1],'C')

title('Measures of Eccentricity')
xlabel('Eccentricity Angle \nu=arctan(b/a)')
