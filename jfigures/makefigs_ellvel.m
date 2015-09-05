function[]=makefigs_ellvel
%MAKEFIGS_ELLVEL  Makes a sample figure for ELLVEL.

lambda=(0:.001:1)';
kappa=1+0*lambda;
phi=(1:length(lambda))'/10;
theta=0*lambda;

fact=1e5;dt=1/24;
vm=ellvel(dt,kappa,lambda,theta,phi,fact,'geometric');
vgamma=ellvel(dt,kappa,lambda,theta,phi,fact,'circulation');
veke=ellvel(dt,kappa,lambda,theta,phi,fact,'kineticenergy');
vbar=ellvel(dt,kappa,lambda,theta,phi,fact,'average');

figure
plot(lambda,[vgamma,vm,vbar,veke]./maxmax(vm));
linestyle 2k k-- k 2G
e=[.2486 .967];
axis([0.01 1 0 1.5]),axis square 
vlines(e.^2./(2-e.^2),'k:')
legend('V_\Gamma','V_M','V_{Bar}','V_{EKE}')
title('Velocity measures for constant \kappa and frequency')
xlabel('Ellipse linearity \lambda')
ylabel('Mean velocity measures')
set(gcf,'paperposition', [2 2 3.5 3.5])
xtick(.1),ytick(.1),fixlabels(-1)
