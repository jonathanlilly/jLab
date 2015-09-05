function[]=makefigs_instmom
%MAKEFIGS_INSTMOM  Makes some sample figures for INSTMOM.

%First figure
dt=0.5;
t=(-200:dt:200)';

x=morsewave(length(t),1,24,2,2*pi/20.*dt);
x=x./maxmax(abs(x));
[a,om,rho1,rho2,rho3,rho4]=instmom(dt,x);
eta=om-sqrt(-1)*rho1;


figure,
subplot(231),plot(t,abs(x)),hold on,uvplot(t,x); ylim([-1.1 1.1]),linestyle 2k k k-- k:
title('Analytic Signal \psi_{24,2}(t)'),ytick(-.8:.4:.8),fixlabels([0 -1]), xlim([-100 100])
subplot(232),uvplot(t,eta);ylim([-.45 .45]),linestyle k k-- k:
title('Complex Instantaneous Frequency \eta(t)'),xlim([-100 100])
subplot(233),plot(t,real(rho1)./om),ylim([-.25 .25]),linestyle k k:
title('Modulation Function #1 \rho_1(t)'),xlim([-100 100])
subplot(234),plot(t,abs(rho2)./om.^2),hold on,uvplot(t,rho2./om.^2),ylim([-.1 .1]),linestyle 2k k k-- k:
title('Modulation Function #2 \rho_2(t)'),ytick(-.15:.05:.15),xlim([-100 100])
subplot(235),plot(t,abs(rho3)./om.^3),hold on,uvplot(t,rho3./om.^3),ylim([-.1 .1]),linestyle 2k k k-- k:
title('Modulation Function #3 \rho_3(t)'),ytick(-.15:.05:.15),xlim([-100 100])

subplot(236),plot(t,abs(rho4)./om.^4),hold on,uvplot(t,rho4./om.^4),ylim([-.1 .1]),hlines(0),linestyle 2k k k-- k:
title('Modulation Function #4 \rho_4(t)'),ytick(-.15:.05:.15),xlim([-100 100])


for i=1:6
    subplot(2,3,i),vlines([-22 22],'k:')
    hlines(0,'k:'),
    xlim([-75 75]),xtick(-60:20:60),fixlabels([0 -2]), 
end

letterlabels(4)


%Second figure
dt=0.5;
t=(-100:dt:100)';
x=10*morsexpand(2000,t,24,2,2*pi/20);
x=x./maxmax(abs(x));

[a,om,rhon{1},rhon{2},rhon{3},rhon{4},rhon{5},rhon{6}]=instmom(dt,x);
eta=om-sqrt(-1)*rhon{1};

tmid=round(length(t)/2);

xo=x(tmid);
omo=om(tmid);


xhat=zeros(length(x),7);
xhat(:,1)=xo.*rot(t.*omo);
for n=1:6
    xhat(:,n+1)=xhat(:,1).*frac(1,factorial(n)).*(t.^n).*rhon{n}(tmid);
end

xhat=cumsum(xhat,2);

figure
subplot(121),plot(t,real(x)),hold on,plot(t,real(xhat(:,3:2:end))),xlim([-35 35]),hlines(0)
ylim([-1.1 1.1]),ytick(-.8:.4:.8),xtick(-60:20:60),fixlabels([0 -1]),linestyle 2k k k-- k-. k: 
title('Demodulate Expansion of \Re\{\psi_{24,2}(t)\} about t=0'),vlines([-22 22],'k:')
vlines(0,'k:')

subplot(122),plot(t,imag(x)),hold on,plot(t,imag(xhat(:,3:2:end))),xlim([-35 35]),hlines(0)
ylim([-1.1 1.1]),ytick(-.8:.4:.8),xtick(-60:20:60),fixlabels([0 -1]),linestyle 2k k k-- k-. k:
title('Demodulate Expansion of \Im\{\psi_{24,2}(t)\} about t=0'),vlines([-22 22],'k:')
vlines(0,'k:')

letterlabels(4),packfig(1,2,'columns')


