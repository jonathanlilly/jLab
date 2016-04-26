function[]=makefigs_widgist
%MAKEFIGS_WIDGIST  Makes some sample figures for WIGDIST.

%First figure
N=1000;
fs=0.01;
ga=3;
be=1;
w1=morlwave(N,2*pi*1.502/2/pi,2*pi*fs,'energy');
w2=morsewave(N,1,ga,be,2*pi*fs,'energy');
om1=vdiff(unwrap(imag(log(w1))),1);
om2=vdiff(unwrap(imag(log(w2))),1);

[d1,f1]=wigdist(w1,-15,25);
[d2,f2]=wigdist(w2,-15,25);
d1=d1./maxmax(d1);
d2=d2./maxmax(d2);

t=1:size(d1,1);t=t-mean(t);

ci=logspace(-2,0,10);

figure
set(gcf,'defaulttextinterpreter','latex')
subplot(221)
plot(t,[real(w1) imag(w1) abs(w1)])
linestyle k- k-- 2k
hlines(0,'k:')
ylabel('Amplitude')
xlabel('Time')
title('Morlet Wavelet with $\omega_\nu=1.5$, $P_{\nu}^2=3$')
xlim([-80 80]),ylim([-.15 .15])
xtick((-60:20:60))
fixlabels([0 -2])

subplot(222)
plot(t,[real(w2) imag(w2) abs(w2)])
linestyle k- k-- 2k
hlines(0,'k:')
ylabel('Amplitude')
xlabel('Time')
title('Morse Wavelet with $\gamma=3$ and $\beta=1$, $P_{\beta,\gamma}^2=3$')
xlim([-80 80]),ylim([-.15 .15])
xtick((-60:20:60))


subplot(223)
contourf(t,1000*f1,abs(d1'),ci),hold on
hlines(0,'k:')
ylabel('Cyclic Frequency x $10^3$')
xlabel('Time')
xlim([-80 80])
xtick((-60:20:60))
plot(t,1000*om1/2/pi,'w','linewidth',3)
plot(t,1000*om1/2/pi,'k--','linewidth',1.5)
colormap gray,flipmap,

subplot(224)
contourf(t,1000*f1,abs(d2'),ci),hold on
hlines(0,'k:')
xlabel('Time')
xlim([-80 80])
xtick((-60:20:60))
plot(t,1000*om2/2/pi,'w','linewidth',3)
plot(t,1000*om2/2/pi,'k--','linewidth',1.5)
colormap gray,flipmap,

letterlabels(1)
packfig(2,2)


%Second figure
N=2000;
fs=0.01;
ga=3;
be=10;
w1=morlwave(N,2*pi*5.47./2./pi,2*pi*fs,'energy');
w2=morsewave(N,1,ga,be,2*pi*fs,'energy');
om1=vdiff(unwrap(imag(log(w1))),1);
om2=vdiff(unwrap(imag(log(w2))),1);

[d1,f1]=wigdist(w1,10,30);
[d2,f2]=wigdist(w2,10,30);
d1=d1./maxmax(d1);
d2=d2./maxmax(d2);

t=1:size(d1,1);t=t-mean(t);

ci=logspace(-2,0,10);

figure
set(gcf,'defaulttextinterpreter','latex')
subplot(221)
plot(t,[real(w1) imag(w1) abs(w1)])
linestyle k- k-- 2k
hlines(0,'k:')
ylabel('Amplitude')
xlabel('Time')
title('Morlet Wavelet with $\omega_\nu=5.5$, $P_{\nu}^2=30$')
xlim([-80 80]*2.75),ylim([-.1 .1]),xtick((-60:20:60)*3)
fixlabels([0 -2])

subplot(222)
plot(t,[real(w2) imag(w2) abs(w2)])
linestyle k- k-- 2k
hlines(0,'k:')
ylabel('Amplitude')
xlabel('Time')
title('Morse Wavelet with $\gamma=3$ and $\beta=10$, $P_{\beta,\gamma}^2=30$')
xlim([-80 80]*2.75),ylim([-.1 .1]),xtick((-60:20:60)*3)

subplot(223)
contourf(t,1000*f1,abs(d1'),ci),hold on
hlines(0,'k:')
ylabel('Cyclic Frequency x $10^3$')
xlabel('Time')
xlim([-80 80]*2.75),xtick((-60:20:60)*3),ytick((6:2:14))
plot(t,1000*om1/2/pi,'w','linewidth',3)
plot(t,1000*om1/2/pi,'k--','linewidth',1.5)
colormap gray,flipmap,

subplot(224)
contourf(t,1000*f1,abs(d2'),ci),hold on
hlines(0,'k:')
xlabel('Time')
xlim([-80 80]*2.75),xtick((-60:20:60)*3),ytick((6:2:14))
plot(t,1000*om2/2/pi,'w','linewidth',3)
plot(t,1000*om2/2/pi,'k--','linewidth',1.5)
colormap gray,flipmap,

letterlabels(1)
packfig(2,2)


%Third figure
N=1000;
fs=0.01;
w1=morsewave(N,1,2,6,2*pi*fs,'energy');
w2=morsewave(N,1,3,2*6./3,2*pi*fs,'energy');
w3=morsewave(N,1,4,3,2*pi*fs,'energy');

[d1,f1]=wigdist(w1,3,20);
[d2,f2]=wigdist(w2,3,20);
[d3,f3]=wigdist(w3,3,20);
d1=d1./maxmax(d1);
d2=d2./maxmax(d2);
d3=d3./maxmax(d3);
om1=vdiff(unwrap(imag(log(w1))),1);
om2=vdiff(unwrap(imag(log(w2))),1);
om3=vdiff(unwrap(imag(log(w3))),1);

t=1:size(d1,1);t=t-mean(t);
ci=logspace(-2,0,10);

figure
set(gcf,'defaulttextinterpreter','latex')
subplot(231)
plot(t,[real(w1) imag(w1) abs(w1)])
linestyle k- k-- 2k
hlines(0,'k:')
ylabel('Amplitude')
title('Morse with $\gamma=2$ and $\beta=6$')
ylim([-.13 .13]),xlim([-190 190]),fixlabels([0 -2])
xtick((-200:40:200))

subplot(232)
plot(t,[real(w2) imag(w2) abs(w2)])
linestyle k- k-- 2k
hlines(0,'k:')
ylabel('Amplitude')
title('Morse with $\gamma=3$ and $\beta=4$')
ylim([-.13 .13]),xlim([-190 190])
xtick((-150:50:150))

subplot(233)
plot(t,[real(w3) imag(w3) abs(w3)])
linestyle k- k-- 2k
hlines(0,'k:')
ylabel('Amplitude')
title('Morse with $\gamma=4$ and $\beta=3$')
ylim([-.13 .13]),xlim([-190 190])
xtick((-150:50:150))

subplot(234)
contourf(t,1000*f1,abs(d1'),ci),hold on
hlines(0,'k:')
ylabel('Cyclic Frequency x $10^3$')
ylim([0 19])
xlim([-190 190])
xtick((-150:50:150))
plot(t,1000*om1/2/pi,'w','linewidth',3)
plot(t,1000*om1/2/pi,'k--','linewidth',1.5)
colormap gray,flipmap,
xlabel('Time')

subplot(235)
contourf(t,1000*f2,abs(d2'),ci),hold on
hlines(0,'k:')
ylim([0 19])
xlim([-190 190])
xtick((-150:50:150))
plot(t,1000*om2/2/pi,'w','linewidth',3)
plot(t,1000*om2/2/pi,'k--','linewidth',1.5)
xlabel('Time')

subplot(236)
contourf(t,1000*f3,abs(d3'),ci),hold on
hlines(0,'k:')
ylim([0 19])
xlim([-190 190])
xtick((-150:50:150))
plot(t,1000*om3/2/pi,'w','linewidth',3)
plot(t,1000*om3/2/pi,'k--','linewidth',1.5)
xlabel('Time')

letterlabels(1)
packfig(2,3)
