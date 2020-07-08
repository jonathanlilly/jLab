function[]=makefigs_maternspec
%MAKEFIGS_MATERNSPEC  Makes some sample figures for MATERNSPEC.

%First figure
N=1000;
alpha=[1.1 1.5 2 4 8]';
h=1;

figure
subplot(1,3,1)
[f,spp,snn]=maternspec(1,N,1,alpha,h);
plot([-flipud(f);f]/pi,[flipud(snn);spp])
xlim([-1 1])
title('Matern spectrum')
xlabel('Frequency (cycles/point)')

subplot(1,3,2)
[to,R]=materncov(0.01,10000,1,alpha,h,'full');
plot(to,R)

title('Matern autocovariance')
xlabel('Time'),xlim([-7 7 ])

subplot(1,3,3)
to=[-2:.001:14]';
g=maternimp(to,alpha,h);
plot(to,g)
title('Matern Green''s function')
xlabel('Time'),xlim([-1.5 12.5 ])


for i=1:3
    subplot(1,3,i)
    ytick off
    %linestyle 2K G k-- K 2G
    vlines(0,'k:')
end
packfig(1,3,'columns')
letterlabels(2)
legend('$\alpha=1.1$','  $\alpha=1.5$','  $\alpha=2$','  $\alpha=4$','  $\alpha=8$','interpreter','latex')
fontsize 18 14 14 14
set(gcf,'paperposition',[1 1 6 3])


%Second figure
N=1000;
alpha=1;
h=1./[3 8]';
nu=1;

figure
subplot(1,3,1)
[f,spp,snn]=maternspec(1,N,1,alpha,h,nu);
plot([-flipud(f);f]/pi,[flipud(snn);spp])
xlim([-1 1])
title('Complex OU spectrum')
xlabel('Frequency (cycles/point)')
linestyle k 2G
vlines(1/pi,'D')

subplot(1,3,2)
[to,R]=materncov(0.01,10000,1,alpha,h,nu,'full');
uvplot(to,R+6*vrep([0 1]*(1+sqrt(-1)),size(R,1),1))
title('Complex OU autocovariance')
xlabel('Time'),xlim([-3.75 3.75 ]*10),ylim([-1.75 10.25])
linestyle k 2G k-- G--

subplot(1,3,3)
to=10*[-2:.01:12]';
g=maternimp(to,alpha,h-sqrt(-1)*nu);
uvplot(to,4*g+6*vrep([0 1]*(1+sqrt(-1)),size(g,1),1))
title('Complex OU Green''s function')
xlabel('Time'),xlim([-1.75 5.75 ]*10),ylim([-1.75 10.25])
linestyle k 2G k-- G--

for i=1:3
    subplot(1,3,i)
    ytick off
    vlines(0,'k:')
end

packfig(1,3,'columns')
letterlabels(2)
legend(' $\omega_o/h=3$',' $\omega_o/h=8$','interpreter','latex')
fontsize 18 14 14 14
set(gcf,'paperposition',[1 1 6 3])

