function[]=makefigs_morlfreq
%MAKEFIGS_MORLFREQ  Makes a sample figure for MORLFREQ.

fm=linspace(0,4,1000);
[fn,fmin]=morlfreq(fm);

figure
plot([fn' fn'],[fm' fmin']),axis equal,
dlines(1,'k--')
axis([0 4 -1 4]),%vlines(1,'k:'),
hlines(0,'k:')

xlabel('Carrier Frequency (Radian)')
ylabel('Frequency of Max and Min (Radian)')
title('Frequencies of Morlet Maximum and Minimum')
