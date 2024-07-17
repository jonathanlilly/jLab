function[]=makefigs_jhermfun
%MAKEFIGS_JHERMFUN  Makes a sample figure for JHERMFUN.

t=(-5:0.1:5.1);
tnorm=(t-t(1))./(t(length(t))-t(1));
h=jhermfun(t,5);
figure,plot(t,h),
title('Hermite functions H_0 -- H_5')
axis tight

