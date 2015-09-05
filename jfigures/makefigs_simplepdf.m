function[]=makefigs_simplepdf
%MAKEFIGS_SIMPLEPDF  Makes a sample figure for SIMPLEPDF.

x=(-100:.1:100)';
mu=25;
sig=10;
f=simplepdf(x,mu,sig,'gaussian');
figure,plot(x,f),vlines(mu,'r')
title('Gaussian with mean 25 and standard deviation 10')
