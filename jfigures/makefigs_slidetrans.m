function[]=makefigs_slidetrans
%MAKEFIGS_SLIDETRANS  Makes a sample figure for SLIDETRANS.

figure
M=3000;
t=(0:M-1)';
N=500;
w=hermfun((-N:N)'./(N/4),0);
w=w./sqrt(w'*w);
fs=2*pi*(1:30)./1000;
clear x
x(1:M/2,1)=sin(2*pi*t(1:M/2)./70/3);
x(M/2:M,1)=sin(2*pi*t(M/2:M)./70);
y=slidetrans(x,w,fs,'zeros');
h=wavespecplot(t,x,1./fs,abs(y),1/2);
hlines(70*3/2/pi),hlines(70/2/pi)