function[]=makefigs_doublen
%MAKEFIGS_DOUBLEN   Makes a sample figure for DOUBLEN.

%Illustration of length doubling
M=100;
N=1;
x=randn(M,N);
y=doublen(x);

t1=(0:100-1)';
t2=(0:1/2:100-1/2)';

figure
subplot(121)
plot(t2,y,'r'),hold on, plot(t1,x,'b')
plot(t2(1:2:end),y(1:2:end,1),'ro')
title('Red circles should fall along blue curve')
xlabel('Time')
subplot(122)
plot(t1,abs(fft(x(:,1)))),hold on,plot(t2,abs(fft(y(:,1))))
linestyle b r
xlabel('Frequency')


