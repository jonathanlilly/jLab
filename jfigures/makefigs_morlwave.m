function[]=makefigs_morlwave
%MAKEFIGS_MORLWAVE  Makes a sample figure for MORLWAVE.

N=1000; % nombre de points
fs=2*pi.*flipud(logspace(log10(6./N),log10(60./N),10)');   %change to log space
[w,W]=morlwave(N,2*pi./4,fs,'energy');
t=1:length(w);
t=t-mean(t);

f=(0:N-1)'./N;
index=find(f>1/2);
f(index)=f(index)-1;

figure,
subplot(221)
plot(t,abs(w))
title('Modulus of Morlet wavelets in time')
 xlim([-100 100])

subplot(222)
plot(f,abs(W)),vlines(fs/2/pi)
title('Modulus of Morlet wavelets in frequency')
 xlim([-.15 .15])
 
subplot(223)
for i=1:size(w,2)
  t1=t.*fs(i)/2/pi;
  plot(t1,real(w(:,i))./max(abs(w(:,i))),'b.'),hold on
  plot(t1,imag(w(:,i))./max(abs(w(:,i))),'g.'),hold on
end
axis([-1 1 -1 1])
title('Stretched and rescaled Morlet wavelets')

i=5;
subplot(2,2,4)
uvplot(t,w(:,i)),hold on,
plot(t,abs(w(:,i)),'k');
dw=diff(abs(w(:,i)));
vlines(t(dw==max(dw)|dw==min(dw)),'k--')
title('Morlet wavelet with f_{max}=1')
axis([-40 40 -.2 .2])


