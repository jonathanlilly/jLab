function[]=makefigs_morsewave
%MAKEFIGS_MORSEWAVE  Makes some sample figures for MORSEWAVE.


%First figure
figure

N=256*4;
be=5;
ga=2;
K=3;
fs=2*pi/8/4;

[x,X]=morsewave(N,K,ga,be,fs,'energy');
f=(0:1:N-1)'./N;

t=(1:length(x))'-length(x)/2;
ax=[-60 60 -maxmax(abs(x))*1.05 maxmax(abs(x))*1.05];
subplot 321
  uvplot(t,x(:,1));axis(ax)
  title('Morse wavelets, time domain')
subplot 323
  uvplot(t,x(:,2));axis(ax)
subplot 325
  uvplot(t,x(:,3));axis(ax)

ax=[0 120./N -maxmax(abs(X))*1.05 maxmax(abs(X))*1.05];
subplot 322
  plot(f,abs(X(:,1))),axis(ax),vlines(fs);
title('Morse wavelets, frequency domain')
subplot 324
  plot(f,abs(X(:,2))),axis(ax),vlines(fs);
subplot 326
plot(f,abs(X(:,3))),axis(ax),vlines(fs);
  
%Second figure
figure

N=2^16;
fs=1/512;
t=1:N;t=t-mean(t);t=t.*fs;
clear psi psif 
P=zeros(4,5);
ga=[1/3 1 3 9 27]';
be=flipud(ga(1:end));
for i=1:5
    for j=1:5
        P(i,j)=sqrt(be(j)*ga(i));
        [psi(:,i,j) psif(:,i,j)]=morsewave(N,1,ga(i),be(j),2*pi*fs,'bandpass');  
    end
end
for i=1:5
    for j=1:5
        subplot(5,5,i+(j-1)*5)
        x=[squeeze(real(psi(:,i,j))) squeeze(imag(psi(:,i,j))) squeeze(abs(psi(:,i,j)))]./max(abs(psi(:,i,j)));
        plot(t./(P(i,j)./pi),x),xlim([-3 3]),ylim([-0.85 1.1])
        %linestyle b r k
        hlines(0,'k:'),%vlines([-1/2 1/2],'k:')
        %if iseven(i+(5-j+1)*5),set(gca,'color',[1 1 1]*0.9),end
        if j==5,xlabel(['$\gamma=' num2str(ga(i),3) '$'],'interpreter','latex'),end
        if i==1,ylabel(['$\beta=' num2str(be(j),3) '$'],'interpreter','latex'),end
        set(gca,'xtick',[]),set(gca,'ytick',[]),%set(gca,'xcolor','w'),set(gca,'ycolor','w')
        set(get(gca,'ylabel'),'color','k')
    end
end
subplot(5,5,3),title('The Generalized Morse Wavelets $\psi_{\beta,\gamma}(t)$','interpreter','latex')
packfig(5,5)

fontsize 14 10 10 10 
orient tall
set(gcf,'paperposition',[1 1 8 6])

if 0
    currentdir=pwd;
    cd([whichdir('jlab_license') '/figures'])
    print -dpng morsewave
    crop morsewave.png
    cd(currentdir)
end
