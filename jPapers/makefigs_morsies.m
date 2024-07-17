function[varargout]=makefigs_morsies(str)
%MAKEFIGS_MORSIES  Make figures for Lilly and Olhede (2009a).
% 
%   This function makes all figures for the paper 
%
%   Lilly, J. M., and S. C. Olhede (2009). Higher-order properties of 
%      analytic wavelets. IEEE Transactions on Signal Processing, 57 (1),
%      146--160.
%
%   Usage: makefigs_morsies
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2023 J.M. Lilly --- type 'help jlab_license' for details

if nargin==0
    str='noprint';
end

%/************************************************************
%Linear chirp
N=2000;
t=(0:N-1)';t=t-mean(t);
x=cos(t.^2/5000).*exp(-(t./250).^2);
%plot(t,x)
 
om=vdiff(t.^2/5000,1);

fs=(1./(logspace(log10(10),log10(2000),50)'));
%fs=[-flipud(fs);fs];

psi1=morsewave(N,1,3,1,2*pi*fs,'energy');
psi2=morlwave(N,2*pi*1.5./2./pi,2*pi*fs,'energy');
wx1=wavetrans(x,psi1,'zeros');
wx2=wavetrans(x,psi2,'zeros');
wx1=wx1./maxmax(abs(wx1));
wx2=wx2./maxmax(abs(wx2));
figure,
subplot(3,1,1),plot(t,[x exp(-(t./250).^2)]),linestyle k 2k,xlim([-400 400]),hlines(0,'k:'),ylim([-1.2 1.2]),ytick(-1:.5:1),fixlabels([0 -1])
subplot(3,1,2),contourf(t,fs*100,(abs(wx1)'),(0:.1:1)),xlim([-400 400]),ylim([0 .1]/2*100), colormap gray,flipmap,nocontours
subplot(3,1,3),contourf(t,fs*100,(abs(wx2)'),(0:.1:1)),xlim([-400 400]),ylim([0 .1]/2*100), colormap gray,flipmap,nocontours

for i=2:3
    subplot(3,1,i)
    if i>1
        hold on,
        plot(t,100*om./2./pi,'w','linewidth',1)
        plot(t,100*om./2./pi,'k--','linewidth',1),
        plot(t,-100*om./2./pi,'w','linewidth',1)
        plot(t,-100*om./2./pi,'k--','linewidth',1),
        ytick(0:1:3),ylim([0 3.2])
        ylabel('Frequency $\omega/(2 \pi) \times 100$')
    end
    xtick(-400:200:400)
end
subplot(3,1,1),title('A Gaussian-Enveloped Chirp')

packfig(3,1,'rows')

fontsize 8 8 8 8
letterlabels(2)

orient tall
set(gcf,'paperposition',[1 1 2.75 6])
if strcmpi(str,'print')
    %
    print -deps morsies_linear_chirp.eps
end
%\************************************************************

%/************************************************************
figure

N=512*2*4;
fs=1/8/4/2;
t=1:N;t=t-mean(t);t=t.*fs;
clear psi psif 
P=zeros(4,4);
for i=1:4
    for j=0:3
        P(i,j+1)=sqrt(i*j);
        if j==0
           %This is due to the definition of the beta=0 case in morsefreq
           [psi(:,i,j+1) psif(:,i,j+1)]=morsewave(N,1,i,j,2*pi*fs*frac(morsefreq(i,0),morsefreq(i,1)), 'bandpass');
        else 
           [psi(:,i,j+1) psif(:,i,j+1)]=morsewave(N,1,i,j,2*pi*fs,'bandpass');
        end
        
    end
end
%P=P(:);
%P(1:4)=P(5:8);
P(:,1)=P(:,2);

for i=1:4
    for j=1:4
        subplot(4,5,i+(j-1)*5)
        x=[squeeze(real(psi(:,i,j))) squeeze(imag(psi(:,i,j))) squeeze(abs(psi(:,i,j)))]./max(abs(psi(:,i,j)));
        plot(t./(P(i,j)./pi),x),xlim([-2 2]),box on,ylim([-0.85 1.1])
        linestyle k k-- 1.5k
        hlines(0,'k:')
        set(gca,'xtick',[]), set(gca,'ytick',[])
        set(gca,'xcolor','w')
        set(gca,'ycolor','w')
        %xlim([-1.5 1.5]/2)
    end
end


for j=1:4
    subplot(4,5,5+(j-1)*5)
    plot(squeeze(psif(:,:,j)))
    xtick([0 65 210]),xlim([0 200]),ytick([-1 3]),%vlines(65,'k:')
    set(gca,'xticklabel',[]),set(gca,'yticklabel',[]),linestyle k k-- k-. k:
    boxoff
    set(gca,'ticklen',[0.025 0.025]*2)
end


subplot(4,5,1),title('Cauchy Family ($\gamma=1$)')
subplot(4,5,2),title('Gaussian Family ($\gamma=2$)')
subplot(4,5,3),title('Airy Family ($\gamma=3$)')
subplot(4,5,4),title('Hyper-Gaussian Family ($\gamma=4$)')
subplot(4,5,5),title('Frequency Domain')
subplot(4,5,1),text(-1.8,.9,'$\beta=0$')
subplot(4,5,6),text(-1.8,.9,'$\beta=1$')
subplot(4,5,11),text(-1.8,.9,'$\beta=2$')
subplot(4,5,16),text(-1.8,.9,'$\beta=3$')


packfig(4,5)
legend(gca,'1','2','3','4')
legend(gca,'Mod','Re','Im')

fontsize 10 10 10 10
orient landscape
set(gcf,'paperposition',[1 1 10 6])

if strcmpi(str,'print')
    print -deps morsie_families.eps
end
%\*********************************************************


%/********************************************************************
%The code for this figure is in MAKEFIGS_WIGDIST
h=gcf;
wigdist --f

figure(h.Number+1)
set(gcf,'paperposition',[0.25 1 7 5])
fontsize 12 10 10 10
if strcmpi(str,'print')
    print -deps morsie_morlet_wigdist.eps
end

figure(h.Number+2)
set(gcf,'paperposition',[0.25 1 7 5])
fontsize 12 10 10 10
if strcmpi(str,'print')
    print -deps morsie_morlet_wigdist_long.eps
end

figure(h.Number+3)
orient landscape
set(gcf,'paperposition',[1 1 10 6])
fontsize 12 10 8 10
if strcmpi(str,'print')
    print -deps morsie_wigdist_three.eps
end
%/*****************************************

%/*****************************************
p1=[ (1:.01:2) (2.1:.1:50) (55:5:100)];
alpha=(-1.5:.025:1.5)*2.3./1.5;

[p,alpha]=meshgrid(p1,alpha);
ga=p.*alpha+3;
be=squared(p)./ga;
index=find(ga<=0);
ga(index)=nan;
be(index)=nan;

%index=find(be>frac(ga-1,2));

[a,dt,dom]=morsebox(ga,be);
index=find(ga<1);
[fm,fe,fi,cf] = morsefreq(ga,be);

figure
subplot(2,3,1),h=contour(p./pi,alpha,ga,(3:10).^2,'k');hold on
h=contour(p./pi,alpha,ga,[2 3 4],'k','linewidth',2); 
h=contour(p./pi,alpha,ga,[2 3 4],'w','linewidth',1);xlog, 
[p1,a1]=morseprops(1,[(0:.1:100) (101:1:1000) (1000:10:10000)] );
plot(p1./pi,a1,'k','linewidth',2);
ylabel('Demodulate Skewness $\Im\{\alpha_{3;\beta,\gamma}\}$')
text(0.55,-2,'(a)  $\gamma$ parameter')

subplot(2,3,2),contour(p./pi,alpha,be,(1:10).^2,'k');xlog,hold on
h=contour(p./pi,alpha,be,[2 3 4],'k','linewidth',2);xlog,hold on
h=contour(p./pi,alpha,be,[2 3 4],'w','linewidth',1);xlog,hold on
contour(p./pi,alpha,be,[1 1],'k','linewidth',2);
title('Properties of Generalized Morse Wavelets')
text(0.55,-2,'(b)  $\beta$ Parameter')

subplot(2,3,3),contour(p./pi,alpha,a,(.52:.01:.59),'k');
xlog,hold on
contour(p./pi,alpha,a,[.51 0.51],'k','linewidth',2);
text(0.55,-2,'(c)  Heisenberg Area')

subplot(2,3,4),
x=(frac(fe-fm,fm));x(x>0.2)=0.2;
contour(p./pi,alpha,x,(0.025:.025:.2),'k'),hold on,xlog
x=-(frac(fe-fm,fm));x(x>0.2)=0.2;
contour(p./pi,alpha,x,(0.025:.025:.2),'k','linewidth',2),%
contour(p./pi,alpha,x,(0.025:.025:.2),'w:','linewidth',1),%
contour(p./pi,alpha,x,[0 0],'k','linewidth',2),%
ylabel('Demodulate Skewness $\Im\{\alpha_{3;\beta,\gamma}\}$')
xlabel('Duration $P_{\beta,\gamma}/\pi$')
%hold on,nocontours
text(0.55,-2,'(d)  Energy Freq. / Peak Freq. - 1 ')

subplot(2,3,5),
x=(frac(fi-fm,fm));x(x>0.2)=0.2;
contour(p./pi,alpha,x,(0.025:.025:.2),'k'),hold on,xlog
x=-(frac(fi-fm,fm));x(x>0.2)=0.2;
contour(p./pi,alpha,x,(0.025:.025:.2),'k','linewidth',2),%
contour(p./pi,alpha,x,(0.025:.025:.2),'w:','linewidth',1),%
contour(p./pi,alpha,x,[0 0],'k','linewidth',2),%
xlabel('Duration $P_{\beta,\gamma}/\pi$')
text(0.55,-2,'(e)  Inst. Freq. / Peak Freq. - 1')

subplot(2,3,6),
x=frac(1,2*pi)*(cf);x(x>0.2)=0.2;
contour(p./pi,alpha,x,(0.025:.025:.2),'k'),hold on,xlog
x=-frac(1,2*pi)*(cf);x(x>0.2)=0.2;
contour(p./pi,alpha,x,(0.025:.025:.2),'k','linewidth',2),%
contour(p./pi,alpha,x,(0.025:.025:.2),'w:','linewidth',1),%
contour(p./pi,alpha,x,[0 0],'k','linewidth',2),%
xlabel('Duration $P_{\beta,\gamma}/\pi$')
text(0.55,-2,'(f)  Frequency Curvature')

for i=1:6
    subplot(2,3,i)
    xlim([10^(-1/2) 10^1.5])
    xtick([1/2 1 2 4 8 16 32])
    %if ~verLessThan('matlab','8.4.0')
    %    %Fix weird labelling 
    %    h=get(gca,'xticklabel');
    %    h{2}='1';
    %    set(gca,'xticklabel',h)
    %end
    hlines(0,'k:')
    [p1,a1]=morseprops(1/100,[(0:.1:100) (101:1:1000) (1000:100:1000000)] );
    plot(p1./pi,a1,'k--','linewidth',1);
end

packfig(2,3)

orient landscape
set(gcf,'paperposition',[2  2 8 5])
fontsize 12 10 10 9
if strcmpi(str,'print')
    print -depsc morsie_parameters.eps
end
%\*****************************************