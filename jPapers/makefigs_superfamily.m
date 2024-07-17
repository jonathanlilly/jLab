function[varargout]=makefigs_superfamily(str)
%MAKEFIGS_SUPERFAMILY  Make figures for Lilly and Olhede (2012b).
% 
%   This function makes all figures for the paper 
%
%   Lilly, J. M., and S. C. Olhede (2012). Generalized Morse wavelets as a
%      superfamily of analytic wavelets.  IEEE Transactions on Signal
%      Processing 60 (11), 6036--6041.
%
%   Usage: makefigs_superfamily
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2012--2023 J.M. Lilly --- type 'help jlab_license' for details

if nargin==0
    str='noprint';
end

%/************************************************************
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
        plot(t./(P(i,j)./pi),x),xlim([-2 2]),ylim([-0.85 1.1])
        linestyle k k-- 1.5k
        hlines(0,'k:'),vlines([-1/2 1/2],'k:')
        if iseven(i+(5-j+1)*5),set(gca,'color',[1 1 1]*0.9),end
        if j==1,title(['$\gamma=$' num2str(ga(i),3)]),end
        if i==1,ylabel(['$\beta=$' num2str(be(j),3)]),end
        set(gca,'xtick',[]),set(gca,'ytick',[]),set(gca,'xcolor','w'),set(gca,'ycolor','w')
        set(get(gca,'ylabel'),'color','k')
        if j==5&&i==1,text(-1.8,-0.6,'(a)'),end
    end
end
packfig(5,5)

fontsize 14 14 14 14
orient tall
set(gcf,'paperposition',[1 1 8 8])

if strcmpi(str,'print')
    print -deps morsie_families_new.eps
end
%\*********************************************************


%/************************************************************
figure

N=2^15;
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


f=frac(1,2*pi)*fourier(N);
for i=1:5
    for j=1:5
        subplot(5,5,i+(j-1)*5)
        psigauss=simplepdf(f,fs,fs./P(i,j),'gaussian');
        psigauss=psigauss./maxmax(psigauss).*2;
        fbg=-frac(1,6)*P(i,j).^2.*(f./fs-1).^3.*(ga(i)-3)-frac(1,24)*P(i,j).^2.*(f./fs-1).^4.*((ga(i)-3).^2+2);
        plot(f,[2*psif(1:length(f),i,j) psigauss psigauss.*exp(fbg)]),xlim([0 f(129+64)]),ylim([0 2.1])
        linestyle 2k k-.  k 
        vlines(1/512,'k:')
        if iseven(i+(5-j+1)*5),set(gca,'color',[1 1 1]*0.9),end
        if j==1,title(['$\gamma=$' num2str(ga(i),3)]),end
        if i==5,ylabel(['$\beta=$' num2str(be(j),3)]),end
        %if i==1,vlines(0,'k'),end
        if j==5&&i==1,text(.0005,0.25,'(b)'),end
    end
end
h=packfig(5,5);
h=reshape(h,5,5);
for i=1:5
    for j=1:5
        axes(h(i,j))
        if i==5,ylabel(['$\beta=$' num2str(be(j),3)]),end
        set(gca,'xtick',[]),set(gca,'ytick',[]),set(gca,'xcolor','w'),set(gca,'ycolor','w')
        set(get(gca,'ylabel'),'color','k'),set(gca,'yaxislocation','right')
    end
end

fontsize 14 14 14 14
orient tall
set(gcf,'paperposition',[1 1 8 8])

if strcmpi(str,'print')
    print -deps morsie_families_fourier.eps
end
%\*********************************************************




%/*********************************************************************
ga1=logspace(log10(.1),2,100);
be1=logspace(log10(1/4+.001),2,101);

[ga,be]=meshgrid(ga1,be1);
[fm,fe,fi,cf] = morsefreq(ga,be);
[a,sigt,sigo,skew]=morsebox(ga,be);
p=sqrt(be.*ga);

figure
contourf(ga1,be1,a,[1/2:.05:1.5]),xlog,ylog,nocontours
hold on, contour(ga1,be1,a,(.500:.005:.55),'k'),colormap gray,flipmap
hold on, contour(ga1,be1,skew,[0 0],'k--'),colormap gray,flipmap

vlines([1 2 3],'k:'),
plot(ga1,frac(ga1-1,2),'k','linewidth',2)
for i=-1:3
    h=plot(ga1,(3.^(2*i))./ga1);linestyle -h h 2w
    h=plot(ga1,(3.^(2*i))./ga1);linestyle -h h E
end



%contour(ga1,be1,(k3./k2.^(3/2)),[0 0],'k')
caxis([1/2 1.5])
axis([1/10 100 1/4 100])

xtick([1/2 1 2 3 4 5 6 8 16 32])
ytick([1/2 1 2 4 8 16 32 64 128])

if verLessThan('matlab','8.4.0')
    xtl=get(gca,'xticklabel');xtl(1,:)='1/2';set(gca,'xticklabel',xtl)
    ytl=get(gca,'yticklabel');ytl(1,:)='1/2';set(gca,'yticklabel',ytl)
else
    xtl=get(gca,'xticklabel');xtl{1}='1/2';set(gca,'xticklabel',xtl)
    ytl=get(gca,'yticklabel');ytl{1}='1/2';set(gca,'yticklabel',ytl)
end
outticks

title('Generalized Morse Wavelet Phase Diagram')
xlabel('Gamma Parameter'),ylabel('Beta Parameter')

text(.12+.25,17.86,'B')
text(1.25,90,'C')
text(2.25,90,'G')
text(3.25,90,'A')
text(1/2,85,'e')
text(64,0.28,'S')
text(.30,.28,'a')

ga=[1/3 1 3 9 27]';
for i=1:5,for j=1:5,plot(ga(i),ga(j),'wo','markersize',4,'markerfacecolor','k'),end,end
plot(1/10,22,'ks','markersize',5,'markerfacecolor','k')

hc=colorbar;hc=discretecolorbar(hc,[.5 1.5],(1/2:.05:1.5));nocontours
axes(hc),hlines(.500:.005:.55,'k')
ylabel('Heisenberg Area'),ytick(1/2:.1:1.5),fixlabels([0 -1])



orient tall
fontsize 12 12 12 12
set(gcf,'paperposition',[1 1 5 5])


if strcmpi(str,'print')
   print -depsc analytic-morsephase.eps
end
%/************************************************************************


%/********************************************************************
%Compute Morlet area and projection onto Gaussian
disp('Sorry, this next bit takes a while....')

om=logspace(log10(1.1),log10(2.5*pi),100);
nu=morlfreq(om);
p_morlet=om.*sqrt(om.*(om-nu)+1);
N=2^13;
t=[0:N-1]';t=t-mean(t);
L=50;

clear sigmat sigmao psi g
for i=1:length(om)
        [psii,psifi]=morlwave(N,om(i),1./L,'bandpass');
        [mut,sigmat(i)]=pdfprops(t,abs(psii).^2); 
        [muo,sigmao(i)]=pdfprops(2*pi*fftshift(t./length(t)),abs(psifi).^2); 
        %Demodulated version
        psi(:,i)=morlwave(N,om(i),1/L,'energy').*rot(-1/L.*t);
        g(:,i)=simplepdf(t,0,p_morlet(i)*L,'gaussian');
        g(:,i)=g(:,i)./sqrt(sum(abs(g(:,i)).^2));
end
proj_morlet=abs(squeeze(vsum(conj(psi).*g,1)));


%ga1=[2.99 3 3.01];
ga1=(1:1:6);
be1=logspace(-10,2,200);
[ga,be]=meshgrid(ga1,be1);
[fm,fe,fi,cf] = morsefreq(ga,be);
a=morsebox(ga,be);
[p,skew,kurt]=morseprops(ga,be);

clear psi_morse g
for i=1:length(be1)
    for j=1:length(ga1)
        psi_morse(:,i,j)=morsewave(N,1,ga1(j),be1(i),1/L,'energy').*rot(-1/L.*t);
        g(:,i,j)=simplepdf(t,0,p(i,j)*L,'gaussian');
        g(:,i,j)=g(:,i,j)./sqrt(sum(abs(g(:,i,j)).^2));
    end
end
proj_morse=abs(squeeze(vsum(conj(psi_morse).*g,1)));

figure, 
%subplot(2,1,1)
plot(p/pi,1./vswap(a,inf,nan))
hold on, plot(p_morlet/pi,1./real(sigmat.*sigmao))
linestyle k k-- 2k 2k-- k-- k-. 4D
xlim([0 2.5]),ylim([0 2])
vlines(sqrt([1:6]/2)/pi,'D:')
fixlabels(-1)

hlines(1./0.55,'D')
title('Wavelet Inverse Heisenberg Area')
xlabel('Duration $P_{\beta,\gamma}/\pi$')
ylabel('Inverse Area $1/A_{\beta,\gamma}$')
text(0.07,0.1,'(a)')
legend('$\gamma=1$','$\gamma=2$','$\gamma=3$','$\gamma=4$','$\gamma=5$','$\gamma=6$','Morlet','location','southeast');

set(gcf,'paperposition',[1 1 5 4])
fontsize 12 12 12 12

if strcmpi(str,'print')
     print -depsc wavelets_heisenberg_area.eps
end

figure
plot(p./pi,proj_morse.^2),xlim([0 2.5]),ylim([.92 1])
hold on, plot(p_morlet./pi,proj_morlet.^2)
linestyle k k-- 2k 2k-- k-- k-. 4D

title('Wavelet Envelope Versus Gaussian')
xlabel('Duration $P_{\beta,\gamma}/\pi$')
ylabel('Squared Inner Product')
ytick([.92:.01:1]),fixlabels([-1 -2])

text(0.1,0.925,'(b)')

legend('$\gamma=1$','$\gamma=2$','$\gamma=3$','$\gamma=4$','$\gamma=5$','$\gamma=6$','Morlet','location','southeast');

set(gcf,'paperposition',[1 1 5 4])
fontsize 12 12 12 12

if strcmpi(str,'print')
     print -depsc wavelets_versus_gaussian.eps
end
%\********************************************************************


%END of jlab_makefigs_superfamily