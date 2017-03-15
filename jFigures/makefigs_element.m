function[varargout]=makefigs_element(str)
%MAKEFIGS_ELEMENT  Makes all figures for Lilly (2017), "Element Analysis."

%cd /Users/lilly/Desktop/Dropbox/Projects/impulses

if nargin==0
  str='print';
end

if strcmp(str,'--f')
     makefigs_impulses('noprint');return
end

%/************************************************************************
be=[0 1/4 1/2 1 2 4];
ga=[1 2 3 4];

M=length(ga);
N=length(be);
t=[1:5000]';
t=t-mean(t);

figure
for i=1:length(be)
    for j=1:length(ga)    
        subplot(M,N,i+(j-1)*N)
        psi=morsewave(length(t),1,ga(j),be(i),1/100,'bandpass');
        p=morseprops(ga(j),be(i));
        if be(i)==0
            p=sqrt(2)*pi./morsefreq(ga(j),be(i))/(2*sqrt(2));
        end
        psi=psi./maxmax(psi);
        uvplot(t./p/50,psi),hold on
        plot(t./p/50,abs(psi))
        linestyle k k-- 2k
        axis off,axis([-10 10 -0.85 1.1])
        hlines(0,'k:')
        vlines(100*sqrt(2)/50*[1 -1],'k:')  %This is L_{\beta,\gamma}
        %plot(t./p/50,squared(psi),'r')
        if i==1
            vlines(10,'G')
        end
        if i==2
            vlines(-10,'G')
        end
        if i==1 
            text(-9,.9,['$\gamma$ = ' int2str(ga(j))])
        end
        if j==1
            text(4,.9,['$\beta$ = ' num2str(be(i))])
        end
        if i==3 && j==1
            title('$\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,$ Time-Domain Forms of Narrow Generalized Morse Wavelets')
            set(gca,'TitleFontWeight','Normal')
        end
    end
end
packfig(M,N,'both')
fontsize 12 10 10 10
orient landscape
set(gcf,'paperPosition',[.25 .25 10.5 6])  

if strcmpi(str,'print')
    %print -djpeg -r1000 timedomainwavelets.jpg
    print -djpeg timedomainwavelets.jpg
    crop('timedomainwavelets.jpg')
end
%\************************************************************************


%/*********************************************************************
N=12000;
t=[1:N]';t=t-mean(t);
ga=2;be=2;mu=1;
fo=logspace(log10(2*pi/100),log10(2*pi/1000),6)';
rho=morsefreq(ga,mu)./fo;
psi=morsewave(N,ga,mu,fo);
tau=[100:200:1200]'*10;
phi=linspace(0,pi/2,6)';
c=zeros(size(phi));

for i=1:size(psi,2)
    c(i)=2*rot(phi(i))./maxmax(abs(psi(:,i)));
    psi(:,i)=real(c(i).*vshift(psi(:,i),N/2-tau(i),1));
end

xo=sum(psi,2);
fs=morsespace(ga,be,{0.05,pi},{3,N});
%--------------------------------------------------------------------------
%Clean example
wo=wavetrans(xo,{ga,be,fs},'mirror');
[index,ww,ff]=transmax(fs,wo);
[ii,jj]=ind2sub(size(wo),index);
[chat,rhohat,fhat]=maxprops(ww,ff,ga,be,mu);
make impulses_cleanexample N t xo wo index ii jj ww ff chat fhat rhohat
%--------------------------------------------------------------------------
%Noisy examples
rng(1);
for alpha=[0 1]
    if alpha==0
        xe=randn(N,1);
    elseif alpha==1
        xe=cumsum(randn(N,1));
        xe=xe./std(xe);
        xe=xe-mean(xe);
    end
    x=xo+xe;    
    we=wavetrans(xe,{ga,be,fs},'mirror');
    w=wavetrans(x,{ga,be,fs},'mirror');
    [index,ww,ff,rr]=transmax(fs,w,{ga,be,false(size(x))});
    [ii,jj]=ind2sub(size(w),index);
    [chat,rhohat,fhat]=maxprops(ww,ff,ga,be,mu);
    [bool,z]=isomax(size(w),index,ww,ff,ga,be,mu,1/2);
    
    if alpha==0
        clear impulses_whitenoiseexample
        make impulses_whitenoiseexample N t x xo xe w wo we index ii jj ww ff rr chat fhat rhohat bool z
    elseif alpha==1
        clear impulses_rednoiseexample
        make impulses_rednoiseexample N t x xo xe w wo we index ii jj ww ff rr chat fhat rhohat bool z
    end
end
%--------------------------------------------------------------------------
%Transform maxima of a clean signal
use impulses_cleanexample

psihat=morsewave(N,ga,mu,fhat);
for i=1:length(ii)
    psihat(:,i)=real(chat(i).*vshift(psihat(:,i),N/2-ii(i),1)).*rhohat(i);
end
xhat=sum(psihat,2);

% figure, [h,hl]=wavespecplot(t,[xhat xo],2*pi./fs,wo);colormap lansey
% linestyle -hh hl k: k
% axes(h(1)),ylabel('Signal Amplitude'),
% title('Wavelet Transform of Isolated Events')
% text(-5700,-1.80,'(a)')
% axes(h(2)),ylabel('Transform Period $2\pi/\omega_s$'),xlabel('Time')
% plot(t(ii),2*pi./ff,'ws','markerfacecolor',[1 1 1]*0.7,'markersize',6)
% text(-5700,2.6*10^3,'(b)','color','w'),caxis([0 maxmax(abs(wo))])
% 
% [bool,z]=isomax(size(wo),index,ww,ff,ga,be,mu,1/2);
% plot(interp1([1:length(t)]',t,real(z))+1i*2*pi./imag(z),'color',[1 1 1]*0.7)
% 
% orient portrait
% fontsize 11 11 11 11
% set(gcf,'paperposition',[1 1 8 5])
% if strcmpi(str,'print')
%    %print -djpeg -r1000 impulses-morsetrain-clean.jpg
%    print -djpeg impulses-morsetrain-clean.jpg
%    crop('impulses-morsetrain-clean.jpg')
% end
%--------------------------------------------------------------------------
%Transform maxima of a clean signal plus background signals
use impulses_cleanexample

%f1=2*pi./10.^[1 1.5 2 2.5 3];
%f1=2*pi./10.^[1 2 3];
f1=2*pi./10.^[2 3];
iieps=[];feps=[];
for i=1:length(f1)
    %ii1=round(linspace(1,N,1000*f1(i)/f1(1)))';
    ii1=round(linspace(1,N,100*f1(i)/f1(1)))';
    iieps=[iieps;ii1];
    feps=[feps;f1(i)+0*ii1];
end
rhoeps=morsefreq(ga,mu)./feps;

xeps=zeros(size(x));
for i=1:length(iieps)
    psi=morsewave(N,ga,mu,feps(i));
    psi=2*psi./maxmax(abs(psi))/10;
    xeps=xeps+real(rot(rand(1)*2*pi)*vshift(psi,N/2-iieps(i),1));
end

x=xo+xeps;

w=wavetrans(x,{ga,be,fs},'mirror');
[index,ww,ff]=transmax(fs,w);
[ii,jj]=ind2sub(size(w),index);
[chat,rhohat,fhat]=maxprops(ww,ff,ga,be,mu);

figure, [h,hl]=wavespecplot(t,[x xo],2*pi./fs,w);colormap lansey
linestyle -h hl D U
axes(h(1)),ylabel('Signal Amplitude'),
title('Wavelet Transform of Isolated Events Plus Densely Spaced Events')
%hold on,plot(t,xhat,'w--')
text(-5700,-1.80,'(a)')
axes(h(2)),ylabel('Transform Period $2\pi/\omega_s$'),xlabel('Time')
fact=frac(morsefreq(ga,mu),morsefreq(ga,be)).*frac(be,mu+1).^(1./ga);
plot(t(iieps),2*pi./(feps/fact),'wo','markersize',4)
plot(t(ii),2*pi./ff,'.','markeredgecolor',[1 1 1]*0.7)
plot(t(ii(abs(ww)>1)),2*pi./ff(abs(ww)>1),'ko','markerfacecolor','k')
%plot(t(ii),2*pi./ff,'ws','markerfacecolor',[1 1 1]*0.7,'markersize',6)
text(-5700,2.6*10^3,'(b)','color','w'),caxis([0 maxmax(abs(wo))])

[bool,z]=isomax(size(wo),index,ww,ff,ga,be,mu,1/2);
plot(interp1([1:length(t)]',t,real(z))+1i*2*pi./imag(z),'color',[1 1 1]*0.7)
[bool,z]=isomax(size(wo),index(abs(ww)>1),ww(abs(ww)>1),ff(abs(ww)>1),ga,be,mu,1/10);
plot(interp1([1:length(t)]',t,real(z))+1i*2*pi./imag(z),'k:')

L=2*sqrt(2)*sqrt(ga*be)./fs;
h=plot([t(1)+L/2 t(end)-L/2],2*pi./fs);linestyle -h h 2E


orient portrait
fontsize 11 11 11 11
set(gcf,'paperposition',[1 1 8*1.5 5*1.5])
if strcmpi(str,'print')
    %print -djpeg -r1000 impulses-morsetrain-clean.jpg
    print -djpeg impulses-morsetrain-dense.jpg
    crop('impulses-morsetrain-dense.jpg')
end
%--------------------------------------------------------------------------
%Figure of noisy transforms
for alpha = [0 1]
    if alpha==0
        use impulses_whitenoiseexample
    elseif alpha==1
        use impulses_rednoiseexample
    end
    
    sigma=sqrt(vmean(vcolon(squared(w(:,1))),1));
    wwtilde=ww./(sigma.*sqrt(ff./fs(1)));
    
    psihat=morsewave(N,ga,mu,fhat);
    for i=1:size(psihat,2)
        psihat(:,i)=real(chat(i).*vshift(psihat(:,i),N/2-ii(i),1)).*rhohat(i);
    end
    
    if alpha==0
        bool1=abs(wwtilde)>4;
    elseif alpha==1
        bool1=abs(ww)>1; 
    end
    
    sigbool=bool&bool1;
    xhat=sum(psihat(:,sigbool),2);
    xhat2=sum(psihat(:,bool1),2);
 
    figure
    if alpha==0
        [h,hl]=wavespecplot(t,[x xo xhat2 xhat],2*pi./fs,w);colormap lansey
        linestyle -h hl D 3k k: U
    elseif alpha==1
        [h,hl]=wavespecplot(t,[xo xhat2 xhat x],2*pi./fs,w);colormap lansey
        linestyle -h hl 3k k: U D
    end
    
    axes(h(1)),ylabel('Signal Amplitude'),ylim([-4.2 4.2])
    if alpha ==0
        title('Wavelet Transform of Isolated Events Plus White Noise')
    else
        title('Wavelet Transform of Isolated Events Plus Red Noise')
    end
    text(-5700,-3.7,'(a)'),
    axes(h(2)),ylabel('Transform Period $2\pi/\omega_s$'),xlabel('Time')
    %plot(t(ii),2*pi./ff,'w.','markersize',8)
    use impulses_cleanexample
    plot(t(ii),2*pi./ff,'ws','markersize',4)
    
    if alpha==0
        use impulses_whitenoiseexample
    elseif alpha==1
        use impulses_rednoiseexample
    end
    
    vindex(index,ii,jj,ff,ww,sigbool,1);
    plot(t(ii),2*pi./ff,'wo','color','k','markerfacecolor','k','markersize',4)
    if alpha==0
        text(-5700,2.6*10^3,'(b)','color','w'),
    elseif alpha==1
        text(-5700,2.6*10^3,'(b)','color','k'),
    end
    
    caxis([0 maxmax(abs(wo))])
    plot(interp1([1:length(t)]',t,real(z(:,sigbool)))+1i*2*pi./imag(z(:,sigbool)),'color',[1 1 1]*0)
    
    if alpha==0
        use impulses_whitenoiseexample
    elseif alpha==1
        use impulses_rednoiseexample
    end
    
    plot(t(ii(~sigbool)),2*pi./ff(~sigbool),'.','color',[1 1 1]*0.7)
    L=2*sqrt(2)*sqrt(ga*be)./fs;
    h=plot([t(1)+L/2 t(end)-L/2],2*pi./fs);linestyle -h h 2E
    
    orient portrait
    fontsize 11 11 11 11
    set(gcf,'paperposition',[1 1 8*1.5 5*1.5])

    if strcmpi(str,'print')
        if alpha ==0  
            %print -djpeg -r1000 impulses-morsetrain-noisy.jpg
            print -djpeg impulses-morsetrain-whitenoisy.jpg
            crop('impulses-morsetrain-whitenoisy.jpg')
            %print -dpng impulses-morsetrain-noisy-bw.png
            %crop('impulses-morsetrain-noisy-bw.png')
        elseif alpha ==1
            %print -djpeg -r1000 impulses-morsetrain-noisy.jpg
            print -djpeg impulses-morsetrain-rednoisy.jpg
            crop('impulses-morsetrain-rednoisy.jpg')
            %print -dpng impulses-morsetrain-noisy-bw.png
            %crop('impulses-morsetrain-noisy-bw.png')
        end
    end
end

%--------------------------------------------------------------------------
%Fourier and Wavelet Spectra
[psi,lambda]=sleptap(N,20); 
[f,s]=mspec([x xo xe],psi,lambda,'adaptive');

figure
subplot(1,2,1)
plot(2*pi./f,s),xlog,ylog,hlines(2)
xlim(2*pi./[max(f) f(2)]),ylim([10^(-2) 100])
title('Fourier Spectrum')
xlabel('Period')
linestyle k 3k 4D k--
h=legend('Signal + Noise','Signal Only','Noise Only','Predicted Noise');
plot(2*pi./f,s(:,1),'k')
%pos=get(h,'Position');
set(h,'Position',[0.3213    0.2550    0.1215    0.1502])
text(3,10^(-1.75),'(a)')
ylabel('Fourier Power Spectral Density')

subplot(1,2,2)
[mo,ffun]=morsemom(0,ga,be);
a=ffun.*morsefreq(ga,be).^(-1);
plot(2*pi./fs,[mean(squared(w),1)' mean(squared(wo),1)' mean(squared(we),1)'  a*fs])
xlog,ylog,xlim(2*pi./[max(f) f(2)]),ylim([10^(-4) 1])
linestyle k 3k 4D k--
hold on,plot(2*pi./fs,mean(squared(w),1),'k')
title('Wavelet Spectrum with 1/s Normalization')
xlabel('Transform Period $2\pi/\omega_s$')
text(3,10^(-3.75),'(b)')
set(gca,'yaxislocation','right')
packfig(1,2,'columns')
set(gca','YTickLabelMode','auto')
ylabel('Averaged Squared Wavelet Transform Modulus')

orient portrait
fontsize 14 12 12 12 
set(gcf,'paperposition',[1 1 11 4.5])

if strcmpi(str,'print')
   %print -djpeg -r1000 impulses-spectra.jpg
   print -djpeg impulses-spectra.jpg
   crop('impulses-spectra.jpg')
end
%--------------------------------------------------------------------------
%Computation of transform maxima for (2,2) wavelet only
%Note:  This figure takes about 5 minutes on a powerful machine
M=10000;  %Reduce this to 1000 if it is too big for your machine
N=12000;
ga=2;be=2;mu=1;
for alpha=[0 1]

    fs=morsespace(ga,be,{0.05,pi},{3,N});
    tic;[count,bins,rate]=transmaxdist(ga,be,alpha,fs,fs(1)./fs(2),N,M);toc
    
    %Scale these by the wavelet footprint L
    counttilde=count;ratetilde=rate;
    for i=1:length(fs)
        L=2*sqrt(2)*sqrt(be*ga)./fs(i);
        ratetilde(:,i)=rate(:,i)./(N./L);
        counttilde(:,i)=counttilde(:,i)./(N./L);
    end
    %[mu,sigma,skew]=pdfprops(bins,counttilde(:,1));  %mu=1.36, skew =0.34
    
    %Same thing but for noise simulation
    M=2000;
    if alpha==0
        x=randn(N*M,1);
    elseif alpha==1
        x=cumsum(randn(N*M,1));
        x=x./std(x);
    end
    tic;w=wavetrans(x,{ga,be,fs(1:3)});toc
    [index,ww]=transmax(fs(1:3),w);[ii,jj]=ind2sub(size(w),index);
    sigma=vstd(w(:,2),1);
    n=hist(abs(ww)./sigma,bins)';
    L=2*sqrt(2)*sqrt(be*ga)./fs(2);
    ntilde=n./(N*M./L);
    nrate=cumsum(n,1,'reverse')./(N*M./L);
    
    figure,subplot(1,2,1)
    h1=plot(bins,counttilde(:,26:end-1)*1e4);linestyle -h h1 G,hold on
    h2=plot(bins,counttilde(:,3:25)*1e4);linestyle -h h2 C,hold on
    h3=plot(bins,counttilde(:,2)*1e4);linestyle -h h3 4w,hold on
    h4=plot(bins,counttilde(:,2)*1e4);linestyle -h h4 k,hold on
    h5=plot(bins,ntilde*1e4,'k.');
    xlim([0 2.999]),set(gca,'ytickmode','auto')
    if alpha==0
        ylim([0 13])
        title('Event Histograms for White Noise with $\psi_{2,2}(t)$')
    elseif alpha==1
        ylim([0 17])
        title('Event Histograms for Red Noise with $\psi_{2,2}(t)$')
    end
    xlabel('Normalized Event Value')
    ylabel('Normalized Event Density ($\times 10^{-4}$)')
    
    subplot(1,2,2)
    h=plot(bins,ratetilde(:,26:end-1)*1e3);linestyle -h h G,hold on
    h=plot(bins,ratetilde(:,3:25)*1e3);linestyle -h h C,hold on
    h=plot(bins,ratetilde(:,2)*1e3);linestyle -h h 4w,hold on
    h=plot(bins,ratetilde(:,2)*1e3);linestyle -h h k,hold on
    h=plot(bins,nrate*1e3,'k.')
    xlim([0 3])
    index=find(nrate(:,1)<1/100,1,'first');
    if alpha==0
        ylim([0 45])
        title('Detection Rates for White Noise with $\psi_{2,2}(t)$')
    elseif alpha==1
        ylim([0 58])
        title('Detection Rates for Red Noise with $\psi_{2,2}(t)$')
    end
    vlines(bins(index),'k:'),hlines(1/100*1e3,'k:')
    xlabel('Normalized Event Value')
    letterlabels(2)
    subplot(1,2,1),legend([h4(1) h2(1) h1(1) h5],'Scale $\#2$','Scales $\#3-25$','Scales $\#26-53$','Scale $\#2$ W.T.')
    subplot(1,2,2)
    h=packfig(1,2,'columns');
    set(gca,'YAxisLocation','Right')
    set(gca,'YTickLabelMode','Auto'),ylabel('Detection Rate per $L_{\beta,\gamma}(s)$  ($\times 10^{-3}$)')
    axes(h(1)),set(gca,'ytickmode','auto')
    
    orient portrait
    fontsize 14 12 12 12
    set(gcf,'paperposition',[1 1 11 5])
    
    if strcmpi(str,'print')
        if alpha==0
            %print -djpeg -r1000 impulses-noisedist-22.jpg
            print -djpeg impulses-noisedist-22.jpg
            crop('impulses-noisedist-22.jpg')
        elseif alpha==1
            %print -djpeg -r1000 impulses-noisedist-22.jpg
            print -djpeg impulses-rednoisedist-22.jpg
            crop('impulses-rednoisedist-22.jpg')
        end
    end
end  
%--------------------------------------------------------------------------
%Next, determine noise distribution for this particular example
%Have to create the contours in a roundabout way because of Matlab colormap issues
%/*********************************************************************
N=12000;
ga=2;be=2;
fs=morsespace(ga,be,{0.05,pi},{3,N});
for alpha=[0 1]

    %tic;[count,bins,rate]=transmaxdist(ga,be,0*ga,fs,fs(1)./fs(2),12000,1000,'parallel');toc
    tic;[count,bins,rate]=transmaxdist(ga,be,alpha+0*ga,fs,fs(1)./fs(2),12000,1000,'parallel','extrap');toc
    
    if alpha==0
        use impulses_whitenoiseexample
    elseif alpha==1
        use impulses_rednoiseexample
    end
    
    ci=[1/1000 1/100 1/10 1 10];
    [m0,ffun]=morsemom(alpha,ga,be);
    [fmat,xmat]=meshgrid(fs,bins);
    a=ffun.*frac(morsefreq(ga,be),fmat).^(2*alpha-1);
    
    frat=frac(fs(1),morsefreq(ga,be)).^(2*alpha-1);
    
    %This is to set the noise amplitude A based on the variance in the 
    %highest wavelet band.  It is only necessary to use for red noise.
    A=sqrt(vmean(squared(we(:,1)),1)*frat./ffun);
    
    figure
    clear xi1 yi1 xi2 yi2
    for i=1:length(ci)
        [h,hc]=contour(2*pi./fmat,xmat.*A.*sqrt(a),rate,[ci(i) ci(i)],'k--');
        npoints=hc.ContourMatrix(2,1);
        xi1{i}=hc.ContourMatrix(1,2:npoints+1)'; yi1{i}=hc.ContourMatrix(2,2:npoints+1)';
        [h,hc]=contour(2*pi./fs,bins,rate,[ci(i) ci(i)],'k--');
        npoints=hc.ContourMatrix(2,1);
        xi2{i}=hc.ContourMatrix(1,2:npoints+1)'; yi2{i}=hc.ContourMatrix(2,2:npoints+1)';
    end
    close
    
    [boolo,z]=isomax(size(w),index,ww,ff,ga,be,0,1/2);
    
    figure,
    subplot(1,2,1)
    contourf(2*pi./fmat,xmat.*A.*sqrt(a),count,50),xlog,hold on,colormap lansey,nocontours
    map=colormap;map(1,:)=[1 1 1];colormap(map);

    wwtilde=ww./sqrt(ffun.*frac(morsefreq(ga,be),ff).^(2*alpha-1));
    
    if alpha==0
         bool=boolo&(abs(wwtilde)>4);
    elseif alpha==1
         bool=boolo&(abs(ww)>1);
    end
  
    plot(2*pi./ff(bool),abs(ww(bool)),'ko','markerfacecolor','k'), hold on
    bool=~boolo;
    plot(2*pi./ff(bool),abs(ww(bool)),'w.','markersize',8),xlim(2*pi./[1.4215 0.0025 ])
    plot(2*pi./ff(bool),abs(ww(bool)),'k.','color',[1 1 1]*0.5),
   
    use impulses_cleanexample
    plot(2*pi./ff,abs(ww),'s','color',[1 1 1]*0.3,'markersize',8)
    xlabel('Transform Period')
    ylabel('Event Magnitude')
    h=cellplot(xi1,yi1);ylim([0 2])
    linestyle -h h 2k k-- k 3D-- 3D
    
    if alpha==0
        title('Event Distribution for Example with White Noise')
    elseif alpha==1
        title('Event Distribution for Example with Red Noise')
    end
    
    text(1.5*10^3,1.9,'(a)')
    
    subplot(1,2,2)
    contourf(2*pi./fmat,xmat,count,50),xlog,hold on,colormap lansey,nocontours
    map=colormap;map(1,:)=[1 1 1];colormap(map);xlim(2*pi./[1.4215 0.0025 ])
    
    if alpha==0
        use impulses_whitenoiseexample
    elseif alpha==1
        use impulses_rednoiseexample
    end
    
    wwtilde=ww./sqrt(ffun.*frac(morsefreq(ga,be),ff).^(2*alpha-1))./A;

    if alpha==0
        bool=boolo&(abs(wwtilde)>4);
    elseif alpha==1
        bool=boolo&(abs(ww)>1);
    end
    
    plot(2*pi./ff(bool),abs(wwtilde(bool)),'ko','markerfacecolor','k'), hold on
    bool=~boolo;
    plot(2*pi./ff(bool),abs(wwtilde(bool)),'w.','markersize',8),xlim(2*pi./[1.4215 0.0025 ])
    plot(2*pi./ff(bool),abs(wwtilde(bool)),'k.','color',[1 1 1]*0.5),
    
    use impulses_cleanexample
    wwtilde=ww./sqrt(ffun.*frac(morsefreq(ga,be),ff).^(2*alpha-1))./A;
    plot(2*pi./ff,abs(wwtilde),'s','color',[1 1 1]*0.3,'markersize',8)
    if alpha==0 
        ylim([0 20])
        text(1.5*10^3,19,'(b)')
    else
        ylim([0 15])
        text(1.5*10^3,14.3,'(b)')
    end
    
    h=cellplot(xi2,yi2);
    linestyle -h h 2k k-- k 3D-- 3D
    xlabel('Transform Period')
    ylabel('Normalized Event Magnitude')
    title('Event Distribution, Normalized Magnitude')
    legend(h,' rate: 1/1000',' rate: 1/100',' rate: 1/10','rate: 1',' rate: 10','location','best')
    
    packfig(1,2,'columns')
    set(gca,'YAxisLocation','Right')
    set(gca,'YTickLabelMode','Auto')
    ylabel('Normalized Event Magnitude')
    
    orient portrait
    fontsize 14 12 12 12
    set(gcf,'paperposition',[1 1 11 5])
    
    if strcmpi(str,'print')
        if alpha==0
            %print -djpeg -r1000 impulses-examplewhitenoisedist.jpg
            print -djpeg impulses-examplewhitenoisedist.jpg
            crop('impulses-examplewhitenoisedist.jpg')
        elseif alpha==1
            %print -djpeg -r1000 impulses-examplerednoisedist.jpg
            print -djpeg impulses-examplerednoisedist.jpg
            crop('impulses-examplerednoisedist.jpg')
        end
    end
end

%Verify that 10x matches prediction
% temp=xmat./sqrt(2).*sqrt(fmat.*a);
% use impulses_noiseexample 
% clear Ni
% for i=1:size(nmatcum,2)
%     index=find(nmatcum(:,i)<=10,1,'first');
%     Ni(i)=length(find((round(jj)==i)&(abs(ww)>temp(index,i))));
% end
% plot(fs/2/pi,Ni,'.'),hlines(10)%Great, looks fine! 
%--------------------------------------------------------------------------
%\*************************************************************************

%/*************************************************************************
figure
clear fs w 
mu=0;
N=1000;
frho=2*pi/100;
alpha=[0.25 0.5 0.75 0.85 0.95];
betas=[0.25 0.5 1 2 4];
A=[1 10 20 40 100];
t=[1:N]';t=t-mean(t);


for i=1:5
    be=betas(i);
    for j=1:4
        ga=j;
        rho=morsefreq(ga,mu)./frho;
        fs{i,j}=morsespace(ga,be,4,0.5/N)';
        psi=morsewave(N,ga,mu,frho).*rho;
        w{i,j}=wavetrans(psi,{ga,be,fs{i,j}});
    end
end

for i=1:5
    be=betas(i);
    for j=1:4
        ga=j;
        subplot(4,5,i+(j-1)*5)
        fact=frac(morsefreq(ga,be),morsefreq(ga,mu)).*frac(mu+1,be).^(1./ga);
        for k=1:length(alpha)
            [t,f]=morseregion(A(k),ga,be,frho*fact);
            h=plot(t,f,'color',[1 1 1]*0.5); hold on
            [t,f]=morseregion(alpha(k),ga,be,mu,frho);    
            h=plot(t,f,'k:');hold on
        end
        xlim([-245 245]),ylim([0.8/1000 3]),hold on,ylog,
        %plot(0,frho,'.','color',[1 1 1]*0.5)
        plot(0,frho*fact,'wo','markerfacecolor','k')
        xtick([-200:100:200])
        ytick([1/1000 1/100 1/10 1])
        t=1:1000;t=t-mean(t);
        %set(gca,'TickLength',[0.03 0.03])
        set(gca,'xticklabel',[])
        set(gca,'yticklabel',[])
        contour(t,fs{i,j},abs(w{i,j})'./maxmax(abs(w{i,j})),alpha,'k')
        if i==1
            text(-200,1.5,['$\gamma=$' int2str(ga)])
        end
        if j==1
            text(80,1.5,['$\beta=$' num2str(be)])
        end
        if (i==3)&&(j==1)
            title('Generalized Morse Wavelet Regions of Influence vs. Localization Regions')
            set(gca,'TitleFontWeight','Normal')
        end
    end
end
h=packfig(4,5,'both');
orient landscape
set(gcf,'paperPosition',[.25 .25 10.5 6])
fontsize 12 10 10 10

if strcmpi(str,'print')
    print -dpng impulses-regionsofinfluence.png
    %print -djpeg impulses-regionsofinfluence.jpg
    crop('impulses-regionsofinfluence.png')
    %print -djpeg -r1000 impulses-regionsofinfluence.jpg
    %print -djpeg impulses-regionsofinfluence.jpg
    %crop('impulses-regionsofinfluence.jpg')
end
%\*************************************************************************


%/*************************************************************************
%Application to the Labrador Sea
load labseatpjaos%,load labseacensus
use labseatpjaos
figure,axis([ -64.9 -42 52 64]),hold on,topoplot continents
h=plot(lon,lat,'k');linestyle(h,'D')
h=plot(lon(:,31),lat(:,31));linestyle(h,'3D')
plot(-52.5,56.75,'wo','markerfacecolor',0*[1 1 1],'markersize',8),latratio(58.5)
jj=31;a=57;b=223;
plot(lon([a b],jj),lat([a b],jj),'wo'),
[topo,topolat,topolon]=readtopo(axis);
contour(topolon,topolat,-topo,[0:5],'k')
%h=plot(labseacensus.interiorcontour);linestyle -h h 2k--
orient portrait
title('Altimeter Tracks in the Labrador Sea')

orient portrait
fontsize 20 16 16 16

if strcmpi(str,'print')
    %print -djpeg -r1000 labseatracks.jpg
    print -dpng  labseatracks.png
    crop('labseatracks.png')
end

load labseatpjaos,use labseatpjaos
jj=31;a=57;b=223;
vindex(lat,lon,num,ssh,mss,atd,depth,dfc,jj,2);
ssh=squeeze(ssh);
badfrac=sum(0+~isfinite(ssh),2)./size(ssh,2);
wasfilled=~isfinite(ssh);
ssh=vswap(fillbad(ssh,nan,inf),nan,0);%Fill gaps and swap remaining nans
%figure,plot(badfrac,'.'),vlines([a b])  
%These are the first datapoints for which not all values are missing
vindex(lat,lon,ssh,mss,atd,depth,dfc,wasfilled,a:b,1);
ssh=ssh-vrep(vmean(ssh,1),size(ssh,1),1);  %Remove mean

%timeindexkk=find(yearfrac(num)>=2007&yearfrac(num)<2008);
%length(find(~wasfilled(:,timeindexkk))) 
%   5216 valid data points

d=spheredist(labseatpjaos.lat(:,jj),labseatpjaos.lon(:,jj),56.75,-52.5); 
[mind,bravoindex]=min(d);%12 km
bravoatd=labseatpjaos.atd(bravoindex,jj);

%min(find(depth>2)) %33
%max(find(depth>2)) %153

be=1;ga=2;mu=0;
fs=morsespace(ga,be,{0.1,pi},{2,size(ssh,1)},8);
tic;w=wavetrans(ssh,{ga,be,fs},'mirror','parallel');toc
sizew=size(w);

%Locating transform maxima
[index,ww,ff,rr]=transmax(fs,w,{ga,be,wasfilled});
[ii,jj,kk]=ind2sub(size(w),index);
lat=lat(ii);lon=lon(ii);num=num(kk);


%Estimating noise standard deviation
sigma=sqrt(vmean(vcolon(squared(w(:,1,:))),1));
wwtilde=ww./(sigma.*sqrt(ff./fs(1)));
%sigeps=sigma(1)*sqrt(frac(morsefreq(ga,be),fs(1)*morseffun(ga,be,0)));

%Forming false detection rates
N=length(find(~wasfilled(:)));
tic;[count,bins,rate]=transmaxdist(ga,be,0,fs,fs(1)./fs(2),N,2000,'extrap');toc
fdr=interp2(fs,bins,rate,ff,abs(wwtilde));
fdr(isnan(fdr))=0;

[C,rho,frho]=maxprops(ww,ff,ga,be,mu);
[A,R,Ro]=max2eddy(5.72,lat,C,rho);

clear bravotrackcensus
matsave bravotrackcensus sizew index ii jj ff kk ww rr lat lon num fdr C rho frho A R Ro sigma wasfilled

load  bravotrackcensus
use bravotrackcensus
bool=(fdr<1/1000&rr<0.1);
vindex(index,ii,jj,ff,kk,ww,rr,lat,lon,num,fdr,C,rho,frho,A,R,Ro,sigma,bool,1);
[bool,z]=isomax(sizew,index,ww,ff,ga,be,mu,1/2);
vindex(index,ii,jj,ff,kk,ww,rr,lat,lon,num,fdr,C,rho,frho,A,R,Ro,sigma,bool,1);

psihat=morsewave(size(ssh,1),ga,mu,frho);
parfor i=1:length(ii)
    psihat(:,i)=real(C(i).*vshift(psihat(:,i),round(size(ssh,1)/2-ii(i)),1)).*rho(i);
end

ssh(wasfilled)=nan;
sshhat=zeros(size(ssh));
for k=1:size(ssh,2)
    index=find(kk==k);
    if ~isempty(index)
        sshhat(:,k)=sum(psihat(:,index),2);
    end
end
sshres=ssh-sshhat;

offsetter=15;
timeindex=find(yearfrac(labseatpjaos.num)>=2007&yearfrac(labseatpjaos.num)<2008);
timeindexkk=find(yearfrac(num)>=2007&yearfrac(num)<2008);
offsetmatrix=vrep(offsetter*[0:length(timeindex)-1],size(ssh,1),1);
atdaxis=labseatpjaos.atd(a:b,31)-bravoatd;

figure
%subplot(1,3,1),plot(atdaxis,ssh(:,timeindex)+offsetmatrix,'k');linestyle 2E k
subplot(1,3,1),plot(atdaxis,ssh(:,timeindex)+offsetmatrix);
subplot(1,3,2),h=plot(atdaxis,ssh(:,timeindex)+offsetmatrix);linestyle(h,'C');
%hold on, h=plot(atdaxis,sshhat(:,timeindex)+offsetmatrix,'k');
hold on, h=plot(atdaxis,sshhat(:,timeindex)+offsetmatrix);linestyle(h,'thick')
%linestyle(h(1:2:end),'3w'),linestyle(h(2:2:end),'2w'),
%h=plot(atdaxis,sshhat(:,timeindex)+offsetmatrix,'k');
%linestyle(h(1:2:end),'2E'),linestyle(h(2:2:end),'k'),
%subplot(1,3,3),h=plot(atdaxis,sshres(:,timeindex)+offsetmatrix,'k');linestyle 2E k
subplot(1,3,3),h=plot(atdaxis,sshres(:,timeindex)+offsetmatrix);
for i=1:3,subplot(1,3,i),axis([min(atdaxis) max(atdaxis) -20 560]),
xlabel('Along-Track Distance (km)'),ytick([0:50:600]),end
subplot(1,3,2),
kko=min(find(yearfrac(labseatpjaos.num)>=2007));
for i=1:length(timeindexkk)
    plot(atdaxis(ii(timeindexkk(i))),sshhat(ii(timeindexkk(i)),kk(timeindexkk(i)))+...
       offsetter*(kk(timeindexkk(i))-kko),'wo','markerfacecolor','k')
end
subplot(1,3,1),title('Labrador Sea SSH Anomalies in 2007')
ylabel('Sea Surface Height Anomaly (cm)')
subplot(1,3,2),title('Wavelet Element Reconstruction')
subplot(1,3,3),title('Residual, SSH-Reconstruction')
subplot(1,3,1),text(-340,550,'(a)')
subplot(1,3,2),text(-340,550,'(b)')
subplot(1,3,3),text(-340,550,'(c)')
packfig(1,3,'columns')


orient portrait
fontsize 16 14 14 14
set(gcf,'paperposition',[1 1 16 10])

if strcmpi(str,'print')
    %print -djpeg -r1000 impulses-labseathreepanel-bw.png
    %print -djpeg impulses-labseathreepanel.jpg
    %crop('impulses-labseathreepanel-bw.png')
    print -djpeg -r1000 impulses-labseathreepanel.jpg
    %print -djpeg impulses-labseathreepanel.jpg
    crop('impulses-labseathreepanel.jpg')
end
%--------------------------------------------------------------------------
%\*************************************************************************




