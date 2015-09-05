function[]=makefigs_morsefreq()
%MAKEFIGS_MORSEFREQ  Makes a sample figure for MORSEFREQ.

ga=[2 3 4 8];
be=[(1:.1:10) (10:.5:100)];

[ga,be]=meshgrid(ga,be);
be(be<frac(ga-1,2))=nan;

[fm,fe,fi,cf]=morsefreq(ga,be);
p=frac(sqrt(be.*ga),pi);

figure
subplot(221)
plot(p,fe./fm),linestyle k 2k k-- k-. 
title('Energy Frequency / Peak Frequency')
ylim([0.90 1.15]),xlim([0.25 3/2]*2),
hlines(1,'k:'),xtick((0:10)/2);ytick(0.8:0.05:1.2),
xlabel('Duration / Period Ratio'),ylabel('Frequency Ratio')
fixlabels([-1 -2])

subplot(222)
plot(p,fi./fm),linestyle k 2k k-- k-. 
title('Instantaneous / Peak Frequency')
ylim([0.90 1.15]),xlim([0.25 3/2]*2),
hlines(1,'k:'),xtick((0:10)/2);ytick(0.8:0.05:1.2),
xlabel('Duration / Period Ratio'),ylabel('Frequency Ratio')
fixlabels([-1 -2])

subplot(223)
plot(p,cf/2/pi),linestyle k 2k k-- k-. 
title('Instantaneous Frequency Curvature')
ylim([-.15 .15]),xlim([0.25 3/2]*2),
hlines(0,'k:'),xtick((0:10)/2);
ytick(-.15:.05:.15)
xlabel('Duration / Period Ratio'),ylabel('Dimensionless Curvature')
fixlabels([-1 -2])

[a,dt,dom]=morsebox(ga,be);
subplot(224)
plot(p,a),linestyle k 2k k-- k-. 
title('Heisenberg Box Area')
ylim([.5 .55]),xlim([0.25 3/2]*2),
hlines(0,'k:'),xtick((0:10)/2);
ytick(.5:.01:.6)
xlabel('Duration / Period Ratio'),ylabel('Area')
fixlabels([-1 -2])


%/********************************************************************
%The rest is to compare with the numerical computation of wavelet properties
%m=(0.3:.05:2);
be=1:.1:6;
f=2*pi*(0.33:.05:2)./2;

fs=2*pi/100;

t=(-10000:10000)';

clear fm fe fi cf p a 
for i=1:length(f)    
    psi=morlwave(length(t),f(i),fs);
    [fm(i),fe(i),fi(i),cf(i),p(i),a(i)]=morsefreq_numerical(t,psi,fs);
end
%h=gcf;figure,plot(m,fm./fs,'+'),figure(h)

clear fm2 fe2 fi2 cf2 p2 a2 
for i=1:length(be)    
    psi=morsewave(length(t),1,3,be(i),fs,'bandpass');
    [fm2(i),fe2(i),fi2(i),cf2(i),p2(i),a2(i)]=morsefreq_numerical(t,psi,fs);
end


letterlabels(2)

subplot(221)
plot(p2,fe2./fm2,'k.'),plot(p,fe./fm,'k:'),plot(p,fe./fm,'k+')
subplot(222)
plot(p2,fi2./fm2,'k.'),plot(p,fi./fm,'k:'),plot(p,fi./fm,'k+')
subplot(223)
plot(p2,cf2,'k.'),plot(p,real(cf),'k:'),plot(p,real(cf),'k+')
subplot(224)
plot(p2,a2,'k.'),plot(p,a,'k:'),plot(p,a,'k+')


function[fm,fe,fi,cf,p,a]=morsefreq_numerical(t,psi,fs)
om_axis=vshift(2*pi*fftshift(t./length(t)),-1,1);
vswap(om_axis,0,1e-10);

W=fft(fftshift(psi(:,:,1)));
psi=psi.*2./maxmax(abs(W));

[maxpsi,index]=max(abs(fft(psi)));
om=om_axis(index);    
%om=2*pi*fs;
fm=om./(2*pi);

[mut,sigma]=pdfprops(t,psi.*rot(-t.*fm*2*pi)); 
p=real(sigma).*(fm)*2;

dpsi=vdiff(imlog(psi),1);      

omi=dpsi((length(t)+1)/2,1);
fi=omi/(2*pi);

ddom=vdiff(vdiff(dpsi,1),1);

[mut,sigmat]=pdfprops(t,abs(psi).^2); 
[ome,sigmao]=pdfprops(om_axis,abs(fft(psi)).^2);
fe=ome/(2*pi);
cf=frac(ome.^3,ome.^2).*ddom((length(t)+1)/2,1)./sigmao.^3;
a=sigmao.*sigmat;
