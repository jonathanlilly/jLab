function[]=makefigs_mspec
%MAKEFIGS_MSPEC  Makes some sample figures for MSPEC.

%First figure
load bravo94
use bravo94
use bravo94.rcm

figure
cv=cv(:,5);

clear psi f spp snn spn
psi{1}=ones(size(cv))./sqrt(length(cv));
psi{2}=sleptap(length(cv),4); 
psi{3}=sleptap(length(cv),32); 

for i=1:length(psi)
    [f(:,i),spp(:,i),snn(:,i),spn(:,i)]=mspec(cv,conj(cv),psi{i}); 
end

h=twospecplot(f,spp,snn);
axes(h(1)),vlines(abs(corfreq(lat))), linestyle D b r k:
ax=axis;axis([10^-2.95 ax(2) 10^-3 10^5]),xtick(10.^[-3 -2 -1 0])
xlabel('Frequency (rad/hour)'),ylabel('Power Spectral Density')
axes(h(2)),vlines(abs(corfreq(lat))), linestyle D b r k:
ax=axis;axis([10^-2.95 ax(2) 10^-3 10^5]),xtick(10.^[-3 -2 -1 0])
xlabel('Frequency (rad/hour)')

%To print
if 0
    set(gcf,'paperposition',[1 1 10 5.5])
    currentdir=pwd;
    cd([whichdir('jlab_license') '/figures'])
    print -dpng mspec
    crop mspec.png
    cd(currentdir)
end

%Former example
load bravo94
x=bravo94.rcm.cv;
vswap(x,nan,0);
[psi,lambda]=sleptap(length(x),16);
[f,sp,sn,spn]=mspec(x,psi);
[f,su,sv,suv]=mspec(real(x),imag(x),psi);

figure,plot(f,[sp sn]),xlog,ylog,axis tight
title('Counterclockwise (blue) and clockwise (red) spectra'),
linestyle b b b b b b r r r r r r

load bravo94
x=bravo94.rcm.cv;
vswap(x,nan,0);
[psi,lambda]=sleptap(length(x),8);

for i=1:size(x,2)
    [f,Suu,Suu3,Cuu]=mspec(real(x(:,i)),real(x(:,3)),psi);
    gammauu(:,i)=Cuu./sqrt(Suu.*Suu3);
end

figure,
plot(f,abs(gammauu)),xlog,yoffset 1,axis tight
title('Coherence of u(t) at each depth vs. u(t) at \#3','interpreter','latex')


