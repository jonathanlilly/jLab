function[]=makefigs_morsebox
%MAKEFIGS_MORSEBOX  Makes a sample figure for MORSEBOX.

ga1=(1/3:.1:11);
be1=(1/3:.1:10);

[ga,be]=meshgrid(ga1,be1);
[fm,fe,fi,cf] = morsefreq(ga,be);
a=morsebox(ga,be);

figure
contourf(ga1,be1,a,(.5:.01:.6))
hold on, contour(ga1,be1,a,(.500:.002:.51),'k:')
colormap gray,flipmap,
axis([1 10 1 10])
xtick(1:10),ytick(1:10)
ax=gca;
hc=colorbar;

contour(ga1,be1,(fm-fe)./(2*pi),[ 0 0],'k','linewidth',2)
contour(ga1,be1,(fm-fi)./(2*pi),[ 0 0],'k','linewidth',2)
contour(ga1,be1,cf./(2*pi),[ 0 0],'k','linewidth',2)
caxis([.5 .6])
vlines(3,'k--')
plot(ga1,12./ga1,'k')
 
title('Morse Wavelet Area and Transitions')
xlabel('Gamma Parameter')
ylabel('Beta Parameter')
plot([1+sqrt(-1)*0 10+sqrt(-1)*9/2],'k','linewidth',3)
plot([1+sqrt(-1)*0 10+sqrt(-1)*9/2],'w--','linewidth',2)



