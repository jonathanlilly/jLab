function[]=makefigs_ridgewalk
%MAKEFIGS_RIDGEWALK  Makes a sample figure for RIDGEWALK.
 
load ebasnfloats
use ebasnfloats

figure
num=num{33};lat=lat{33};lon=lon{33};
dt=num(2)-num(1);

ga=3;be=3;P=sqrt(ga*be);
fs=morsespace(ga,be,{0.05,2*pi/3},2*pi/100,8);

%Compute wavelet transforms using generalized Morse wavelets
cx=fillbad(latlon2xy(lat,lon,30,-25));
cv=latlon2uv(num,lat,lon);

[wx,wy]=wavetrans(real(cx),imag(cx),{ga,be,fs,'bandpass'},'mirror');
[wrx,wry,ir,jr,fr,er]=ridgewalk(dt,wx,wy,fs,P,4./(2*P./pi));  
fr=ridgemap(length(cx),fr,ir);

ci=(0:5:65);
numo=datenum(1986,1,1)-1;
[h,hl]=wavespecplot(num-numo,cv,2*pi./fs,sqrt(abs(wx).^2+abs(wy).^2),1,ci);
linestyle -h hl k k--
axes(h(1)),ylim([-18 18]),ylabel('Current Speed (cm/s)'),title('Bivariate Ridge Method Example')
text(-90,15,'(a)')

axes(h(2)),caxis([0 40]),colormap gray,flipmap,ylim([3.6 60]),hold on

plot(num-numo,2*pi./fr,'w','linewidth',4)
plot(num-numo,2*pi./fr,'k','linewidth',2)

xlabel('Day of Year 1986'),ylabel('Period in Days')
set(gca,'ytick',2.^(2:.5:5.5))
set(gca,'yticklabel',[' 4';'  ';' 8';'  ';'16';'  ';'32';'  '])
inticks
text(-90,4.5,'(b)')

fontsize 12 10 10 10
set(gcf,'paperposition',[1 1 9 5])

%For printing
if 0
  currentdir=pwd;
  cd([whichdir('jlab_license') '/figures'])
  print -dpng ridgewalk_example
  crop ridgewalk_example.png
  cd(currentdir)
end