function[]=makefigs_ridgewalk
%MAKEFIGS_RIDGEWALK  Makes a sample figure for RIDGEWALK.
 
load ebasnfloats
use ebasnfloats

figure
len=cellength(lat);
index=find(len>200);
lato=30;

id=id(index);num=num(index);lat=lat(index);lon=lon(index);
p=p(index);t=t(index);

index=24;
id=id(index);num=num{index};lat=lat{index};lon=lon{index};
p=p{index};t=t{index};
dt=num(2)-num(1);

ga=3;be=3;
fs=morsespace(ga,be,{0.05,2*pi/3},2*pi/100,8);

%Compute wavelet transforms using generalized Morse wavelets
cx=fillbad(latlon2xy(lat,lon,30,-25));
cv=latlon2uv(num,lat,lon);

wx=wavetrans(real(cx),{1,ga,be,fs,'bandpass'},'mirror');
wy=wavetrans(imag(cx),{1,ga,be,fs,'bandpass'},'mirror');

[ir,jr,xr,yr,fxr,fyr]=ridgewalk(dt,wx,wy,fs,{2*morseprops(ga,be),0});  
[xhat,yhat,fxhat,fyhat]=ridgemap(length(cx),xr,yr,fxr,fyr,ir);
fbar=vmean([fxhat fyhat],2,squared([xhat yhat]));

ci=(0:5:65);
numo=datenum(1986,1,1)-1;
[h,hl]=wavespecplot(num-numo,cv,2*pi./fs,sqrt(abs(wx).^2+abs(wy).^2),1,ci);
linestyle -h hl k k--
axes(h(1)),ylim([-18 18]),ylabel('Current Speed (cm/s)'),title('Bivariate Ridge Method Example')
text(-90,15,'(a)')

axes(h(2)),caxis([0 40]),colormap gray,flipmap,ylim([3.6 60]),hold on

plot(num-numo,2*pi./fbar,'w','linewidth',4)
plot(num-numo,2*pi./fbar,'k','linewidth',2)

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