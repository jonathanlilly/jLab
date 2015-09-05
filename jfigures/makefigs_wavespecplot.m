function[]=makefigs_wavespecplot
%MAKEFIGS_WAVESPECPLOT  Makes a sample figure for WAVESPECPLOT.


load npg2006 
use npg2006
cv=npg2006.cv; 

fs=morsespace(2,4,2*pi/10,2*pi/500,8);

vindex(num,lat,lon,p,t,cv,cx,1:size(cx,1)-3,1); %Strip off a few NaNs at the end

%Compute wavelet transforms using generalized Morse wavelets
wx=wavetrans(real(cx),{1,2,4,fs,'bandpass'},'mirror');
wy=wavetrans(imag(cx),{1,2,4,fs,'bandpass'},'mirror');

[wp,wn]=vectmult(tmat,wx,wy);

h=wavespecplot(num,cv,dt./fs,wp,wn,0.5);
axes(h(1)),title('Positive and Negative Wavelet Transforms of Lilly and Gascard (2006) Time Series')
axes(h(3)),cax=caxis;
axes(h(2)),caxis(cax);
