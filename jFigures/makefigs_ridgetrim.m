function[]=makefigs_trim
%MAKEFIGS_RIDGETRIM  Makes a sample figure for RIDGETRIM.
 
t=(1:1000)';t=t-mean(t);
t=t/10;dt=t(2)-t(1);
z=rot(2*pi*t./5);
z(1:100)=0;z(end+1-100:end)=0;

fs=2*pi./(logspace(log10(10),log10(100),50)');

beta=3;gamma=3;P=sqrt(beta*gamma);

%Compute wavelet transforms using generalized Morse wavelets
[wx,wy]=wavetrans(real(z),imag(z),{gamma,beta,fs,'bandpass'},'periodic');
[wxr,wyr,ir,jr,fr,er]=ridgewalk(dt,wx,wy,fs,P);

%Time along ridges
tr=nan*ir;tr(~isnan(ir))=t(ir(~isnan(ir)));

%Trim at the alpha=1 level
[fr2,tr2,er2]=ridgetrim(dt,sqrt(beta*gamma),1,fr,tr,er);

figure,
subplot(3,1,1),uvplot(t,z),vlines(t([100 900]),'k:')
title('Example of Ridge Trimming with ALPHA=1')

subplot(3,1,2),
plot(tr,er,'r'), hold on,plot(tr2,er2,'k','linewidth',3)
ylim([10^-5 100]),ylog,vlines(t([100 900]),'k:')
ylabel('Error Estimate')

subplot(3,1,3),
contourf(t,2*pi./fs*dt,sqrt(squared(wx)+squared(wy))',100),
nocontours,colormap lansey,hold on
ylim([15 max(2*pi./fs)]*dt),ylog,flipy
plot(tr,2*pi./fr,'w','linewidth',2)
plot(tr2,2*pi./fr2,'k','linewidth',3)
hlines(50,'k:'),vlines(t([100 900]),'k:')
set(gca,'ytick',[20 30 40 50 60 70 80]/10),ylabel('Period')

packfig(3,1,'rows')

%Example of ridge trimming for a truncated complex exponential signal.
%White is before trimming and includes, black is after trimming. 
%Some spurious ridges are completely eliminated as they are too short. 
%Alpha=1 corresponds almost exactly to the start and end of the oscilation.

