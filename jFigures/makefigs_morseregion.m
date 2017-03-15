function[]=makefigs_morseregion
%MAKEFIGS_MORSEREGION  Makes somes sample figures for MORSEREGION.

%Recreating the Morse wavelet concentrations regions presented in the 
%second and fourth columns of Fig. 1 from Olhede and Walden (2003a). 
%
%Note the definition of area used in jLab differs from Olhede and Walden
%by a factor of 1/2. 

ga=[1 3 1 3];
be=[1 20 1 20];
A=[10 10 150 150]'/2;
fm=morsefreq(ga,be);
clear ax
axlim=[5 4 27 14];
figure
for i=1:4
    subplot(2,2,i)
    [t,f]=morseregion(A(i),ga(i),be(i),fm(i));
    plot(t,f),hold on,plot(t,-f)
    axis([-1 1 -1 1]*axlim(i)),axis square
    title(['$\gamma$ = ' int2str(ga(i)) ', $\beta$ = ' int2str(be(i)), ', A = ' num2str(A(i))])
end
 
%Localization regions for GAMMA=[1 3 9 27] and BETA=[1 3 9 27].
%Different values of GAMMA are nested, while different values of BETA are 
%offset horizontally.  The peak frequency is shown with a horizontal line.

figure
ga1=[1 3 9 27];
be1=[1 3 9 27];
[ga,be]=meshgrid(ga1,be1);
vcolon(ga,be);
[t,f]=morseregion(10,ga,be,1);
dt=200*oprod(ones(size(t,1),1),sqrt(sqrt(be)));
h=plot(t+dt,f,'k');

linestyle(h(ga==3),'2k-')
linestyle(h(ga==9),'2k--')
linestyle(h(ga==27),'k--')
title('Morse wavelet localization regions')
xlabel('Normalized Time')
ylabel('Normalized Frequency'),
hlines(1,':')

