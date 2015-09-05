function[]=makefigs_morsearea
%MAKEFIGS_MORSEAREA  Makes a sample figure for MORSEAREA.

C1=logspace(0.1,3,100)';
be1=linspace(1,10,101)';
[be,C]=meshgrid(be1,C1);
ga=2;

figure
A = morsearea(C,ga(1),be);
contourf(C1,be1,log10(A'),20)
xlog
xlabel('Parameter C')
ylabel('Parameter \beta')
title('Log10 Area with \gamma =2')
colorbar
