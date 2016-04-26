function[]=makefigs_jpcolor
%MAKEFIGS_JPCOLOR  Makes a sample figure for JPCOLOR.

figure
%mat=[1 2 3; 4 5 6; 7 8 9]';
mat=zeros(3,3)+randn(3,3)/5;
mat(2,2)=1;
subplot(1,2,1)
pcolor([1:3],[1:3],mat),axis equal,axis tight
title('A 3x3 matrix plotted by PCOLOR')
subplot(1,2,2)
jpcolor([1:3],[1:3],mat),axis equal,axis tight
title('A 3x3 matrix plotted by JPCOLOR')