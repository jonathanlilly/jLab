function[]=makefigs_dlines
%MAKEFIGS_DLINES  Makes a sample figure for DLINES.

figure
subplot(2,2,1)
axis([-1 0 0 1]),dlines(-1,'b'),dlines(1,'r')
subplot(2,2,2),
axis([0 1 0 1]),dlines(-1,'b'),dlines(1,'r')
subplot(2,2,3)
axis([-1 0 -1 0]),dlines(-1,'b'),dlines(1,'r')
subplot(2,2,4)
axis([0 1 -1 0]),dlines(-1,'b'),dlines(1,'r')
