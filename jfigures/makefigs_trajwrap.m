function[]=makefigs_trajwrap
%MAKEFIGS_TRAJWRAP  Makes a sample figure for TRAJWRAP.

load vortex
use vortex.drifters
cx=x+1i*y;
cx_unwrap=trajunwrap(cx,[L*pi,nan]);
cx(abs(vdiff(real(cx),1))>100)=nan;

figure,
subplot(2,1,1),plot(cx_unwrap),axis equal, xlim([-2000  -2000+4*pi*L]),ylim([-1000 1000])
vlines([-1 1]*L*pi,'k')
title('Unwrapped float trajectories')
subplot(2,1,2),plot(cx),axis equal,  xlim([-2000  -2000+4*pi*L]),ylim([-1000 1000])
vlines([-1 1]*L*pi,'k')
title('Wrapped float trajectories')

