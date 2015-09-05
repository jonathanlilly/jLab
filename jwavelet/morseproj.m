function[b]=morseproj(ga,be1,be2)
%MORSEPROJ  Projection coefficient for two generalized Morse wavelets.
%
%   MORSEPROJ is a low-level function called by MORSEWAVE.
%
%   B=MORSEPROJ(GA,BE1,BE2) returns the projection coefficient B for the 
%   projection of one generalized Morse wavelet onto another.
%
%   The first generalized Morse has parameters GA and BE1, and the second
%   has parameters GA and BE2.
%
%   'morseproj --t' runs a test.
%
%   Usage: b=morseproj(ga,be1,be2);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(ga, '--t')
    morseproj_test,return
end
 
c=frac(morseafun(ga,be1,'energy').*morseafun(ga,be2,'energy'),2.^frac(be1+be2+1,ga));
b=c.*frac(1,2*pi*ga).*gamma(frac(be1+be2+1,ga));

function[]=morseproj_test
 
ga1=(2:1:9);
be1=(1:1:10);
[ga,be]=meshgrid(ga1,be1);
om=morsefreq(ga,be);

dom=0.01;
omgrid=permute((0:dom:20)',[3 2 1]);
omgrid=vrep(omgrid,length(ga1),2);
omgrid=vrep(omgrid,length(be1),1);

omgrid=omgrid.*vrep(om,size(omgrid,3),3);
a=morseafun(ga,be,'energy');

begrid=vrep(be,size(omgrid,3),3);
gagrid=vrep(ga,size(omgrid,3),3);
agrid=vrep(a,size(omgrid,3),3);

psi1=agrid.*omgrid.^begrid.*exp(-omgrid.^gagrid);

psi2=morseafun(gagrid,2+0*begrid,'energy').*omgrid.^2.*exp(-omgrid.^gagrid);

b1=vsum(psi1.*psi2,3).*dom.*om./(2*pi);
b2=morseproj(ga,be,2);

reporttest('MORSEPROJ',aresame(b1,b2,1e-3))
