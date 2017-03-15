function[a]=morseafun(varargin)
%MORSEAFUN  Returns the generalized Morse wavelet amplitude or a-function.
%
%   MORSEAFUN is a low-level function called by many a number of the Morse
%   wavelet functions.
%
%   A=MORSEAFUN(GAMMA,BETA) returns the generalized Morse wavelet 
%   amplitude, called "A_{BETA,GAMMA}" by Lilly and Olhede (2009).
%
%   By default, A is chosen such that the maximum of the frequency-
%   domain wavelet is equal to 2, the ``bandpass normalization.''
%
%   A=MORSEAFUN(GAMMA,BETA,'energy') instead returns the coefficient
%   giving the wavelet unit energy.  
%
%   A=MORSEAFUN(K,GAMMA,BETA,'energy') returns the unit energy coefficient 
%   appropriate for the Kth-order wavelet.  The default choice is K=1.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2016 J.M. Lilly --- type 'help jlab_license' for details

if strcmpi(varargin{1},'--t')
      morseafun_test;return
end

str='band';
if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
end

k=1;
if length(varargin)==3
    k=varargin{1};
    varargin=varargin(2:end);
end

ga=varargin{1};
be=varargin{2};
if strcmpi(str(1:3),'ban')
    om=morsefreq(ga,be);     
%    a=frac(2,(om.^be).*exp(-om.^ga));
    a=frac(2,exp(be.*log(om)-om.^ga));
    a(be==0)=2;
elseif strcmpi(str(1:3),'tes')
    a=sqrt(frac(2 *pi*ga.*2.^frac(2*be+1,ga),gamma(frac(2*be+1,ga)))); 
elseif strcmpi(str(1:3),'ene')
    r=frac(2*be+1,ga);
    a=double((2*pi*ga.*(2.^r).*exp(gammaln(k)-gammaln(k+r-1))).^(1/2));
end

function[]=morseafun_test

ga1=(2:1:9);
be1=(1:1:10);
[ga,be]=meshgrid(ga1,be1);
om=reshape(morsefreq(ga,be),10,8);

dom=0.01;
omgrid=permute((0:dom:20)',[3 2 1]);
omgrid=vrep(omgrid,length(ga1),2);
omgrid=vrep(omgrid,length(be1),1);

omgrid=omgrid.*vrep(om,size(omgrid,3),3);
a=morseafun(ga,be,'energy');

begrid=vrep(be,size(omgrid,3),3);
gagrid=vrep(ga,size(omgrid,3),3);
agrid=vrep(a,size(omgrid,3),3);

psi=agrid.*omgrid.^begrid.*exp(-omgrid.^gagrid);
psiint=vsum(psi.^2,3).*dom.*om./(2*pi);

reporttest('MORSEAFUN unit energy', allall(abs(psiint-1)<1e-2))

a2=morseafun(ga,be,'test');
reporttest('MORSEAFUN unit energy, alternate formulation', aresame(a,a2,1e-6));


