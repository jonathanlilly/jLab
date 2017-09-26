function[f]=morsehigh(varargin)
%MORSEHIGH  High-frequency cutoff of the generalized Morse wavelets.
%
%   MORSEHIGH is a low-level function called by MORSESPACE in order to 
%   set the high-frequency cutoff of the wavelet transform.
%
%   FHIGH=MORSEHIGH(GAMMA,BETA,ETA) returns the high-frequency cutoff
%   FHIGH of the generalized Morse wavelet specified by GAMMA and BETA, 
%   with cutoff level ETA.
%
%   Specifically, FHIGH will give the highest possible radian frequency 
%   for which the generalized Morse wavelet MORSEWAVE(N,GAMMA,BETA,FHIGH)
%   will have a greater than ETA times its maximum value at the Nyquist.
%   Here N is any choice of wavelet length.
%
%   If F=MORSEFREQ(GAMMA,BETA) is the wavelet peak frequency, then 
%
%      PSI(PI*F/FHIGH) <= ETA * PSI(F) 
%
%   is the cutoff condition. See Appendix C of Lilly (2017) for details. 
%  
%   The input parameters may either all be scalars, or GAMMA and BETA
%   may be arrays of the same size with scalar ALPHA.
%  
%   See also MORSEFREQ, MORSEWAVE.
%
%   Usage: falpha=morsehigh(ga,be,alpha);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2017 J.M. Lilly --- type 'help jlab_license' for details
 
ga=varargin{1};
be=varargin{2};
alpha=varargin{3};

N=10000;
omhigh=linspace(0,pi,N)';

f=0*ga;
for i=1:length(ga(:))
    om=morsefreq(ga(i),be(i)).*pi./omhigh;
    
    %Use logs to avoid errors for really small gammas
    %Note that ln(2)'s cancel
    lnpsi1=frac(be(i),ga(i))*log(frac(exp(1)*ga(i),be(i)));
    lnpsi2=be(i).*log(om)-om.^ga(i);
    lnpsi=lnpsi1+lnpsi2;
    index=find(log(alpha)-lnpsi<0,1,'first');
    
    %psi=frac(1,2)*morseafun(ga(i),be(i)).*(om.^be(i)).*exp(-om.^ga(i));
    %index=find(alpha-psi<0,1,'first');
    %psi(1)
    %figure,plot(lnpsi),hlines(log(alpha))
    %index
    f(i)=omhigh(index);
end

% Faster but only marginally
% f2=0*ga;
% om=pi./omhigh;
% for i=1:length(ga(:))
%     psi=(om.^be(i)).*exp(-frac(be(i),ga(i))*(om.^ga(i)-1));
%     index=find(alpha-psi<0,1,'first');
%     f2(i)=omhigh(index);
% end
% aresame(f,f2)



