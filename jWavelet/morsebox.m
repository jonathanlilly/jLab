function[a,sigt,sigo,skew]=morsebox(ga,be)
%MORSEBOX  Heisenberg time-frequency box for generalized Morse wavelets.
%
%   A=MORSEBOX(GAMMA,BETA) returns the time-bandwidth product of the
%   a generalized Morse wavelet, that is, the area of its Heisenberg box.
%
%   [A,SIGT,SIGO,SKEW]=MORSEBOX(GAMMA,BETA) also returns the box width in 
%   time SIGT and in frequency SIGO, with A=SIGT*SIGO, as well as the 
%   frequency-domain skewness SKEW.
%
%   Both SIGT and SIGO are non-dimensionalized with respect to the (radian)
%   peak frequency MORSEFREQ, as in Lilly and Olhede (2009) given below.
%
%   Not for all values of BETA and GAMMA are these standard derivations
%   defined as real-valued quantities.  For locations where these would 
%   be imaginary, their values are set to NAN.
%
%   For details see
%
%       Lilly and Olhede (2009).  Higher-order properties of analytic 
%           wavelets.  IEEE Trans. Sig. Proc., 57 (1), 146--160.
%
%   In this calculation, care is taken to avoid potential problems caused
%   by large arguments to the gamma function, by utilizing GAMMALN.  This
%   lets the Heisenberg area be computed for small values of GAMMA/BETA.  
%
%   See also MORSEFREQ, MORSEPROPS, MORSEWAVE.
%
%   'morsebox --t' runs a test.
%   'morsebox --f' generates a sample figure.
%
%   Usage: a=morsebox(ga,be);
%          [a,sigt,sigo]=morsebox(ga,be);
%          [a,sigt,sigo,skew]=morsebox(ga,be);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2014 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(ga, '--t')
    morsebox_test,return
end
if strcmpi(ga, '--f')
    type makefigs_morsebox
    makefigs_morsebox;
    return
end

be(be==1)=1+1e-10;  %Correction for unstable ratios and beta=1

% [mo,no]=morsemom(0,ga,be);
% [m1,n1]=morsemom(1,ga,be);
% [m2,n2]=morsemom(2,ga,be);
% [m3,n3]=morsemom(3,ga,be);
% 
% [ma,na]=morsemom(0,ga,be-1);
% [mb,nb]=morsemom(0,ga,be-1+ga);
% [mc,nc]=morsemom(0,ga,be-1+ga/2);
% 
% rata=frac(morseafun(ga,be),morseafun(ga,be-1)).^2;
% ratb=frac(morseafun(ga,be),morseafun(ga,be-1+ga)).^2;
% ratc=frac(morseafun(ga,be),morseafun(ga,be-1+ga/2)).^2;
%
% sig2a=  rata.*be.^2.*frac(na,no);
% sig2b=  ratb.*ga.^2.*frac(nb,no);
% sig2c=2*ratc.*be.*ga.*frac(nc,no);
% sigt=real(om.*sqrt(sig2a+sig2b-sig2c));  %Earlier way

om = morsefreq(ga,be);
be(be<1/2)=nan;

logsigo1=frac(2,ga).*log(frac(ga,2*be))+gammaln(frac(2*be+1+2,ga))-gammaln(frac(2*be+1,ga));
logsigo2=frac(2,ga).*log(frac(ga,2*be))+2.*gammaln(frac(2*be+2,ga))-2.*gammaln(frac(2*be+1,ga));

sigo=sqrt(exp(logsigo1)-exp(logsigo2));  %This works, yay,
%sigo=frac(1,om).*sqrt(frac(n2,no)-frac(n1,no).^2);  %Earlier way

ra=2*morse_loga(ga,be)-2*morse_loga(ga,be-1)+morse_loga(ga,2*(be-1))-morse_loga(ga,2*be);
rb=2*morse_loga(ga,be)-2*morse_loga(ga,be-1+ga)+morse_loga(ga,2*(be-1+ga))-morse_loga(ga,2*be);
rc=2*morse_loga(ga,be)-2*morse_loga(ga,be-1+ga./2)+morse_loga(ga,2*(be-1+ga./2))-morse_loga(ga,2*be);

logsig2a=ra+frac(2,ga).*log(frac(be,ga))+2*log(be)+gammaln(frac(2*(be-1)+1,ga))-gammaln(frac(2*be+1,ga));
logsig2b=rb+frac(2,ga).*log(frac(be,ga))+2*log(ga)+gammaln(frac(2*(be-1+ga)+1,ga))-gammaln(frac(2*be+1,ga));
logsig2c=rc+frac(2,ga).*log(frac(be,ga))+log(2)+log(be)+log(ga)+gammaln(frac(2*(be-1+ga./2)+1,ga))-gammaln(frac(2*be+1,ga));

sig2a=exp(logsig2a);
sig2b=exp(logsig2b);
sig2c=exp(logsig2c);
sigt=sqrt(sig2a+sig2b-sig2c);

sigt(be<=1/2)=nan;
sigt=real(sigt);
sigo=real(sigo);

a=sigt.*sigo;

skew=morseskew(ga,be);

function[loga]=morse_loga(ga,be)
loga=frac(be,ga).*(1+log(ga)-log(be));

function[logmom]=morse_logmom(p,ga,be)
logmom=morse_loga(ga,be)-log(2*pi*ga)+gammaln(frac(be+1+p,ga));

function[skew]=morseskew(ga,be)
%Calculate frequency-domain skewness using GAMMALN
mo=morse_logmom(0,ga,be);
m1=morse_logmom(1,ga,be);
m2=morse_logmom(2,ga,be);
m3=morse_logmom(3,ga,be);
logskewa=m3-mo;
logskewb=log(3)+m2+m1-2*mo;
logskewc=log(2)+3*m1-3*mo;
logskewd=m2-mo;
logskewe=2*m1-2*mo;

skew=frac(exp(logskewa)-exp(logskewb)+exp(logskewc),(exp(logskewd)-exp(logskewe)).^(3/2));



function[]=morsebox_test
 
t=(-500:500)';
n=0;
ga1=(2:2:10);
be1=(2:2:10);
clear sigmat sigmao
for i=1:length(ga1)
    for j=1:length(be1)
        n=n+1;
        fs=2*pi/20;
        psi=morsewave(length(t),1,ga1(i),be1(j),fs,'bandpass');
        [mut,sigmat(j,i)]=pdfprops(t,abs(psi).^2); 
        sigmat(j,i)=sigmat(j,i).*fs;
    
        fs=2*pi/5;
        psi=morsewave(length(t),1,ga1(i),be1(j),fs,'bandpass');
        [muo,sigmao(j,i)]=pdfprops(2*pi*fftshift(t./length(t)),abs(fft(psi)).^2); 
        sigmao(j,i)=sigmao(j,i)./(fs);
    end
end


[ga,be]=meshgrid(ga1,be1);
[a,sigt,sigo]=morsebox(ga,be);

reporttest('MORSEBOX numerical trials for SIGT',aresame(sigt,sigmat,1e-5))
reporttest('MORSEBOX numerical trials for SIGO',aresame(sigo,sigmao,1e-3))



