function [fm,fe,fi,cf] = morsefreq(ga,be)
%MORSEFREQ  Frequency measures for generalized Morse wavelets. [with F. Rekibi]
%
%   [FM,FE,FI]=MORSEFREQ(GAMMA,BETA) calculates three different measures of
%   the frequency of the lowest-order generalized Morse wavelet specified 
%   by parameters GAMMA and BETA.
%
%   FM is the modal or peak, FE is the "energy" frequency, and FI is the 
%   instantaneous frequency at the wavelet center.
%
%   [FM,FE,FI,CF]=MORSEFREQ(GAMMA,BETA) also computes the curvature CF of 
%   the instantaneous frequency at the wavelet center. 
%
%   Note that all frequency quantities here are *radian* as in cos(omega t)
%   and not cyclic as in cos(2 pi f t).
% 
%   The input parameters must either be matrices of the same size,
%   or some may be matrices and the others scalars.   
%
%   For BETA=0, the "wavelet" becomes an analytic lowpass filter, and FM 
%   is not defined in the usual way.  Instead, FM is defined as the point
%   at which the filter has decayed to one-half of its peak power. 
%
%   For details see
%
%       Lilly and Olhede (2009).  Higher-order properties of analytic 
%           wavelets.  IEEE Trans. Sig. Proc., 57 (1), 146--160.
%
%   See also MORSEBOX, MORSEPROPS, MORSEWAVE.
%
%   Usage: fm = morsefreq(ga,be);
%          [fm,fe,fi] = morsefreq(ga,be);  
%          [fm,fe,fi,cf] = morsefreq(ga,be);  
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2023 J. M. Lilly and F. Rekibi
%                         --- type 'help jlab_license' for details    

%   'morsefreq --f' generates a sample figure, but I'm hiding this since
%   it's uncommented

if strcmpi(ga,'--f')
  type makefigs_morsefreq
  makefigs_morsefreq;
  return
end

%fm=frac(be,ga).^frac(1,ga)
fm=exp(frac(1,ga).*(log(be)-log(ga)));
fm(be==0)=(log(2)).^frac(1,ga(be==0)); %Half-power point

%e^-omega^gamma

%fm(be==0)=sqrt(3)*sqrt(frac(gamma(frac(3,ga(be==0))),gamma(frac(1,ga(be==0)))));
%Instead the square root of the second
%%   frequency-domain moment is used.  For the Gaussian or GAMMA=2 case,   
%   this gives MORSEFREQ(2,0)=1/SQRT(2)
%fm(be==0)=(log(2)).^frac(1,ga(be==0)); %Half-power point
%length(find(be==0))
%fm(be==0)=sqrt(frac(gamma(frac(3,ga(be==0))),gamma(frac(1,ga(be==0)))));
%sqrt(frac(gamma(frac(3,ga1)),gamma(frac(1,ga1))))
%Oher possibilities for zero beta case
%fm(be==0)=frac(be(be==0)+1,ga(be==0)).^frac(1,ga(be==0));%Peak of next wavelet
%    fm=frac(ga-1,ga).^frac(1,ga);%Most rapidly decreasing point ---this moves to higher frequencies!
%    fm=exp(frac(1,ga).*(log(be)-log(ga)));
%    fm(be==0)=sqrt(frac(gamma(frac(3,ga(be==0))),gamma(frac(1,ga(be==0)))));%Second moment

%Instead the frequency-domain half-
%   power point is used, reflecting the bawidth of the filter.  
% 
if nargout>1
    fe=frac(1,2.^frac(1,ga)).*frac(gamma(frac(2*be+2,ga)),gamma(frac(2*be+1,ga)));
end

if nargout>2
    fi=frac(gamma(frac(be+2,ga)),gamma(frac(be+1,ga)));
end

if nargout>2
    [m2,n2,k2]=morsemom(2,ga,be);
    [m3,n3,k3]=morsemom(3,ga,be);
    cf=-frac(k3,sqrt(k2.^3));
end

if 0
%    C = morsecfun(A,ga,be);
%   MORSEFREQ uses the formula of Olhede and Walden (2002),
%   "Generalized Morse Wavelets", the first formula on page 2665.
    r=(2*be+1)./ga;
    coeff=gamma(r+(1./ga)).*2.^(-1./ga);
    fmin=coeff./(2*pi.*gamma(r).*(C+sqrt(C.^2-1)).^(1./ga));
    fmax=coeff./(2*pi.*gamma(r).*(C-sqrt(C.^2-1)).^(1./ga));

    if nargout >4
        fo=frac(1,2).*(fmax+fmin);
    end
    if nargout >5
        fw=frac(1,2).*(fmax-fmin);
    end
end
