function [w,W]=morlwave(N,fo,fs,str)
%MORLWAVE  Morlet wavelet.   
%  
%   The Morlet wavelet is not recommended.  Instead use a generalized Morse
%   wavelet with GAMMA set to 3.  This function is included for testing and
%   comparison purposes.  See Lilly and Olhede (2009, 2012) for details.
%
%   PSI=MORLWAVE(N,FMAX,F) returns the complex-valued Morlet wavelet of 
%   length N, with a peak in the frequency domain at *radian* F, based on a
%   mother wavelet having a peak frequency FMAX>1.
%  
%   If F is a scalar, PSI is a column vector of length M.  If F is an
%   array, PSI is a matrix of size M x LENGTH(F), with the columns in order  
%   of decreasing frequency.
%
%   The peak frequency FMAX should be greater than one.  As FMAX increases
%   from unity, the number of oscillations in the wavelet increases.
%
%   Note that the wavelets are centered at the midpoint in time, row 
%   number ROUND(SIZE(PSI,1)/2).
%
%   [PSI,PSIF]=MORLWAVE(N,FMAX,F) also returns a size SIZE(PSI) matrix 
%   PSIF which is the frequency-domain version of the wavelets.
%
%   For further details on this representation of the Morlet wavelet, see
%   Appendix A of
%
%       Lilly and Olhede (2009).  Higher-order properties of analytic 
%           wavelets.  IEEE Trans. Sig. Proc., 57 (1), 146--160.
%   _________________________________________________________________
%
%   Normalization
%
%   MORLWAVE supports two kinds of normalization for the wavelets.
%
%   MORLWAVE(..., 'bandpass') uses "bandpass normalization", the 
%   default.  See Lilly and Gascard (2006).  This implies that the 
%   FFT of wavelet has a maximum value of 2. 
%
%   MORLWAVE(..., 'energy') uses the unit energy normlization.  The 
%   time-domain wavelet energy SUM(ABS(PSI).^2) is then unity. 
%
%   Note MORLWAVE now uses bandpass-normalization by default.
%   _________________________________________________________________
%
%   See also MORLFREQ, MORLPROPS, MORSEWAVE.
%   
%   Usage: psi=morlwave(n,fmax,f);
%          [psi,psif]=morlwave(n,fmax,f);  
%
%   'morlwave --t' run a test.
%   'morlwave --f' generates a sample figure.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 F. Rekibi and J. M. Lilly 
%                         --- type 'help jlab_license' for details  

if strcmpi(N,'--t')
    morlwave_test;return
end

if strcmpi(N,'--f')
    type makefigs_morlwave
    makefigs_morlwave;
    return
end

if nargin<4
    str='bandpass';
end
%/********************************************************
%Enforce convention of high frequencies first
fs=fs(:);
if length(fs)>1
  if fs(2)-fs(1)>0
    fs=flipud(fs);
  end
end
%\********************************************************

for i=1:length(fs)
    w(:,i)= morlwave1(N,abs(fs(i)),abs(fo),str);
    if fs(i)<0
        w(:,i)=conj(w(:,i));
    end
end

wshifted=0*w;
for i=1:size(w,2)
  wshifted(:,i)=fftshift(w(:,i));
end
W=fft(wshifted);

%---------------------------------------------------------------

function  [w]=morlwave1(N,fs,fo,str)

if fo <=1
    error('MORLWAVE frequency FMAX should be > 1')
end


s=frac(fo,fs);
a=morla(fo,str);

if ~isempty(strfind(str,'ene'))
    a=a.*frac(1,sqrt(s));
elseif ~isempty(strfind(str,'ban'))
    a=a.*frac(1,s);
end

t=(1:N)';   % Choose dt = 1
t=t-mean(t);
t=t./s;

nu=morlfreq(fo);
w=a.*exp(-frac(1,2).*t.^2).*(rot(nu.*t)-exp(-frac(1,2).*nu.^2));


%---------------------------------------------------------------

       
function[a]=morla(fm,str)
%MORLA  Returns the Morlet wavelet amplitude coefficient "a".
%
%   A=MORLA(FMAX), given a peak frequency FMAX, returns the Morlet
%   wavelet amplitude called "A_{NU}" by Lilly and Olhede (2010).
%
%   By default, A is chosen such that the maximum of the frequency-
%   domain wavelet is equal to 2, the ``bandpass normalization.''
%
%   A=MORLA(FMAX,'energy') instead returns the coefficient giving 
%   the wavelet unit energy.  
%
%   Usage: a=morla(fmax);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008 J.M. Lilly --- type 'help jlab_license' for details

fnu=morlfreq(fm);

if strcmpi(str(1:3),'ban')
    a=2*frac(fm,fnu).*exp(frac(1,2).*(fm-fnu).^2)*frac(1,sqrt(2*pi));
elseif strcmpi(str(1:3),'ene')
    a=frac(1,sqrt(sqrt(pi))).*frac(1,sqrt(1+exp(-(fnu).^2)-2.*exp(-frac(3,4).*(fnu).^2)));  
end


%\**************************************************************************

function []= morlwave_test
  
N=1028; % nombre de points
fs=logspace(log10(1/4),log10(10*1/4),10)';   %change to log space
[w,W]=morlwave(N,10/4,fs,'energy');
t=1:length(w);
t=t-mean(t);

tol=5e-3;

b=aresame(mean(w,1),0*w(1,:),tol);
reporttest('MORLWAVE zero mean',b)

b=aresame(sum(abs(w).^2,1),1+0*w(1,:),tol);
reporttest('MORLWAVE unit energy for energy normalization',b)

[w,W]=morlwave(N,10/4,fs(end),'bandpass');
b=aresame(maxmax(abs(W)),2,tol);
reporttest('MORLWAVE maximum of 2 for bandpass normalization',b)

fs=1./logspace(log10(5),log10(40))'; 
N=1023;
[psi,Psi]=morlwave(N,2*pi/4,fs,'bandpass');
reporttest('MORLWAVE Fourier transform, N odd, FS positive',aresame(abs(fft(psi,[],1)),abs(Psi),1e-8))
[psi,Psi]=morlwave(N+1,2*pi/4,fs,'bandpass');
reporttest('MORLWAVE Fourier transform, N even, FS positive',aresame(abs(fft(psi,[],1)),abs(Psi),1e-8))
[psi,Psi]=morlwave(N,2*pi/4,-fs,'bandpass');
reporttest('MORLWAVE Fourier transform, N odd, FS negative',aresame(abs(fft(psi,[],1)),abs(Psi),1e-8))
[psi,Psi]=morlwave(N+1,2*pi/4,-fs,'bandpass');
reporttest('MORLWAVE Fourier transform, N even, FS negative',aresame(abs(fft(psi,[],1)),abs(Psi),1e-8))
