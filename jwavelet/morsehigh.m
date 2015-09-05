function[f]=morsehigh(varargin)
%MORSEHIGH  High-frequency cutoff of the generalized Morse wavelets.
%
%   MORSEHIGH is a low-level function called by MORSESPACE.
%
%   FALPHA=MORSEHIGH(GAMMA,BETA,ALPHA) returns the high frequency cutoff
%   FALPHA of the generalized Morse wavelet specified by GAMMA and BETA, 
%   with cutoff level FALPHA.
%
%   Specifically, if PSI is the wavelet and PSIMAX is its maximum value, 
%   then FALPHA is the highest *radian* frequency at which 
%
%      PSI(FALPHA)/PSIMAX > ALPHA.
%
%   This gives a way to choose the high-frequency cutoff in the wavelet
%   transform.  See Lilly and Olhede (2009d) for details. 
%  
%   The input parameters may either all be scalars, or GAMMA and BETA
%   may be scalars of the same size with scalar ALPHA.
%   ___________________________________________________________________
%
%   Precision vs. speed
%
%   MORSEHIGH(..., N) uses 1/N times the peak frequency MORSEFREQ as the
%   numerical interval.  N=100 is the default; choose a smaller value
%   for faster speed but diminished precision. 
%   ___________________________________________________________________
%  
%   See also MORSEFREQ, MORSEWAVE.
%
%   Usage: falpha=morsehigh(ga,be,alpha);
%          falpha=morsehigh(ga,be,alpha,N);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2015 J.M. Lilly --- type 'help jlab_license' for details
 

gamma=varargin{1};
beta=varargin{2};
alpha=varargin{3};
if nargin>3
    N=varargin{4};
end
    
if nargin~=3&&nargin~=4;
    error('MORSEHIGH takes either three or four input arguments.')
end

ompeak=morsefreq(gamma,beta);
N=100;

dom=vrep(ompeak/N,10*N,3);
dom(:,:,1)=ompeak;

ommat=cumsum(dom,3);

amat=vrep(morseafun(gamma,beta),10*N,3);
betamat=vrep(beta,10*N,3);
gammamat=vrep(gamma,10*N,3);

morse=frac(1,2)*amat.*(ommat.^betamat).*exp(-ommat.^gammamat);

%The "+0" is to convert the logical into a numerical value
kk=vsum(0+(morse>alpha),3);
%temp=zeros(size(morse));
%temp(morse>alpha)=1;
%kk=vsum(temp,3);

ii=vrep((1:size(gamma,1))',size(gamma,2),2);
jj=vrep(1:size(gamma,2),size(gamma,1),1);

index=sub2ind(size(morse),ii,jj,kk);

%f=frac(ommat(index),ompeak);
f=ommat(index);


