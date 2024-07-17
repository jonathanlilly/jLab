function[fs]=morsespace(varargin)
%MORSESPACE  Logarithmically-spaced frequencies for generalized Morse wavelets.
%
%   F=MORSESPACE(GAMMA,BETA,N) generates a frequency array for the 
%   generalized Morse wavelet transform of an N-point time series. The 
%   wavelets are specified by GAMMA and BETA.  
%
%   LOG(F) is uniformly spaced, following convention for wavelet analysis.
%
%   F has units of *radians* per sample point, and is a column vector with 
%   the frequencies arranged in decending order. 
%   
%   In this usage, the frequencies F are determined using default settings,
%   described below, which should be appropriate for most applications.  
%
%   Additional control over the frequency array F can be obtained using the
%   following alternate usages.
%
%   For details on this calculation, see Lilly (2016).
%   __________________________________________________________________
%
%   High- and low-frequency specification
%
%   F=MORSESPACE(GAMMA,BETA,HIGH,LOW) explicitly sets the high-frequency
%   and low-frequency cutoffs for the frequency array.  
%
%   The first (largest) value of F is then just smaller than HIGH and the 
%   smallest is just larger than LOW.
%
%   HIGH and LOW have units of *radian* per unit sample point.
%   __________________________________________________________________
%
%   High-frequency cutoff
%
%   The highest frequency can be set to be the minimum of a specified value
%   and a cutoff frequency based on a Nyquist overlap condition.
%
%   F=MORSESPACE(GAMMA,BETA,{ETA,HIGH},LOW) sets the highest frequency 
%   to be the minimum of the specified value HIGH, and the largest 
%   frequency for which the wavelet will satisfy the threshold level ETA. 
%
%   Here ETA be a number between zero and one specifying the ratio of a
%   frequency-domain wavelet at the Nyquist frequency to its peak value.
%
%   Note that in this usage, {ETA,HIGH} is a cell array with two entries.
%
%   The simplified usage F=MORSESPACE(GAMMA,BETA,N) corresponds to the 
%   choice ETA=0.1, so that by default, the highest-frequency wavelet 
%   will decay to at least 10% of its peak value at the Nyquist frequency.
%   __________________________________________________________________
%
%   Low-frequency cutoff
%
%   The lowest frequency can be set to a cutoff frequency based on an 
%   endpoint overlap condition.
%
%   F=MORSESPACE(GAMMA,BETA,HIGH,{P,N}) sets the lowest frequency such that
%   the lowest-frequency wavelet will reach some number P times its central 
%   window width at the ends of the time series. 
%  
%   P is called the packing number.  A choice of P=1 corresponds to 
%   roughly 95% of the time-domain wavelet energy being contained within 
%   the time series endpoints for a wavelet at the center of the domain.
%
%   F=MORSESPACE(GAMMA,BETA,HIGH,{P,N,LOW}) alternately chooses the maximum 
%   of the R-level cutoff frequency, and a specified low frequency LOW.
%   
%   The simplified usage F=MORSESPACE(GAMMA,BETA,N) corresponds to the
%   default value P=5.  At the lowest frequency, five wavelets will then
%   fit into the time series, with 5% energy overlap between them.
%   __________________________________________________________________
%
%   Wavelet density
%
%   F=MORSESPACE(GAMMA,BETA,HIGH,LOW,D) controls the number of points in 
%   the frequency array through the 'density' D. 
%
%   Higher values of D mean more overlap in the frequency domain. The
%   default value of the density is D=4.
%
%   When D=1, the peak of one wavelet is located at the half-power points 
%   of the adjacent wavelet. D=4 means that four other wavelets will occur 
%   between the peak of one wavelet and its half-power point. 
%   __________________________________________________________________
%
%   See also MORSEWAVE, WAVETRANS.
%  
%   Usage: f=morsespace(ga,be,N);
%          f=morsespace(ga,be,high,low,D);
%          f=morsespace(ga,be,{eta,high},low,D);
%          f=morsespace(ga,be,{eta,high},{p,N},D);
%          f=morsespace(ga,be,{eta,high},{p,N,low},D);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2023 J.M. Lilly --- type 'help jlab_license' for details


%   Note then when HIGH and LOW are explicitly input, the time series
%   length N is no longer input as the first argument. 

if strcmpi(varargin{1},'--t')
  morsespace_test;return
end

if strcmpi(varargin{1},'--f')
  type makefigs_morsespace;
  makefigs_morsespace
  return
end


%/***************************************************************
%Sort out input arguments 
% if ischar(varargin{end})
%     str=varargin{end};
%     varargin=varargin(1:end-1);
% else
%     str='peak';
% end

ga=varargin{1};
be=varargin{2};
D=4;    

if length(varargin)==3
    N=varargin{3};
    
    %Default choices
    low={5,N};
    high={0.1,pi};
else
    high=varargin{3};
    low=varargin{4};
    if length(varargin)==5
        D=varargin{5};
    end
end

%\***************************************************************


if iscell(high)
    %Recall pi is Nyquist
    high=min(high{2},morsehigh(ga,be,high{1}));
    %high
end

if iscell(low)
    if length(low)==2
        low{3}=0;
    end
    low=max(low{3},morsespace_low(ga,be,low{1},low{2}));
end

r=1+frac(1,D*morseprops(ga,be));
N=floor(frac(log(frac(high,low)),log(r)));
fs=high*ones(N+1,1)./r.^[0:N]';

% if length(fs)>1000
%     warning('F is longer than 1000 points... ')
% end

% if strfind(str,'ene')
%      [fm,fe,fi] = morsefreq(ga,be);
%      fs=frac(fe,fm).*fs;
% end
 
function[fmin]=morsespace_low(ga,be,r,N)

p=morseprops(ga,be);
fmin=frac(2*sqrt(2)*p*r,N);

function[]=morsespace_test
morsespace_hightest;
morsespace_lowtest;

function[]=morsespace_hightest
%Test for high-frequncy cutoff

ga=3;
be=4;
eta=0.1;
ompeak=morsefreq(ga,be);
N=100;

dom=ompeak/N+zeros(4*N,1);
om=cumsum(dom,1);

%morse=frac(1,2)*morseafun(ga,be).*(om.^be).*exp(-om.^ga);

fhigh=morsehigh(ga,be,eta);
[psi,morse] = morsewave(10000,ga,be,fhigh);
reporttest('MORSESPACE high-frequency cutoff',aresame(morse(end/2),eta,1e-3))

function[]=morsespace_lowtest
%Test for low-frequency cutoff

N=1001;
ga1=(1/3:1:11);
be1=(1:1:10);
%ga1=(1/3:.5:11);
%be1=(1:.5:10);

[ga,be]=meshgrid(ga1,be1);
vcolon(ga,be);
%[a,sigt,sigo]=morsebox(ga,be);

psi=zeros(N,length(ga));
for i=1:length(ga)
    fs{i}=morsespace(ga(i),be(i),{0.95,pi},{1,N/10},4);
    psi(:,i)=morsewave(N,ga(i),be(i),fs{i}(end),'energy');
end
%figure,plot(sum(squared(psi((end+1)/2-50:(end+1)/2+50,:)),1))
bool=aresame(median(sum(squared(psi((end+1)/2-50:(end+1)/2+50,:)))),0.95,0.01);

reporttest('MORSESPACE low-frequency cutoff, P approximates 95% energy',bool)

%fs=morsespace(2,2,10000);
%psi=morsewave(10000,2,2,fs(end),'energy');
%Can make a figure with the wavelet shifted by 2000 points each time


