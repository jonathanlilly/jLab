function[varargout]=maxprops(varargin)
%MAXPROPS  Returns properties of wavelet transform maxima.
%
%   This function is part of 'element analysis' described in Lilly (2017), 
%   "Element analysis: a wavelet-based method for analyzing time-localized
%   events in noisy time series", available at www.jmlilly.net.
%
%   [C,RHO,FRHO]=MAXPROPS(WW,FF,GAMMA,BETA,MU) returns estimated signal
%   properties as inferred from values of the wavelet the transfrom maxima.
%
%   A signal believed to consist of rescaled versions of a (GAMMA,MU)
%   generalized Morse wavelet is transformed with a (GAMMA,BETA) wavelet 
%   at frequency levels FS, as output by WAVETRANS.  The transform maxima
%   have values WW and scale frequencies FF, as output by TRANSMAX.
%   
%   C is the estimated complex-valued coefficients of the (GAMMA,MU)
%   wavelets that are taken to correspond to these maxima, RHO are their
%   scales, and FRHO are their scale frequencies.
%
%   The output fields are all arrays having the same length as WW and FF.
%
%   For details, see Lilly (2017).
%
%   See also TRANSMAX, ISOMAX, MAX2EDDY.
%
%   Usage: [C,rho,frho]=maxprops(ww,ff,ga,be,mu);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2017--2023 J.M. Lilly --- type 'help jlab_license' for details
  
str='band';
if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
end


ww=varargin{1};
ff=varargin{2};
ga=varargin{3};
be=varargin{4};
mu=varargin{5};

%fhat=ff*frac(morsefreq(ga,mu),morsefreq(ga,be)).*frac(be,mu+1).^(1./ga);
%fhat=ff*frac(mu,mu+1).^(1./ga); %This form does not work if mu=0
%chat=2*ww./morsezeta(ga,be,mu);

if strcmpi(str(1:3),'ene')
    fhat=ff*frac(morsefreq(ga,mu),morsefreq(ga,be)).*frac(be+1/2,mu+1/2).^(1./ga);
elseif strcmpi(str(1:3),'ban')
    fhat=ff*frac(morsefreq(ga,mu),morsefreq(ga,be)).*frac(be,mu+1).^(1./ga);
end

%fhat=ff*frac(mu,mu+1).^(1./ga); %This form does not work if mu=0
chat=ww./morsezeta(ga,be,mu,str);

varargout{1}=chat;
varargout{2}=morsefreq(ga,mu)./fhat;
varargout{3}=fhat;


function[zeta]=morsezeta(varargin)
%MORSEZETA  The maximum of the Morse wavelet transform of another wavelet.
%
%   ZETAMAX=MORSEZETA(GAMMA,BETA,MU) returns the maximum value ZETAMAX of 
%   the generalized Morse wavelet transform of an order MU wavelet with an
%   order BETA wavelet, both of which are in the GAMMA family.
%
%   This is used in converting wavelet transform values into the weighting
%   coefficients of the analyzed wavelets, see Lilly (2017).
%
%   Usage: zetamax=morsezeta(ga,be,mu);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2017 J.M. Lilly --- type 'help jlab_license' for details

str='band';
if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
end

ga=varargin{1};
be=varargin{2};
mu=varargin{3};

if strcmpi(str(1:3),'ene')
    numer=frac(be+1/2,mu+1/2).^((be+1/2)./ga);
    denom=(frac(be+1/2,mu+1/2)+1).^((be+mu+1)./ga);
elseif strcmpi(str(1:3),'ban')
    numer=frac(be,mu+1).^(be./ga);
    denom=(frac(be,mu+1)+1).^((be+mu+1)./ga);
end

B=frac(numer,denom);
zeta=frac(1,2*pi*ga).*morseafun(ga,be,str).*morseafun(ga,mu,str)...
    .*gamma(frac(be+mu+1,ga)).*B;

