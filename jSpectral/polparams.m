function[trS,p,alpha,beta,preal,gamma]=polparams(varargin)
% POLPARAMS  Spectral matrix polarization parameters.
%
%   [E,P,ALPHA,BETA]=POLPARAMS(S) where S is a 2x2 spectral matrix returns
%   the energy E and the polarization paramters P, ALPHA, and BETA.
%
%   If S is a 2 x 2 x M matrix, then all output variables are Mx1 column
%   vectors.  More generally, if S is 2 x 2 x M x ... N, then all output
%   variables are M x ... N.
%  
%   [...]=POLPARAMS(SXX,SXY,SXY) also works. In this case all input 
%   variables are matrices of the same size. 
%  
%   [E,P,ALPHA,BETA,PREAL,GAMMA]=POLPARAMS(S) also optionally retuns the
%   polarization ratio of the real part of S and the coherence GAMMA.
%
%   See also SPECDIAG.
%
%   Usage:  [E,p,alpha,beta]=polparams(s);
%           [E,p,alpha,beta,preal,gamma]=polparams(s);
%           [E,p,alpha,beta,preal,gamma]=polparams(sxx,syy,sxy);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2012 J.M. Lilly --- type 'help jlab_license' for details        

breshape=0;
if nargin==1
  S=varargin{1};
  S11=S(1,1,:);
  S22=S(2,2,:);
  S12=S(1,2,:);
  S21=S(2,1,:);
elseif nargin==3
  S11=varargin{1};
  S22=varargin{2};
  S12=varargin{3};
  S21=conj(varargin{3});
end


detS=real(S11.*S22-S12.*S21);
trS=S11+S22;

p= sqrt(1 - frac(4.*detS,trS.*trS));


alpha=frac(S11-S22,trS);
beta=frac(2.*S12,trS);

if nargout>3
  detrS=real(S11.*S22-real(S12).^2);
  preal=sqrt(1 - frac(4.*detrS,trS.*trS));
end

if nargout>4
  gamma=frac(S12,sqrt(S11.*S22));
end
