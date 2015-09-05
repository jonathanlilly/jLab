function[r]=ellrad(varargin)
%ELLRAD  Average and instantaneous ellipse radius. 
%
%   R=ELLRAD(KAPPA,LAMBDA,PHI,STR) where KAPPA and LAMBDA are the amplitude 
%   and linearity of a time-varying ellise, and PHI is its time-varying 
%   phase, returns various measures of the ellipse 'radius'.
%   
%   STR determines the radius quantity to be output.
%
%       STR       Symbol     Description
%       ----------------------------------------------------------------
%      'geo'      RM         Geometric mean radius = SQRT(A*B)    
%      'ave'      RBAR       Period-averaged distance from the origin          
%      'ins'      RI         Instantaneous distance from the origin
%  
%   See Lilly and Gascard (2006) for details.
%
%   If STR is omitted, the geometric mean radius is returned by default.
%
%   RM=ELLRAD(KAPPA,LAMBDA) with PHI omitted also works for computing RM.
%   ____________________________________________________________________
%   
%   Cell array input/output
%
%   If ELLRAD is given cell array input, it returns cell array output.
%
%   Thus KAPPA, LAMBDA, and PHI may each be cell arrays of the same size, 
%   where each element in the cell array is a numerical array.  
%
%   The output variables will be also cell arrays of this size.
%   ____________________________________________________________________
%
%   See also ELLVEL, ELLPARAMS, ELLDIFF.
%
%   Usage:  rm=ellrad(kappa,lambda);
%           rbar=ellrad(kappa,lambda,phi,'average');
%           ri=ellrad(kappa,lambda,phi,'instantaneous');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2005--2014 J.M. Lilly --- type 'help jlab_license' for details    
  

if ischar(varargin{end})
    str=varargin{end}(1:3);
    varargin=varargin(1:end-1);
else
    str='geo';
end

kappa=varargin{1};
lambda=varargin{2};

if length(varargin)==2
    if ~aresame(str,'geo')
           error('Only the geometric mean radius can be computed without PHI input.')
    else
        phi=kappa;  %It doesn't matter what we set this to as it will not be used.
    end
else
    phi=varargin{3};
end

if ~isempty(kappa)
    if ~iscell(kappa)
        r=ellrad_one(kappa,lambda,phi,str);
    else
        r=kappa;
        for i=1:length(kappa)
            r{i}=ellrad_one(kappa{i},lambda{i},phi{i},str);
        end
    end
else
    r=kappa;
end

function[r]=ellrad_one(kappa,lambda,phi,str)
if aresame(str,'geo')
     [a,b]=kl2ab(kappa,lambda);
     r=sqrt(a.*abs(b));   %r=kappa^2*sqrt(1-lambda^2);
elseif aresame(str,'ave')
     ecc=ecconv(lambda,'lin2ecc');
     [K,E]=ellipke(ecc.^2);
     r=frac(2*kappa,pi).*sqrt(1+abs(lambda)).*E;
elseif aresame(str,'ins')
     r=kappa.*sqrt(1+abs(lambda).*cos(2*phi));
end

