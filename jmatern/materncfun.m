function[c]=materncfun(al,ga)
%MATERNCFUN Returns the normalization or C-function for a Matern process.
%
%   MATERNCFUN is a low-level function used in the JMATERN package.
%
%   MATERNCFUN(ALPHA) returns the Matern C-function defined in Eqn. (44) of
%   Lilly et al. (2017).
%
%   MATERNCFUN(ALPHA,GAMMA) returns C-function for the two-parameter 
%   generalized Matern process.  See MATERNSPEC for details.
%
%   Usage: c=materncfun(alpha);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015--2020 J.M. Lilly --- type 'help jlab_license' for details
 
%c=frac(1,2*pi)*frac(gamma(1/2)*gamma(alpha-1/2),gamma(alpha));
if nargin==1
    if al>1/2
        c=frac(1,2*sqrt(pi))*frac(gamma(al-1/2),gamma(al));
    else
        c=frac(2.^(abs(al-1/2)-abs(al)),sqrt(2*pi)).*...
            frac(gamma(abs(al-1/2)),gamma(abs(al)));
    end
else
    c=frac(1,2*pi*ga).*beta(frac(1,2*ga),frac(1,ga).*(al-1/2));
end
%Verified that these two forms are the same for al>1/2
%al1=[1/2:.01:10]';
%gamma1=[0:.01:10]';
%[al,gamma]=meshgrid(al1,gamma1);
%c=materncfun(al,gamma);
%jpcolor(al1,gamma1,log10(c)),caxis([-1 0])