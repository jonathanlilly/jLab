function[monthbar,monthstd]=monthstats(varargin)
%MONTHSTATS  Mean month and standard deviation using circular statistics.
%
%   [MFBAR,MFSTD]=MONTHSTATS(NUM,DIM) where NUM is an ND array of dates in 
%   Matlab's DATENUM format, finds the mean and standard deviation of the
%   month number in NUM over dimension DIM using circular statistics.
%
%   MFBAR and MFSTD are the mean and standard deviations of the month dot
%   fraction, defined such that 1.0 is January first and 12.99999 is the 
%   end of December, as computed from YEARFRAC.
%
%   MONTHSTATS works by converting month number, computed using YEARFRAC, 
%   into angles, finding the circular statistics, and converting back. 
%
%   For a complex-valued, unit magnitude quanity z, one can find the 
%   statistics of its angles through R e^{1i* Theta} = <z> where "<>" is
%   an averaging operator.  Then Theta is the average angle, and the
%   standard deviation of angles is sqrt(-2 ln R).  
%
%   See e.g. https://en.wikipedia.org/wiki/Directional_statistics.
%
%   Usage: [mfbar,mfstd]=monthstats(num,dim);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2020 J.M. Lilly --- type 'help jlab_license' for details
 
%if strcmp(varargin{1}, '--t')
%    monthstats_test,return
%end

num=varargin{1}; 
dim=varargin{2}; 

[yf,mf]=yearfrac(num); 
monthangle=frac(2*pi,360)*(90+360*(1-mf)/12); %North = January 1 
monthangle=angle(rot(monthangle));
mz=vmean(rot(monthangle),dim);


thetabar=angle(mz);
thetastd=sqrt(-2*log(abs(mz)));

monthbar=1-frac(12,360)*(deg360(thetabar*frac(360,2*pi))-90);
monthbar(monthbar<1)=monthbar(monthbar<1)+12;
monthstd=frac(12,360)*thetastd*frac(360,2*pi);

%function[]=monthstats_test
