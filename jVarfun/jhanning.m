function[window]=jhanning(n)
%JHANNING  Hanning window.
%
%   JHANNING(N) returns a length N Hanning window.
%
%   JHANNING is a low-level function used by VFILT.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2019 J.M. Lilly --- type 'help jlab_license' for details    

window=frac(1,2)*(1-cos(frac(2*pi*(1:ceil(n/2)),n+1)))';

if iseven(n)
    window=[window;flipud(window)];
else
    window=[window(1:end-1);flipud(window)];
end


