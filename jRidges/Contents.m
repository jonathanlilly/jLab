% jRidges:  Wavelet ridge analysis of modulated oscillatory signals
%  
%  Top-level functions
%   ridgewalk  - Extract wavelet transform ridges, including bias estimates. 
%   ridgemap   - Maps ridge quantities back onto the time series.            
%   instmom    - Univariate and multivariate instantaneous moments.          
%
%  Ridge utilities
%   ridgetrim   - Trim edge effect regions from wavelet ridges.
%   ridgelen    - Wavelet ridge length expressed as number of full cycles.    
%   periodindex - Returns time index in increments of instantaneous period.  
%   ridgemult   - Ridge multiplicity, the number of simultaneous ridges present.
%
%  See also jEllipse, jWavelet, jOceans

%   Low-level functions
%   isridgepoint - Finds wavelet ridge points using one of several criterion.
%   lininterp  - Fast linear interpolation for arbitrary-sized arrays.       
%   quadinterp - Fast quadratic interpolation for arbitrary-sized arrays.    
%   ridgechains - Forms ridge curves by connecting transform ridge points.   
%   ridgeinterp - Interpolate quantity values onto ridge locations.          

help jRidges
