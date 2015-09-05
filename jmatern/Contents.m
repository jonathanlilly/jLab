% jMatern:  Parametric spectral analysis based on the Matern process
%
% Please note that this toolbox is still in development.  
%
% Top-level functions
%   maternoise - Realizations of the Matern random process and variations.  [with A. Sykulski]      
%   maternspec - Fourier spectrum of the Matern random process and variations.                                
%   materncov  - Autocovariance of the Matern random process and variations.                         
%   maternimp  - Impulse response function for the Matern random process.                      
%
% Other utilities
%   blurspec   - Returns the blurred and aliased spectrum given the autocovariance.
%   fminsearchbnd: - FMINSEARCH, but with bound constraints by transformation. [By J. D'Errico]
%
% Low-level Matern functions
%   materncfun - Returns the normalization function C_ALPHA for a Matern process.
%   maternchol - Cholesky decomposition of the Matern covariance and variations.
%   maternedge - Long-time cutoff edge for the Matern impulse response function.               
%
% See also jSpectral.

help jMatern
