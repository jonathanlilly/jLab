function[varargout]=jlab_allhelp(varargin)
%   jCell:  Tools for operating on cell arrays of column vectors 
%   
%   Basic mathematical operations 
%     cellabs    - Absolute value of each element in a cell array.                     
%     cellmax    - Maximum of each element in a cell array.                            
%     cellmin    - Minimum of each element in a cell array.  
%     cellsum    - Sum of each element a cell array, possibly weighted. 
%     cellmean   - Mean value of each element a cell array, possibly weighted.   
%     cellstd    - Standard deviation of each element a cell array, possibly weighted. 
%     cellmed    - Median value of each element a cell array. 
%     cellreal   - Real part of each element in a cell array.                          
%     cellimag   - Imaginary part of each element in a cell array. 
%     cellpair   - Complex pairing for elements in a cell array. 
%     celllog10  - Base ten logarithm of each element in a cell array. 
%     celladd    - Addition acting on each element in a cell array.                    
%     cellmult   - Multiplication acting on each element in a cell array.    
%     celldiv    - Division acting on each element in a cell array. 
%     celldot    - Dot product for arrays of column vectors. 
%   
%   Reshaping, indexing, and sizes 
%     cell2col   - Converts cell arrays of numeric arrays into 'column-appended' form. 
%     col2cell   - Converts 'column-appended' data into cell arrays of numeric arrays. 
%     cellindex  - Applies a cell array of indices to a cell array of column vectors.  
%     cellchunk  - Converts cell array data into uniform length 'chunks'.              
%     cellength  - Length of each element in a cell array.                             
%     cellsize   - Size of each element in a cell array along specified dimension. 
%     cellget    - Indexes a cell array of numerical arrays by ID number. 
%     cellimit   - Limits the ranges of times in a cell array of numerical arrays. 
%   
%   Data processing 
%     cellstrip  - Strips INF values from the beginnings or ends of cell arrays. 
%     cellsplit  - Splits cell arrays of numeric arrays at data gaps. 
%     cellprune  - Removes all empty cells, or cells less than a specified length. 
%     cellpack   - Removes all INF values from cell arrays of numeric arrays. 
%     cellfill   - Fills missing data marked by INFs in a cell array. 
%     cellgrid   - Interpolate a cell array of numeric arrays onto a regular grid. 
%     cellfirst  - Returns the first element of each entry in a cell array. 
%   
%   Plotting 
%     cellplot   - Rapidly plot all element of a cell array of numeric arrays.  
%   
%   See also jVarfun. 
%   __________________________________________________________________ 
%  
%   jCommon:  Useful general-purpose functions, common to other JLAB modules 
%   
%   Array functions 
%     aresame    - Test whether two N-D arrays are the same.                         
%     arrayify   - Converts a set of scalars or arrays into column arrays.             
%     blocklen   - Counts the lengths of 'blocks' in an array.                         
%     blocknum   - Numbers the contiguous blocks of an array.    
%     lnsd       - Last non-singleton dimension of an array.                           
%     nonnan     - Return all non-NAN elements of an array.                            
%     vectmult   - Matrix multiplication for arrays of vectors.     
%   
%   File and directory tools 
%     ncload       - Load all variables from a NetCDF file and unpack any columns. 
%     commentlines - Returns the comment lines from m-files.                           
%     findfiles    - Returns all files in a directory with a specified extension.        
%     findpath     - Returns the full pathname of a directory on the Matlab search path. 
%     jhelp        - Opens linked JLAB help files in Matlab's internal web browser.      
%     whichdir     - Returns directory name containing file in search path.              
%     standalone   - Create stand-alone version of an m-file, including dependencies. 
%   
%   Mathematical aliases 
%     choose     - Binomial coefficient: CHOOSE(N,K) =  N!K!/(N-K)! 
%     frac       - Fraction: FRAC(A,B)=A./B                                                     
%     imlog      - Imaginary part of log: IMLOG(X)=UNWRAP(IMAG(LOG(X)))                                       
%     iseven     - True for even integer values; false otherwise.                      
%     isodd      - True for odd integer values; false otherwise.  
%     oprod      - Outer product:  OPROD(X,Y)=X*CONJ(Y')                               
%     res        - Residual after flooring:  RES(X)=X-FLOOR(X)                         
%     rot        - Complex-valued rotation:  ROT(X)=EXP(SQRT(-1)*X)   
%     squared    - Squares the modulus of its argument:  SQUARED(X)=ABS(X).^2 
%    
%   Matrices, polynomials, and special functions 
%     jmat2      - 2x2 rotation matrix through specified angle.                        
%     jmat3      - 3x3 rotation matrix through specified angle.   
%     tmat       - 2x2 complex grouping matrix.  TMAT = [1  i; 1 -i] / SQRT(2) 
%     bellpoly   - Complete Bell polynomials.   
%     jhermpoly  - Hermite polynomials. [with F. Rekibi]    
%     jhermfun   - Orthonormal Hermite functions. [with F. Rekibi]    
%     chisquared - The chi-squared distribution. 
%   
%   Dataset organization as structures 
%     make        - Create a structure containing named variables as fields.   
%     matsave     - Create and save structure of variables as a mat-file.  
%     use         - Copies structure fields into named variables in workspace.          
%     catstruct   - Concatenates the array elements of a cell array of structures. 
%     structindex - Applies an index to all array-valued fields in a structure. 
%   
%   Statistics 
%     cum2mom   - Convert cumulants to moments.     
%     mom2cum   - Convert moments to cumulants.                                       
%     pdfprops  - Mean and variance associated with a probability distribution.   
%     simplepdf - Gaussian, uniform, Cauchy, and exponential pdfs.                    
%   
%   Filling bad data points 
%     fillbad    - Linearly interpolate over bad data points.   
%   
%   Date, time, and units 
%     yearfrac    - Converts a DATENUM into 'year.fraction' and 'month.fraction'. 
%     cms2kmd     - Converts centimeters per second to kilometers per day. 
%   __________________________________________________________________ 
%  
%   jData: Routines for generating oceanographic datasets 
%   
%     make_goldaltimetry - Create a synthetic altimeter dataset and gridded statistics. 
%     make_goldocean    - Create quarter-degree gridded files from GOLD simulations. 
%     make_betaeddyone  - Create the BetaEddyOne nonlinear QG eddy simulation. 
%     make_gomed        - Create a drifter-derived eddy database for the Gulf of Mexico. 
%     make_gulfdrifters - Create a drifter dataset for the Gulf of Mexico. 
%     make_gulfflow     - Create a gridded velocity dataset for the Gulf of Mexico. 
%     make_jason        - Create a reformatted version of the Beckley altimeter dataset. 
%   __________________________________________________________________ 
%  
%   jEllipse:  Analysis of modulated elliptical, or bivariate, signals 
%   
%    Elliptical signal properties 
%     ellband    - Bandwidth of modulated elliptical signals in two or three dimensions.  
%     ellparams  - Ellipse parameters of a modulated bivariate or trivariate oscillation. 
%     ellvel     - Average and instantaneous ellipse velocities.     
%     ellrad     - Average and instantaneous ellipse radius.   
%     ellpol     - Polarization parameters of an elliptical signal. 
%   
%    Elliptical signal generation and conversions 
%     ellsig     - Creates a modulated elliptical signal in two or three dimensions.      
%     sig2latlon - Converts an oscillatory signal to lat/lon displacements.   
%     elldiff    - Differentiation of modulated elliptical signals.     
%             
%    Plotting tools 
%     ellcurves   - Returns curves corresponding to specified ellipse properties.          
%     ellipseplot - Plot ellipses. 
%   
%    Basic ellipse properties  
%     ab2kl      - Converts A and B to ellipse parameters Kappa and Lambda.               
%     kl2ab      - Converts ellipse parameters Kappa and Lambda to A and B.               
%     ecconv     - Converts between eccentricity measures.           
%     ellrossby  - Ellipse Rossby number, for oceanographic applications.                 
%   
%    See also jRidges, jSpectral, jWavelet. 
%   __________________________________________________________________ 
%  
%   jFigures:  A low-level directory containing routines for example figures. 
%   __________________________________________________________________ 
%  
%   jGraph:  Plotting tools and refinements for making high-quality figures 
%   
%   High-level graphical post-processing 
%     fixlabels    - Specify precision of axes labels.  
%     fontsize     - Rapidly set title, axes, label, and text fontsizes.               
%     letterlabels - For automatically putting letter labels on subplots.            
%     linecolor    - Set line colors based on a property value within a colormap.                     
%     linestyle    - Rapidly set color, style, and width properties of lines. 
%     linering     - Moves lines through the current line style order.  
%     packfig      - Squeeze together rows and/or columns of the current figure.  
%   
%    Oceanography-specific plotting tools 
%     hodograph  - Generate hodograph plots (simple and fancy).                                    
%     provec     - Generate progressive vector diagrams (simple and fancy).              
%     stickvect  - Plots "stick vectors" for multicomponent velocity time series.   
%     [See also CELLPLOT, REGIONPLOT, TOPOPLOT] 
%   
%   Modified plotting functions 
%     discretecolorbar - Plots a colorbar with discrete variation.   
%     jpcolor          - Modified version of PCOLOR appropriate for cell-centered grids.   
%     patchcontourf    - Generate filled contours using patches, with specified colors. 
%     uvplot           - Plots the real and imaginary parts of a signal on the same axis. 
%   
%   Figure tweaking 
%     jfig       - Shorthand for tweaking figures.    
%     dlines     - Add diagonal lines to a plot.   
%     hlines     - Add horizontal lines to a plot.                                   
%     vlines     - Add vertical lines to a plot.    
%     nocontours - Removes contours from a CONTOURF plot.                            
%     noxlabels  - Remove some or all x-axis tick mark labels.                       
%     noylabels  - Remove some or all y-axis tick mark labels.  
%     xoffset    - Offsets lines in the x-direction after plotting.        
%     yoffset    - Offsets lines in the y-direction after plotting.                  
%     xtick      - Sets locations of x-axis tick marks.                              
%     ytick      - Sets locations of y-axis tick marks.                              
%     ztick      - Sets locations of z-axis tick marks. 
%   
%   Color and colormaps 
%     colorquant  - Sets the color axis according to quantiles of the data. 
%     lansey      - The Lansey modification of Cynthia Brewer's "Spectral" colormap. 
%     haxby       - The Haxby colormap. 
%     seminfhaxby - The seminf-Haxby colormap. 
%   
%   Printing 
%     jprint     - Print to a specified directory and crop the resulting file. 
%     printall   - Print and close all open figures. 
%   
%   Graphics aliases 
%     boxoff     - Sets 'box' property to 'off'.                                     
%     boxon      - Sets 'box' property to 'on'.                                                                   
%     flipmap    - Flips the current colormap upside-down.                           
%     flipx      - Flips the direction of the x-axis                                 
%     flipy      - Flips the direction of the y-axis     
%     inticks    - Sets the 'tickdir' property of the current axis to 'in'. 
%     outticks   - Sets the 'tickdir' property of the current axis to 'out'.          
%     xlin       - Sets x-axis scale to linear.                                       
%     xlog       - Sets x-axis scale to logarithm.                                  
%     ylin       - Sets y-axis scale to linear.                                       
%     ylog       - Sets y-axis scale to log.         
%   
%   Low-level or specialized functions 
%     axeshandles - Returns handles to all axes children. 
%     crop        - Gets rid of whitespace around an image. [by A. Bliss]              
%     linehandles - Finds all line and patch handles from a given set of axes. 
%     linestyleparse - Parses the input string to LINESTYLE. 
%     gulf4plot     - A four-panel circulation plot for the Gulf of Mexico. 
%   __________________________________________________________________ 
%  
%   jMatern:  Parametric spectral analysis based on the Matern process 
%   
%   Top-level functions 
%     maternspec - Fourier spectrum of the Matern random process and variations.                                 
%     materncov  - Autocovariance of the Matern random process and variations.                          
%     maternimp  - Impulse response function for the Matern random process.     
%     maternoise - Realizations of the Matern process and variations, including fBm. [with A. Sykulski] 
%     maternfit  - Parametric spectral fit to the Matern form. [with A. Sykulski] 
%   
%   Other utilities 
%     blurspec   - Returns the blurred and aliased spectrum given the autocovariance. 
%     fminsearchbnd  - FMINSEARCH, but with bound constraints by transformation. [By J. D'Errico] 
%   
%   Low-level Matern functions 
%     materncfun - Returns the normalization or C-function for a Matern process. 
%     maternchol - Cholesky decomposition of Matern and fBm covariances. [with A. Sykulski] 
%     maternedge - Long-time cutoff edge for the Matern impulse response function.                
%   
%   See also jSpectral, makefigs_matern. 
%   __________________________________________________________________ 
%  
%   jMatfiles:  Descriptions of included mat-files, mostly for testing and examples. 
%   
%     bravo94      - Labrador Sea mooring data from Lilly et. al (1999). 
%     drifterbetty - Two-hour resolution drifter trajectory, AOML GDP id #44000. 
%     ebasnfloats  - "Meddy" float data analyzed in Lilly and Olhede (2009,10,11,12). 
%     goldsnapshot - Snapshot from GOLD model, courtesy of Harper Simmons. 
%     labseatpjaos - Altimeter data for the Labrador Sea, analyzed in Lilly (2017). 
%     m1244        - The M1244 mooring from the Labrador Sea boundary current. 
%     npg2006      - Float data analyzed in Lilly and Gascard (2006). 
%     qgsnapshot   - Snapshot of f-plane QG turbulence from a simulation by J. Early. 
%     slasnapshot  - Snapshot of AVISO 1/4 degree sea level anomaly on Jan. 1, 2007. 
%     solomon      - Solomon Islands Earthquake data analyzed in Lilly (2011). 
%     vortex       - QG model of a barotropic jet on a beta plane, by R. K. Scott. 
%   __________________________________________________________________ 
%  
%   jOceans:  Oceanography-specific data and model analysis tools 
%   
%   Conversions for Lagrangian trajectories 
%     latlon2uv   - Converts latitude and longitude to horizontal velocity.   
%     uv2latlon   - Integrates horizontal velocity to give latitude and longitude.   
%   
%   Manipulating Lagrangian trajectories [see also jCell] 
%     trajextract  - Extracts Lagrangian trajectory segments within given region.         
%     trajunwrap   - Unwraps Lagrangian trajectories from a periodic domain.               
%     trajwrap     - Wraps Lagrangian trajectories to fit within a periodic domain.        
%     trajchunk    - Chunks Lagrangian trajectories based on the Coriolis period. 
%     griddrifters - Average drifter velocities onto a space/time 3D grid. 
%   
%   Idealized numerical model tools 
%     psi2fields   - Velocity and other fields from the streamfunction. [with P.E. Isachsen]     
%     periodize    - Returns a doubly periodic version of an input array.    
%   
%   Lagrangian eddy identification and analysis 
%     eddyridges    - Coherent eddy ridges from Lagrangian trajectories.    
%     noisedrifters - Create a noise Lagrangian dataset matching mean and variance. 
%     eddylevels    - Eddy ridge significance levels using the survival function. 
%   
%   Eulerian eddy identification and analysis 
%     inellipse    - Locates points on the interior of ellipses.                           
%     closedcurves - Locate and interpolate closed curves in a possibly periodic domain. 
%     curvemoments - Centroid, area, and many other moments of a closed curve.           
%     divgeom      - Geometric decomposition of eddy vorticity flux divergence.      
%     eddyfit2d    - Least squares fit of 2D velocity data to an eddy profile. 
%     simpleddy     - Streamfunction, velocity, and vorticity for various eddy profiles. 
%   
%   Plotting tools for mooring data 
%     hodograph  - Generate hodograph plots (simple and fancy).                                    
%     provec     - Generate progressive vector diagrams (simple and fancy).              
%     stickvect  - Plots "stick vectors" for multicomponent velocity time series.  
%   
%   Alongtrack altimetry tools 
%     trackextract  - Extracts alongtrack altimetry segments within given region. 
%   
%   Wind-driven surface currents 
%      windtrans  - Ekman-like transfer-functions for the wind-driven response. 
%   
%   NetCDF tools 
%     ncinterp    - One-line interpolation from 3D lat/lon/time field in NetCDF file. 
%   
%   Date and time 
%     monthstats  - Mean month and standard deviation using circular statistics. 
%   
%   Topography tools and data 
%     jtopo.mat   - One-sixth degree global topography, from Smith and Sandwell + IBCAO.               
%     topoplot    - Plot regional or global topography at one-sixth degree resolution. 
%    
%   Low-level functions 
%     besselktilde  - K-type Bessel function after factoring off exponential decay. 
%     besselitilde  - I-type Bessel function after factoring off exponential growth. 
%     curveinterp  - Interpolate a field or its gradient onto a set of curves.            
%     orbitbreaks  - Separate orbit into passes based on turning points.           
%     topo_copyright - Copyright statement for the Smith and Sandwell topography. 
%   
%    See also jData. 
%   __________________________________________________________________ 
%  
%   jPapers: Routines for generating figures from published papers. 
%          
%     makefigs_stokes - Make figures for Lilly et al. (2024).                   
%     makefigs_gulfcensus - Make figures for Lilly and Perez-Brunius (2021b).   
%     makefigs_gulfdrifters - Make figures for Lilly and Perez-Brunius (2021a). 
%     makefigs_transfer - Make figures for Lilly and Elipot (2021).   
%     makefigs_kinematics - Make figures for Lilly (2018).                      
%     makefigs_matern - Make figures for Lilly et al. (2017). 
%     makefigs_element - Make figures for Lilly (2017). 
%     makefigs_superfamily - Make figures for Lilly and Olhede (2012b).         
%     makefigs_multivariate - Make figures for Lilly and Olhede (2012a).        
%     makefigs_ridges - Make figures for Lilly, Scott, and Olhede (2011).       
%     makefigs_trivariate - Make figures for Lilly (2011).  
%     makefigs_analytic - Make figures for Lilly and Olhede (2010b).            
%     makefigs_bandwidth - Make figures for Lilly and Olhede (2010a).    
%     makefigs_asilomar - Make figures for Lilly and Olhede (2009b).            
%     makefigs_morsies - Make figures for Lilly and Olhede (2009a).             
%     makefigs_vortex - Make figures for Lilly and Gascard (2006). 
%   __________________________________________________________________ 
%  
%   jRidges:  Wavelet ridge analysis of modulated oscillatory signals 
%     
%    Top-level functions 
%     ridgewalk  - Extract wavelet transform ridges, including bias estimates.  
%     ridgemap   - Maps ridge quantities back onto the time series.             
%     instmom    - Univariate and multivariate instantaneous moments.           
%   
%    Ridge utilities 
%     ridgetrim   - Trim edge effect regions from wavelet ridges. 
%     ridgelen    - Wavelet ridge length expressed as number of full cycles.     
%     periodindex - Returns time index in increments of instantaneous period.   
%     ridgemult   - Ridge multiplicity, the number of simultaneous ridges present. 
%   
%    See also jEllipse, jWavelet, jOceans 
%   __________________________________________________________________ 
%  
%    jSpectral:  Multitaper spectral analysis, and other time series tools 
%    
%    Multitaper spectral analysis 
%     sleptap    - Calculate Slepian tapers.                                         
%     mspec      - Multitaper power and cross spectra. 
%     mconf      - Confidence intervals for the multitaper spectral estimate. 
%   
%    Multitaper polarization analysis 
%     msvd       - Singular value decomposition for polarization analysis.    
%     polparams  - Spectral matrix polarization parameters.                          
%     specdiag   - Diagonalize a 2 x 2 spectral matrix.         
%   
%   Assorted other transforms  
%     anatrans    - Analytic part of signal.                                          
%     wigdist     - Wigner distribution (alias-free algorithm).    
%    
%   Time series analysis utilities 
%     doublen     - Interpolates a time series to double its length.                  
%     fourier     - The one-sided Fourier frequencies for a given length time series. 
%     sampletimes - Computes mean sampling intervals and their statistics.           
%   
%    Plotting tools  
%     twospecplot - Plots a pair of rotary or Cartesian spectra.   
%   
%    See also jWavelet, jEllipse, jMatern. 
%   __________________________________________________________________ 
%  
%   jSphere:  Tools for latitude, longitude, and spherical geometry 
%   
%   Simple lat / lon conversions 
%     deg180     - Converts degrees to the range [-180,180].                         
%     deg360     - Converts degrees to the range [0, 360]. 
%     degunwrap  - Unwraps arrays given in degrees. 
%     jdeg2rad   - Converts degrees to radians.                                      
%     jrad2deg   - Converts radians to degrees.    
%   
%   Distances and regions  
%     inregion   - Tests whether lat/lon points lie within a specified box.          
%     spheredist - Computes great circle distances on a sphere.             
%   
%   Earth constants and frequencies 
%     radearth   - The radius of the earth in kilometers.    
%     corfreq    - Coriolis frequency in radians per hour.    
%     tidefreq   - Frequencies of the eight major tidal components. 
%   
%   Coordinate transformations  
%     latlon2xy  - Converts latitude and longitude into tangent plane coordinates. 
%     latlon2xyz - Converts latitude and longitude into 3D Cartesian coordinates.    
%     xy2latlon  - Converts tangent plane coordinates into latitude and longitude. 
%     xyz2latlon - Converts 3D Cartesian coordinates into latitude and longitude.    
%     sphere2uvw - Converts a 3D spherical vector to a 3D Cartesian vector.  
%   
%   Differential operators 
%     spherecurl - Curl of a vector field on the surface of a sphere.                
%     spherediv  - Divergence of a vector field on the surface of a sphere.          
%     spheregrad - Gradient of a field on the surface of a sphere.                   
%     spherelap  - Laplacian of a field on the surface of a sphere.    
%   
%   Interpolation 
%     interplatlon  - Interpolation for working with latitude and longitude. 
%     sphereinterp  - Fast linear interpolation on the sphere from non-plaid grids. 
%   
%   Plotting tools 
%     latratio   - Set plot aspect ratio for latitude / longitude plot.     
%     lonshift   - Shifts longitude origin for plotting purposes.                    
%     regionplot - Plots a simple box indicating a latitude / longitude region. 
%   __________________________________________________________________ 
%  
%   jStats:  Fast statistical analysis of large datasets in two dimensions 
%    
%     twodhist   - Two-dimensional histogram.   
%     twodstats  - Mean, variance, and covariance of functions of two variables. 
%     twodmed    - Median value of a function of two variables. 
%   __________________________________________________________________ 
%  
%   jVarfun:  Perform common operations on multiple variables simultaneously. 
%   
%    Sizes and statistics 
%     vsize      - Returns the sizes of multiple arguments.                            
%     vmean      - Mean over non-NaN elements along a specified dimension. 
%     vrms       - Root-mean-square of non-NaN elements along a specified dimension. 
%     vsum       - Sum over non-NaN elements along a specified dimension.               
%     vstd       - Standard deviation over non-NaN elements along a specfied dimension. 
%     vmoment    - Central moment over non-NaN elements along a specfied dimension.     
%     vmedian    - Median over finite elements along a specified dimension.  
%   
%    Reshaping, shifting, swapping 
%     vcolon     - Condenses its arguments, like X(:).                                 
%     vshift     - Cycles the elements of an array along a specified dimension.        
%     vsqueeze   - Squeezes multiple input arguments simultaneously.                   
%     vrep       - Replicates an array along a specified dimension.                    
%     vswap      - VSWAP(X,A,B) replaces A with B in numeric array X.                          
%     vtranspose - Transpose multiple input arguments simultaneously.  
%   
%    Initializing 
%     vzeros     - Initializes multiple variables to arrays of zeros or nans.  
%     vempty     - Initializes multiple variables to empty sets or empty cell arrays.  
%   
%    Operations and indexing 
%     vfilt      - Filtering along rows without change in length.                      
%     vdiff      - Length-preserving first central difference.  
%     vindex     - Indexes an N-D array along a specified dimension.                   
%     vindexinto - Indexes into N-D array along a specified dimension. 
%   
%    Reshaping dataset composed of irregular length segments  
%     col2mat    - Expands 'column-appended' data into a matrix.                       
%     colbreaks  - Insert NANs into discontinuties in a vector.                        
%     mat2col    - Compress NAN-padded matrix data into long columns.     
%   
%    Low-level functions 
%     jhanning   - Hanning window. 
%      
%    See also jCell. 
%   __________________________________________________________________ 
%  
%   jWavelet: Continuous wavelet analysis using generalized Morse wavelets 
%   
%   Top-level wavelet transform functions 
%     wavetrans    - Continuous wavelet transform. 
%     morsewave    - Generalized Morse wavelets of Olhede and Walden (2002).               
%     morsespace   - Logarithmically-spaced frequencies for generalized Morse wavelets.    
%     wavespecplot - Plot of wavelet spectra together with time series. 
%     spheretrans  - Wavelet transform for oscillations on the surface of a sphere. 
%   
%   Details of generalized Morse wavelets  
%     morsearea   - Time-frequency concentration area of Morse wavelets. [with F. Rekibi] 
%     morsebox    - Heisenberg time-frequency box for generalized Morse wavelets.         
%     morsefreq   - Frequency measures for generalized Morse wavelets. [with F. Rekibi]   
%     morseprops  - Properties of the demodulated generalized Morse wavelets.             
%     morsexpand  - Generalized Morse wavelets via time-domain Taylor series. 
%     morsederiv  - Frequency-domain derivatives of generalized Morse wavelets. 
%   
%   Element analysis using generalized Morse wavelets 
%     transmax     - Locates modulus maximum points of wavelet transform. 
%     transmaxdist - Distributions of wavelet transform maxima in noise. 
%     morseregion  - Generalized Morse wavelet time-frequency concentration region. 
%     isomax       - Returns those transform maxima that are isolated from others. 
%     maxprops     - Returns properties of wavelet transform maxima. 
%     max2eddy     - Converts transform maxima into oceanic coherent eddy properties. 
%     
%   Low-level functions 
%     morseafun  - Returns the generalized Morse wavelet amplitude or a-function.                  
%     morsehigh  - High-frequency cutoff of the generalized Morse wavelets. 
%     morseproj  - Projection coefficient for two generalized Morse wavelets.            
%     morsemom   - Frequency-domain moments of generalized Morse wavelets.  
%     jdawson    - The Dawson function. [By P.J. Acklam]           
%   
%   Morlet wavelet  
%     morlwave   - Morlet wavelet.                                                       
%     morlfreq   - Compute Morlet wavelet carrier frequency given peak frequency.    
%   
%   See also jRidges, jSpectral, jEllipse. 
%   __________________________________________________________________ 
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015--2024 J.M. Lilly --- type 'help jlab_license' for details
 
if nargin==0
    help jlab_allhelp
else
    if strcmpi(varargin{1},'--create')
        jlab_allhelp_create;
        return
    end
end

function[]=jlab_allhelp_create

jlab_path=findpath('jlab');

dirlist=dir(jlab_path);
xx=[];

newlinestr=[ '  __________________________________________________________________' char(13) char(13)];
for i=1:length(dirlist)
    if dirlist(i).isdir
        if ~strcmpi(dirlist(i).name(1),'.')
            if ~strcmpi(dirlist(i).name,'doc')&&~strcmpi(dirlist(i).name,'html')&&~strcmpi(dirlist(i).name,'figures')
                 eval(['xx=[xx help(''' dirlist(i).name ''') newlinestr];'])
            end
        end
    end
end
xx=real(xx);
xx(xx==10)=nan;
xx(xx==13)=nan;
xx=xx';
xx=col2mat(xx);
xx(4:end+3,:)=xx;
xx(1,:)=10;
xx(2,:)=real('%');
xx(3,:)=real(' ');
xx=setstr(mat2col(xx))';
xx=xx(1:end-3);
disp(xx)
