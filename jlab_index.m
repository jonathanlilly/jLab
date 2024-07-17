function[]=jlab_index(varargin)
%JLAB_INDEX  Alphabetical index into JLAB contents.
%
%   ab2kl         - Converts A and B to ellipse parameters Kappa and Lambda.                             
%   about_jtopo   - One-sixth degree global topography, from Smith and Sandwell + IBCAO.                 
%   allall        - ALLALL(X)=ALL(X(:))                                                                  
%   anatrans      - Analytic part of signal.                                                             
%   anyany        - ANYANY(X)=ANY(X(:))                                                                  
%   aresame       - Test whether two N-D arrays are the same.                                            
%   arrayify      - Converts a set of scalars or arrays into column arrays.                              
%   axeshandles   - Returns handles to all axes children.                                                
%   bellpoly      - Complete Bell polynomials.                                                           
%   besselitilde  - I-type Bessel function after factoring off exponential growth.                       
%   besselktilde  - K-type Bessel function after factoring off exponential decay.                        
%   bindata       - Rapidly sort data into adjacent bins.                                                
%   blocklen      - Counts the lengths of 'blocks' in an array.                                          
%   blocknum      - Numbers the contiguous blocks of an array.                                           
%   blurspec      - Returns the blurred and aliased spectrum given the autocovariance.                   
%   boxon         - Sets 'box' property to 'off'.                                                        
%   boxon         - Sets 'box' property to 'on'.                                                         
%   catstruct     - Concatenates the array elements of a cell array of structures.                       
%   cell2col      - Converts cell arrays of numeric arrays into 'column-appended' form.                  
%   cellabs       - Absolute value of each element in a cell array.                                      
%   celladd       - Addition acting on each element in a cell array.                                     
%   cellchunk     - Converts cell array data into uniform length 'chunks'.                               
%   celldiv       - Division acting on each element in a cell array.                                     
%   celldot       - Dot product for arrays of column vectors.                                            
%   cellength     - Length of each element in a cell array.                                              
%   cellfill      - Fills missing data marked by INFs in a cell array.                                   
%   cellfirst     - Returns the first (or last) element of each entry in a cell array.                   
%   cellget       - Indexes a cell array of numerical arrays by ID number.                               
%   cellgrid      - Interpolate a cell array of numeric arrays onto a regular grid.                      
%   cellimag      - Imaginary part of each element in a cell array.                                      
%   cellimit      - Limits the ranges of times in a cell array of numerical arrays.                      
%   cellindex     - Applies a cell array of indices to a cell array of column vectors.                   
%   celllog10     - Base ten logarithm of each element in a cell array.                                  
%   cellmax       - Maximum of each element in a cell array.                                             
%   cellmean      - Mean value of each element a cell array, possibly weighted.                          
%   cellmed       - Median value of each element a cell array.                                           
%   cellmin       - Minimum of each element in a cell array.                                             
%   cellmult      - Multiplication acting on each element in a cell array.                               
%   cellpack      - Removes all INF values from cell arrays of numeric arrays.                           
%   cellpair      - Complex pairing for elements in a cell array.                                        
%   cellplot      - Rapidly plot all elements of a cell array of numeric arrays.                         
%   cellprune     - Removes all empty cells, or cells less than a specified length.                      
%   cellreal      - Real part of each element in a cell array.                                           
%   cellsize      - Size of each element in a cell array along specified dimension.                      
%   cellsplit     - Splits cell arrays of numeric arrays at data gaps.                                   
%   cellstd       - Standard deviation of each element a cell array, possibly weighted.                  
%   cellstrip     - Strips INF values from the beginnings or ends of cell arrays.                        
%   cellsum       - Sum of each element a cell array, possibly weighted.                                 
%   chisquared    - The chi-squared distribution.                                                        
%   choose        - Binomial coefficient: CHOOSE(N,K) = N!K!/(N-K)!                                      
%   closedcurves  - Locate and interpolate closed curves in a possibly periodic domain.                  
%   cms2kmd       - Converts centimeters per second to kilometers per day.                               
%   col2cell      - Converts 'column-appended' data into cell arrays of numeric arrays.                  
%   col2mat       - Expands 'column-appended' data into a matrix.                                        
%   colbreaks     - Insert NANs into discontinuties in a vector.                                         
%   colorquant    - Set the color axis according to quantiles of the data.                               
%   commentlines  - Returns the comment lines from m-files.                                              
%   corfreq       - Coriolis frequency in radians per hour.                                              
%   crop          - Gets rid of whitespace around an image. [by A. Bliss]                                
%   crop_license  - License statement for CROP by Andrew Bliss.                                          
%   cum2mom       - Convert cumulants to moments.                                                        
%   curveinterp   - Interpolate a field or its gradient onto a set of curves.                            
%   curvemoments  - Centroid, area, and many other moments of a closed curve.                            
%   deg180        - Converts degrees to the range [-180,180].                                            
%   deg360        - Converts degrees to the range [0, 360].                                              
%   degunwrap     - Unwraps arrays given in degrees.                                                     
%   discretecolorbar - Plots a colorbar with discrete variation.                                         
%   divgeom       - Geometric decomposition of eddy vorticity flux divergence.                           
%   dlines        - Add diagonal lines to a plot.                                                        
%   doublen       - Interpolates a time series to double its length.                                     
%   ecconv        - Converts between eccentricity measures.                                              
%   eddyfit2d     - Least squares fit of 2D velocity data to an eddy profile.                            
%   eddylevels    - Eddy ridge significance levels using the survival function.                          
%   eddyridges    - Coherent eddy ridges from Lagrangian trajectories.                                   
%   ellband       - Bandwidth of modulated elliptical signals in two or three dimensions.                
%   ellcurves     - Returns curves corresponding to specified ellipse properties.                        
%   elldiff       - Differentiation of modulated elliptical signals.                                     
%   ellipseplot   - Plot ellipses.                                                                       
%   ellparams     - Ellipse parameters of a modulated bivariate or trivariate oscillation.               
%   ellpol        - Polarization parameters of an elliptical signal.                                     
%   ellrad        - Average and instantaneous ellipse radius.                                            
%   ellrossby     - Ellipse Rossby number, for oceanographic applications.                               
%   ellsig        - Creates a modulated elliptical signal in two or three dimensions.                    
%   ellvel        - Average and instantaneous ellipse velocities.                                        
%   fillbad       - Linearly interpolate over bad data points.                                           
%   findfiles     - Returns all files in a directory with a specified extension.                         
%   findpath      - Returns the full pathname of a directory on the Matlab search path.                  
%   fixlabels     - Specify precision of axes labels.                                                    
%   flipmap       - Flips the current colormap upside-down.                                              
%   flipx         - Flips the direction of the x-axis                                                    
%   flipy         - Flips the direction of the y-axis                                                    
%   fminsearchbnd - FMINSEARCH, but with bound constraints by transformation. [By J. D'Errico]           
%   fontsize      - Rapidly set title, axes, label, and text fontsizes.                                  
%   fourier       - Returns the Fourier frequencies for a given length time series.                      
%   frac          - Fraction: FRAC(A,B)=A./B                                                             
%   griddrifters  - Average drifter velocities onto a space/time 3D grid.                                
%   gulf4plot     - A four-panel circulation plot for the Gulf of Mexico.                                
%   haxby         - The Haxby colormap.                                                                  
%   haxby_copyright - Copyright statement for the HAXBY colormap.                                        
%   hlines        - Add horizontal lines to a plot.                                                      
%   hodograph     - Generate hodograph plots (simple and fancy).                                         
%   imlog         - Imaginary part of log: IMLOG(X)=UNWRAP(IMAG(LOG(X)))                                 
%   inellipse     - Locates points on the interior of ellipses.                                          
%   inregion      - Tests whether lat/lon points lie within a specified box.                             
%   instmom       - Univariate and multivariate instantaneous moments.                                   
%   interplatlon  - Interpolation for working with latitude and longitude.                               
%   inticks       - Sets the 'tickdir' property of the current axis to 'in'.                             
%   iseven        - True for even integer values; false otherwise.                                       
%   isodd         - True for odd integer values; false otherwise.                                        
%   isomax        - Returns those transform maxima that are isolated from others.                        
%   isridgepoint  - Finds wavelet ridge points using one of several criterion.                           
%   jdawson       - The Dawson function and its derivatives. [With P.J. Acklam]                          
%   jdeg2rad      - Converts degrees to radians.                                                         
%   jfig          - Shorthand for tweaking figures.                                                      
%   jhanning      - Hanning window.                                                                      
%   jhelp         - Opens linked JLAB help files in Matlab's internal web browser.                       
%   jhermfun      - Orthonormal Hermite functions. [with F. Rekibi]                                      
%   jhermpoly     - Hermite polynomials. [with F. Rekibi]                                                
%   jlab_addpath  - Adds JLAB and JDATA subdirectories to your Matlab search path.                       
%   jcell:        - Tools for operating on cell arrays of column vectors                                 
%   jlab_changes  - Changes to JLAB in each release.                                                     
%   jlab_highlights - Introduction to some of the most useful routines.                                  
%   jlab_index    - Alphabetical index into JLAB and JDATA contents.                                     
%   jlab_install  - Instructions for installing JLAB.                                                    
%   jlab_license  - License statement and permissions for JLAB package.                                  
%   jlab_makefigs - Make figures for papers by J. M. Lilly.                                              
%   jlab_runtests - Runs a test suite for JLAB package.                                                  
%   jlab_thanks   - Sources of support for developing this software.                                     
%   jmat2         - 2x2 rotation matrix through specified angle.                                         
%   jmat3         - 3x3 rotation matrix through specified angle.                                         
%   jpcolor       - Modified version of PCOLOR appropriate for cell-centered grids.                      
%   jprint        - Print to a specified directory and crop the resulting file.                          
%   jrad2deg      - Converts radians to degrees.                                                         
%   kl2ab         - Converts ellipse parameters Kappa and Lambda to A and B.                             
%   lansey        - The Lansey modification of Cynthia Brewer's "Spectral" colormap.                     
%   lansey_copyright - Copyright statement for the Lansey colormap.                                      
%   latlon2uv     - Converts latitude and longitude to horizontal velocity.                              
%   latlon2xy     - Converts latitude and longitude into tangent plane coordinates.                      
%   latlon2xyz    - Converts latitude and longitude into 3D Cartesian coordinates.                       
%   latratio      - Set plot aspect ratio for latitude / longitude plot.                                 
%   letterlabels  - For automatically putting letter labels on subplots.                                 
%   linecolor     - Set line colors based on a property value within a colormap.                         
%   linehandles   - Finds all line and patch handles from a given set of axes.                           
%   linering      - Moves lines through the current line style order.                                    
%   linestyle     - Rapidly set color, style, and width properties of lines.                             
%   linestyleparse - Parses the input string to LINESTYLE.                                               
%   lininterp     - Fast linear interpolation for arbitrary-sized arrays.                                
%   lnsd          - Last non-singleton dimension of an array.                                            
%   lonshift      - Shifts longitude origin for plotting purposes.                                       
%   make          - Create a structure containing named variables as fields.                             
%   make_betaeddyone - Create the BetaEddyOne nonlinear QG eddy simulation.                              
%   make_goldaltimetry - Create a synthetic altimeter dataset and gridded statistics.                    
%   make_goldocean - Create quarter-degree gridded files from GOLD simulations.                          
%   make_gomed    - Create a drifter-derived eddy database for the Gulf of Mexico.                       
%   make_gulfdrifters - Create a drifter dataset for the Gulf of Mexico.                                 
%   make_gulfflow - Create a gridded velocity dataset for the Gulf of Mexico.                            
%   make_jason    - Create a reformatted version of the Beckley altimeter dataset.                       
%   makefigs_analytic - Make figures for Lilly and Olhede (2010b).                                       
%   makefigs_asilomar - Make figures for Lilly and Olhede (2009b).                                       
%   makefigs_bandwidth - Make figures for Lilly and Olhede (2010a).                                      
%   makefigs_closedcurves - Makes a sample figure for CLOSEDCURVES.                                      
%   makefigs_curvemoments - Makes a sample figure for CURVEMOMENTS.                                      
%   makefigs_dlines - Makes a sample figure for DLINES.                                                  
%   makefigs_doublen - Makes a sample figure for DOUBLEN.                                                
%   makefigs_ecconv - Makes a sample figure for ECCONV.                                                  
%   makefigs_element - Make figures for Lilly (2017).                                                    
%   makefigs_ellband - Makes a sample figure for ELLBAND.                                                
%   makefigs_ellipseplot - Makes a sample figure for ELLIPSEPLOT.                                        
%   makefigs_ellvel - Makes a sample figure for ELLVEL.                                                  
%   makefigs_gulfcensus - Make figures for Lilly and Perez-Brunius (2021b).                              
%   makefigs_gulfdrifters - Make figures for Lilly and Perez-Brunius (2021a).                            
%   makefigs_inellipse - Makes a sample figure for INELLIPSE.                                            
%   makefigs_instmom - Makes some sample figures for INSTMOM.                                            
%   makefigs_jdawson - Makes a sample figure for DJAWSON.                                                
%   makefigs_jfig - Makes a sample figure for JFIG.                                                      
%   makefigs_jhermfun - Makes a sample figure for JHERMFUN.                                              
%   makefigs_jpcolor - Makes a sample figure for JPCOLOR.                                                
%   makefigs_jtopo - Makes a sample figure for ABOUT_JTOPO.                                              
%   makefigs_kinematics - Make figures for Lilly (2018).                                                 
%   makefigs_lansey - Makes a sample figure for LANSEY.                                                  
%   makefigs_matern - Make all figures for Lilly et al. 2017.                                            
%   makefigs_matern - Make all figures for Lilly et al. 2017.                                            
%   makefigs_maternoise - Makes a sample figure for MATERNOISE.                                          
%   makefigs_maternspec - Makes some sample figures for MATERNSPEC.                                      
%   makefigs_mconf - Makes sample figure for MCONF.                                                      
%   makefigs_morlfreq - Makes a sample figure for MORLFREQ.                                              
%   makefigs_morlwave - Makes a sample figure for MORLWAVE.                                              
%   makefigs_morsearea - Makes a sample figure for MORSEAREA.                                            
%   makefigs_morsebox - Makes a sample figure for MORSEBOX.                                              
%   makefigs_morsefreq - Makes a sample figure for MORSEFREQ.                                            
%   makefigs_morseregion - Makes somes sample figures for MORSEREGION.                                   
%   makefigs_morsespace - Makes some sample figures for MORSESPACE.                                      
%   makefigs_morsewave - Makes some sample figures for MORSEWAVE.                                        
%   makefigs_morsies - Make figures for Lilly and Olhede (2009a).                                        
%   makefigs_mspec - Makes some sample figures for MSPEC.                                                
%   makefigs_multivariate - Make figures for Lilly and Olhede (2012a).                                   
%   makefigs_patchcontourf - Makes a sample figure for PATCHCONTOURF.                                    
%   makefigs_periodize - Makes a sample figure for PERIODIZE.                                            
%   makefigs_polymap - Makes some sample figures for POLYMAP.                                            
%   makefigs_polymap2 - Makes a sample figure for POLYMAP using TPJAOS.MAT.                              
%   makefigs_psi2fields - Makes a sample figure for PSI2FIELDS.                                          
%   makefigs_regionplot - Makes some sample figures for REGIONPLOT.                                      
%   makefigs_ridges - Make figures for Lilly, Scott, and Olhede (2011).                                  
%   makefigs_ridgetrim - Makes a sample figure for RIDGETRIM.                                            
%   makefigs_ridgewalk - Makes a sample figure for RIDGEWALK.                                            
%   makefigs_seminfhaxby - Makes a sample figure for SEMINFHAXBY.                                        
%   makefigs_simplepdf - Makes a sample figure for SIMPLEPDF.                                            
%   makefigs_sphereinterp - Makes two sample figures for SPHEREINTERP.                                   
%   makefigs_stickvect - Makes a sample figure for STICKVECT.                                            
%   makefigs_stokes - Make figures for Lilly et al. (2024).                                              
%   makefigs_superfamily - Make figures for Lilly and Olhede (2012b).                                    
%   makefigs_topoplot - Makes a sample figure for TOPOPLOT.                                              
%   makefigs_trackextract - Makes a sample figure for TRACKEXTRACT.                                      
%   makefigs_trajextract - Makes a sample figure for TRAJEXTRACT.                                        
%   makefigs_trajwrap - Makes a sample figure for TRAJWRAP.                                              
%   makefigs_transfer - Make figures for Lilly and Elipot (2021).                                        
%   makefigs_trivariate - Make figures for Lilly (2011).                                                 
%   makefigs_twodhist - Makes a sample figure for TWODHIST.                                              
%   makefigs_twodmed - Makes a sample figure for TWODMED.                                                
%   makefigs_twodstats - Makes a sample figure for TWODSTATS.                                            
%   makefigs_vortex - Make figures for Lilly and Gascard (2006).                                         
%   makefigs_wavespecplot - Makes a sample figure for WAVESPECPLOT.                                      
%   makefigs_wavetrans - Makes a sample figure for WAVETRANS.                                            
%   makefigs_widgist - Makes some sample figures for WIGDIST.                                            
%   mat2col       - Compress NAN-padded matrix data into long columns.                                   
%   materncfun    - Returns the normalization or C-function for a Matern process.                        
%   maternchol    - Cholesky decomposition of Matern and fBm covariances. [with A. Sykulski]             
%   materncov     - Autocovariance of the Matern random process and variations.                          
%   maternedge    - Long-time cutoff edge for the Matern impulse response function.                      
%   maternfit     - Parametric spectral fit to the Matern form. [with A. Sykulski]                       
%   maternfun     - Returns the Matern function.                                                         
%   maternimp     - Impulse response function for the Matern random process.                             
%   maternoise    - Realizations of the Matern process and variations, including fBm. [with A. Sykulski] 
%   maternspec    - Fourier spectrum of the Matern random process and variations.                        
%   matsave       - Create and save structure of variables as a mat-file.                                
%   max2eddy      - Converts transform maxima into oceanic coherent eddy properties.                     
%   maxmax        - MAXMAX(X)=MAX(X(~ISNAN(X(:))))                                                       
%   maxprops      - Returns properties of wavelet transform maxima.                                      
%   mconf         - Confidence intervals for the multitaper spectral estimate.                           
%   minmin        - MINMIN(X)=MIN(X(~ISNAN(X(:))))                                                       
%   mom2cum       - Convert moments to cumulants.                                                        
%   monthstats    - Mean month and standard deviation using circular statistics.                         
%   morlfreq      - Compute Morlet wavelet carrier frequency given peak frequency.                       
%   morlwave      - Morlet wavelet.                                                                      
%   morseafun     - Returns the generalized Morse wavelet amplitude or a-function.                       
%   morsearea     - Time-frequency concentration area of Morse wavelets. [with F. Rekibi]                
%   morsebox      - Heisenberg time-frequency box for generalized Morse wavelets.                        
%   morsederiv    - Frequency-domain derivatives of generalized Morse wavelets.                          
%   morsefreq     - Frequency measures for generalized Morse wavelets. [with F. Rekibi]                  
%   morsehigh     - High-frequency cutoff of the generalized Morse wavelets.                             
%   morsemom      - Frequency-domain moments of generalized Morse wavelets.                              
%   morseproj     - Projection coefficient for two generalized Morse wavelets.                           
%   morseprops    - Properties of the demodulated generalized Morse wavelets.                            
%   morseregion   - Generalized Morse wavelet time-frequency concentration region.                       
%   morsespace    - Logarithmically-spaced frequencies for generalized Morse wavelets.                   
%   morsewave     - Generalized Morse wavelets of Olhede and Walden (2002).                              
%   morsexpand    - Generalized Morse wavelets via time-domain Taylor series.                            
%   mspec         - Multitaper power and cross spectra.                                                  
%   msvd          - Singular value decomposition for polarization analysis.                              
%   ncinterp      - One-line interpolation from 3D lat/lon/time field in NetCDF file.                    
%   ncload        - Load all variables from a NetCDF file and convert trajectories to cells.             
%   nocontours    - Removes contours from a CONTOURF plot.                                               
%   noisedrifters - Create a noise Lagrangian dataset matching mean and variance.                        
%   nonnan        - Return all non-NAN elements of an array.                                             
%   noxlabels     - Remove some or all x-axis tick mark labels.                                          
%   noylabels     - Remove some or all y-axis tick mark labels.                                          
%   oprod         - Outer product:  OPROD(X,Y)=X*CONJ(Y')                                                
%   orbitbreaks   - Separate orbit into passes based on turning points.                                  
%   outticks      - Sets the 'tickdir' property of the current axis to 'out'.                            
%   packfig       - Squeeze together rows and/or columns of the current figure.                          
%   patchcontourf - Generate filled contours using patches, with specified colors.                       
%   pdfprops      - Mean and variance associated with a probability distribution.                        
%   periodindex   - Returns time index in increments of instantaneous period.                            
%   periodize     - Returns a doubly periodic version of an input array.                                 
%   polparams     - Spectral matrix polarization parameters.                                             
%   printall      - Print and close all open figures.                                                    
%   provec        - Generate progressive vector diagrams (simple and fancy).                             
%   psi2fields    - Velocity and other fields from the streamfunction. [with P.E. Isachsen]              
%   quadinterp    - Fast quadratic interpolation for arbitrary-sized arrays.                             
%   radearth      - The radius of the earth in kilometers.                                               
%   regionplot    - Plots a simple box indicating a latitude / longitude region.                         
%   reporttest    - Reports the result of an m-file function auto-test.                                  
%   res           - Residual after flooring:  RES(X)=X-FLOOR(X)                                          
%   ridgechains   - Forms ridge curves by connecting transform ridge points.                             
%   ridgeinterp   - Interpolate quantity values onto ridge locations.                                    
%   ridgelen      - Wavelet ridge length expressed as number of full cycles.                             
%   ridgemap      - Maps ridge quantities back onto the time series.                                     
%   ridgemult     - Ridge multiplicity, the number of simultaneous ridges present.                       
%   ridgetrim     - Trim edge effect regions from wavelet ridges.                                        
%   ridgewalk     - Extract wavelet transform ridges, including bias estimates.                          
%   rot           - Complex-valued rotation:  ROT(X)=EXP(SQRT(-1)*X)                                     
%   sampletimes   - Computes mean sampling intervals and their statistics.                               
%   seminfhaxby   - The seminf-Haxby colormap.                                                           
%   sig2latlon    - Converts an oscillatory signal to lat/lon displacements.                             
%   simpleddy     - Streamfunction, velocity, and vorticity for various eddy profiles.                   
%   simplepdf     - Gaussian, uniform, Cauchy, and exponential pdfs.                                     
%   sleptap       - Calculate Slepian tapers.                                                            
%   specdiag      - Diagonalize a 2 x 2 spectral matrix.                                                 
%   sphere2uvw    - Converts a 3D spherical vector to a 3D Cartesian vector.                             
%   spherecurl    - Curl of a vector field on the surface of a sphere.                                   
%   spheredist    - Computes great circle distances on a sphere.                                         
%   spherediv     - Divergence of a vector field on the surface of a sphere.                             
%   spheregrad    - Gradient of a field on the surface of a sphere.                                      
%   sphereinterp  - Fast linear interpolation on the sphere from non-plaid grids.                        
%   spherelap     - Laplacian of a field on the surface of a sphere.                                     
%   spheretrans   - Wavelet transform for oscillations on the surface of a sphere.                       
%   squared       - Squares the modulus of its argument:  SQUARED(X)=ABS(X).^2                           
%   standalone    - Create stand-alone version of an m-file, including dependencies.                     
%   stickvect     - Plots "stick vectors" for multicomponent velocity time series.                       
%   structindex   - Applies an index to all array-valued fields in a structure.                          
%   tidefreq      - Frequencies of the eight major tidal components.                                     
%   timeseries_boundary - Apply boundary conditions to data before transform.                            
%   tmat          - 2x2 complex grouping matrix.  TMAT = [1  i; 1 -i] / SQRT(2)                          
%   to_grab_from_caller - Returns a string to grab variable values from caller.                          
%   to_overwrite  - Returns a string to overwrite original arguments.                                    
%   topo_copyright - Copyright statement for the Smith and Sandwell topography.                          
%   topoplot      - Plot regional or global topography at one-sixth degree resolution.                   
%   trackextract  - Extracts alongtrack altimetry segments within given region.                          
%   trajchunk     - Chunks Lagrangian trajectories based on the Coriolis period.                         
%   trajextract   - Extracts Lagrangian trajectory segments within given region.                         
%   trajunwrap    - Unwraps Lagrangian trajectories from a periodic domain.                              
%   trajwrap      - Wraps Lagrangian trajectories to fit within a periodic domain.                       
%   transmax      - Locates modulus maximum points of wavelet transform.                                 
%   transmaxdist  - Distributions of wavelet transform maxima in noise.                                  
%   twodhist      - Two-dimensional histogram.                                                           
%   twodmed       - Median value of a function of two variables.                                         
%   twodstats     - Mean, variance, and covariance of functions of two variables.                        
%   twospecplot   - Plots a pair of rotary or Cartesian spectra.                                         
%   use           - Copies structure fields into named variables in workspace.                           
%   uv2latlon     - Integrates horizontal velocity to give latitude and longitude.                       
%   uvplot        - Plots the real and imaginary parts of a signal on the same axis.                     
%   vcolon        - Condenses its arguments, like X(:).                                                  
%   vdiff         - Length-preserving first central difference.                                          
%   vectmult      - Matrix multiplication for arrays of vectors.                                         
%   vempty        - Initializes multiple variables to empty sets or empty cell arrays.                   
%   vfilt         - Filtering along rows without change in length.                                       
%   vindex        - Indexes an N-D array along a specified dimension.                                    
%   vindexinto    - Indexes into N-D array along a specified dimension.                                  
%   vlines        - Add vertical lines to a plot.                                                        
%   vmean         - Mean over non-NaN elements along a specified dimension.                              
%   vmedian       - Median over finite elements along a specified dimension.                             
%   vmoment       - Central moment over non-NaN elements along a specfied dimension.                     
%   vrep          - Replicates an array along a specified dimension.                                     
%   vrms          - Root-mean-square of non-NaN elements along a specified dimension.                    
%   vshift        - Cycles the elements of an array along a specified dimension.                         
%   vsize         - Returns the sizes of multiple arguments.                                             
%   vsqueeze      - Squeezes multiple input arguments simultaneously.                                    
%   vstd          - Standard deviation over non-NaN elements along a specfied dimension.                 
%   vsum          - Sum over non-NaN elements along a specified dimension.                               
%   vswap         - VSWAP(X,A,B) replaces A with B in numeric array X.                                   
%   vtranspose    - Transpose multiple input arguments simultaneously.                                   
%   vzeros        - Initializes multiple variables to arrays of zeros or nans.                           
%   wavespecplot  - Plot of wavelet spectra together with time series.                                   
%   wavetrans     - Continuous wavelet transform.                                                        
%   whichdir      - Returns directory name containing file in search path.                               
%   wigdist       - Wigner distribution (alias-free algorithm).                                          
%   windtrans     - Ekman-like transfer-functions for the wind-driven response.                          
%   xlin          - Sets x-axis scale to linear.                                                         
%   xlog          - Sets x-axis scale to logarithmic.                                                    
%   xoffset       - Offsets lines in the x-direction after plotting.                                     
%   xtick         - Sets locations of x-axis tick marks.                                                 
%   xy2latlon     - Converts tangent plane coordinates into latitude and longitude.                      
%   xyz2latlon    - Converts 3D Cartesian coordinates into latitude and longitude.                       
%   yearfrac      - Converts a DATENUM into 'year.fraction' and 'month.fraction'.                        
%   ylin          - Sets y-axis scale to linear.                                                         
%   ylog          - Sets y-axis scale to log.                                                            
%   yoffset       - Offsets lines in the y-direction after plotting.                                     
%   ytick         - Sets locations of y-axis tick marks.                                                 
%   ztick         - Sets locations of z-axis tick marks.     
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2024 J.M. Lilly --- type 'help jlab_license' for details


if nargin==0
    help jlab_index
else
    if strcmpi(varargin{1},'--create')
        jlab_index_create;
        return
    end
end

function[]=jlab_index_create
jlab_dir=whichdir('jlab_license');
jdata_dir=[jlab_dir(1:end-4) 'jData'];


for k=1
    switch k
        case 1
            mfiles=findfiles(jlab_dir,'m','recursive');
        case 2
            mfiles=findfiles(jdata_dir,'m','recursive');
    end
    for i=1:length(mfiles)
        if strcmp(mfiles{i},'Contents.m')
            mfiles{i}=[];
        end
    end
    mfiles=mfiles(cellength(mfiles)>0);
    mfiles=strs2mat(mfiles);
    
    for i=size(mfiles,2):-1:1
        [temp,sorter]=sort(mfiles(:,i));
        mfiles=mfiles(sorter,:);
    end
    clear comments
    for i=1:size(mfiles,1)
        %mfiles(i,:)
        comments{i}=commentlines(mfiles(i,:));
    end
    switch k 
        case 1
            comments1=strs2mat(comments);
        case 2
            comments2=strs2mat(comments);
    end
end


% clear comments
% comments{1}=comments2;
% comments{2}=comments1;

% disp('Start here -->')
% disp('%JLAB_INDEX  Alphabetical index into JLAB and JDATA contents.')
% disp('%')
% disp('% JDATA index')
% disp(comments2)
% disp('%')
% disp('% JLAB index')
% disp(comments1)

disp('Start here -->')
disp('%JLAB_INDEX  Alphabetical index into JLAB contents.')
disp('%')
disp(comments1)


function[mat]=strs2mat(strcell)

for i=1:length(strcell)
    strcell{i}=real(strcell{i}');
end

cell2col(strcell);
col2mat(strcell);
vswap(strcell,nan,32);
mat=char(strcell');



function[]=jlab_index_test
 
%reporttest('JLAB_INDEX',aresame())
