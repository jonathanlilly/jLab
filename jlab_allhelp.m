function[varargout]=jlab_allhelp(varargin)
%JLAB_ALLHELP  Displays the help comments for all JLAB modules. 
%
%   jCell:  Tools for operating on cell arrays of column vectors
%  
%   Basic mathematical operations
%     cellabs    - Absolute value of each element in a cell array.                    
%     cellmax    - Maximum of each element in a cell array.                           
%     cellmin    - Minimum of each element in a cell array. 
%     cellmean   - Mean value of each element a cell array.  
%     cellstd    - Standard deviation of each element a cell array.
%     cellreal   - Real part of each element in a cell array.                         
%     cellimag   - Imaginary part of each element in a cell array.    
%     celllog10  - Base ten logarithm of each element in a cell array.
%     celladd    - Addition acting on each element in a cell array.                   
%     cellmult   - Multiplication acting on each element in a cell array.   
%     celldiv    - Division acting on each element in a cell array.
%  
%   Reshaping, indexing, and sizes
%     cell2col    - Converts cell arrays of numeric arrays into 'column-appended' form.
%     col2cell    - Converts 'column-appended' data into cell arrays of numeric arrays.
%     cellindex  - Applies a cell array of indices to a cell array of column vectors. 
%     cellchunk  - Converts cell array data into uniform length 'chunks'.             
%     cellength  - Length of each element in a cell array.                            
%     cellsize   - Size of each element in a cell array along specified dimension.
%     cellget    - Indexes a cell array of numerical arrays by ID number.
%     cellimit   - Limits the ranges of times in a cell array of numerical arrays.
%  
%   Data processing
%     cellstrip  - Strips NaN values from the beginnings or ends of cell arrays.
%     cellsplit  - Splits cell arrays of numeric arrays at data gaps.
%     cellprune  - Removes all empty cells, or cells less than a specified length.
%     cellfill   - Fills missing data marked by NaNs in a cell array.
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
%     matmult    - Matrix multiplication for arrays of matrices.                      
%     nonnan     - Return all non-NAN elements of an array.                           
%     vectmult   - Matrix multiplication for arrays of vectors.    
%  
%   File and directory tools
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
%     hermpoly   - Hermite polynomials. [with F. Rekibi]   
%     hermfun    - Orthonormal Hermite functions. [with F. Rekibi]                    
%  
%   Dataset organization as structures
%     make       - Create a structure containing named variables as fields.  
%     matsave    - Create and save structure of variables as a mat-file. 
%     use        - Copies structure fields into named variables in workspace.         
%  
%   Statistics
%     cum2mom    - Convert cumulants to moments.    
%     mom2cum    - Convert moments to cumulants.                                      
%     pdfprops   - Mean and variance associated with a probability distribution.  
%     simplepdf  - Gaussian, uniform, Cauchy, and exponential pdfs.                   
%  
%   Filling bad data points
%     fillbad    - Linearly interpolate over bad data points.  
%  
%   Date, time, and units
%     yearfrac    - Converts a DATENUM into 'year.fraction' and 'month.fraction'.
%     cms2kmd     - Converts centimeters per second to kilometers per day.
%   __________________________________________________________________
% 
%   jEllipse:  Analysis of modulated elliptical, or bivariate, signals
%  
%    Elliptical signal properties
%     ellband    - Bandwidth of modulated elliptical signals in two or three dimensions. 
%     ellparams  - Ellipse parameters of a modulated bivariate or trivariate oscillation.
%     ellvel     - Average and instantaneous ellipse velocities.    
%     ellrad     - Average and instantaneous ellipse radius.                             
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
%   Graphics aliases
%     boxon      - Sets 'box' property to 'off'.                                    
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
%   Colormaps
%     lansey     - The Lansey modification of Cynthia Brewer's "Spectral" colormap.
%  
%   Low-level functions
%     axeshandles - Returns handles to all axes children.
%     crop        - Gets rid of whitespace around an image. [by A. Bliss]             
%     linehandles - Finds all line and patch handles from a given set of axes.
%     linestyleparse - Parses the input string to LINESTYLE.
%   __________________________________________________________________
% 
%   jMap:  Mapping scattered data using local polynomial fitting
%  
%     polysmooth - Smoothing scattered 2D data with local polynomial fitting. 
%     spheresort - Sorted great circle distances to nearby points on the earth.
%     twodsort   - Distances from data points to nearby grid points.   
%   __________________________________________________________________
% 
%   jMatern:  Parametric spectral analysis based on the Matern process
%  
%   Top-level functions
%     maternoise - Realizations of the Matern process and variations, including fBm.  [with A. Sykulski]
%     maternspec - Fourier spectrum of the Matern random process and variations.                                
%     materncov  - Autocovariance of the Matern random process and variations.                         
%     maternimp  - Impulse response function for the Matern random process.                      
%  
%   Other utilities
%     blurspec   - Returns the blurred and aliased spectrum given the autocovariance.
%     fminsearchbnd  - FMINSEARCH, but with bound constraints by transformation. [By J. D'Errico]
%  
%   Low-level Matern functions
%     materncfun - Returns the normalization function C_ALPHA for a Matern process.
%     maternchol - Cholesky decomposition of Matern and fBm covariances.
%     maternedge - Long-time cutoff edge for the Matern impulse response function.               
%  
%   See also jSpectral.
%   __________________________________________________________________
% 
%   jMatfiles:  Descriptions of included mat-files, mostly for testing and examples.
%  
%     bravo94      - Labrador Sea mooring data from Lilly et. al (1999).
%     drifterbetty - Two-hour resolution drifter trajectory, AOML GDP id #44000.
%     ebasnfloats  - "Meddy" float data analyzed in Lilly and Olhede (2009,10,11,12).
%     npg2006      - Float data analyzed in Lilly and Gascard (2006).
%     qgsnapshot   - Snapshot of f-plane QG turbulence from a simulation by J. Early.
%     qgmodelfit   - Parameter fits of Matern model to QG turbulence data.
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
%     trajextract - Extracts Lagrangian trajectory segments within given region.        
%     trajunwrap  - Unwraps Lagrangian trajectories from a periodic domain.              
%     trajwrap    - Wraps Lagrangian trajectories to fit within a periodic domain.       
%     trajchunk   - Converts cell array data into chunks based on the Coriolis period.
%  
%   Idealized numerical model tools
%     psi2fields   - Velocity and other fields from the streamfunction. [with P.E. Isachsen]    
%     periodize    - Returns a doubly periodic version of an input array.   
%  
%   Eulerian eddy identification and analysis
%     inellipse    - Locates points on the interior of ellipses.                          
%     closedcurves - Locate and interpolate closed curves in a possibly periodic domain.
%     curvemoments - Centroid, area, and many other moments of a closed curve.          
%     divgeom      - Geometric decomposition of eddy vorticity flux divergence.     
%  
%   Plotting tools for mooring data
%     hodograph  - Generate hodograph plots (simple and fancy).                                   
%     provec     - Generate progressive vector diagrams (simple and fancy).             
%     stickvect  - Plots "stick vectors" for multicomponent velocity time series. 
%  
%   Alongtrack altimetry tools
%     trackextract  - Extracts alongtrack altimetry segments within given region.
%  
%   NetCDF tools
%     ncinterp      - Interpolate field from NetCDF file onto specified positions.
%  
%   Topography tools and data
%     jtopo.mat   - One-sixth degree global topography, from Smith and Sandwell + IBCAO.              
%     topoplot    - Plot regional or global topography at one-sixth degree resolution.
%   
%    See also jData.
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
%     ridgelen   - Wavelet ridge length expressed as number of full cycles.    
%     periodindex - Returns time index in increments of instantaneous period.  
%  
%    See also jEllipse, jWavelet.
%   __________________________________________________________________
% 
%    jSpectral:  Multitaper spectral analysis, and other time series tools
%   
%    Multitaper spectral analysis
%     sleptap    - Calculate Slepian tapers.                                        
%     mspec      - Multitaper power and cross spectra.       
%  
%    Multitaper polarization analysis
%     msvd       - Singular value decomposition for polarization analysis.   
%     polparams  - Spectral matrix polarization parameters.                         
%     specdiag   - Diagonalize a 2 x 2 spectral matrix.        
%  
%   Assorted other transforms 
%     slidetrans  - Sliding-window ('moving-window') Fourier transform.   
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
%     latlon2xy  - Converts latitude and longitude into local Cartesian coordinates.
%     latlon2xyz - Converts latitude and longitude into 3D Cartesian coordinates.   
%     xy2latlon  - Converts local Cartesian coordinates into latitude and longitude.
%     xyz2latlon - Converts 3D Cartesian coordinates into latitude and longitude.   
%     sphere2uvw - Converts a 3D spherical vector to a 3D Cartesian vector. 
%  
%   Differential operators
%     spherecurl - Curl of a vector field on the surface of a sphere.               
%     spherediv  - Divergence of a vector field on the surface of a sphere.         
%     spheregrad - Gradient of a field on the surface of a sphere.                  
%     spherelap  - Laplacian of a field on the surface of a sphere.   
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
%    See also jCell.
%   __________________________________________________________________
% 
%   jWavelet: Continuous wavelet analysis using generalized Morse wavelets
%  
%   Top-level wavelet transform functions
%     wavetrans  - Continuous wavelet transform.
%     morsewave  - Generalized Morse wavelets of Olhede and Walden (2002).              
%     morsespace - Logarithmically-spaced frequencies for generalized Morse wavelets.   
%     wavespecplot - Plot of wavelet spectra together with time series.
%  
%   Details of generalized Morse wavelets 
%     morsearea  - Time-frequency concentration area of Morse wavelets. [with F. Rekibi]
%     morsebox   - Heisenberg time-frequency box for generalized Morse wavelets.        
%     morsefreq  - Frequency measures for generalized Morse wavelets. [with F. Rekibi]  
%     morseprops - Properties of the demodulated generalized Morse wavelets.            
%     morsexpand - Generalized Morse wavelets via time-domain Taylor series.
%     morsederiv - Frequency-domain derivatives of generalized Morse wavelets.
%  
%   Morlet wavelet 
%     morlwave   - Morlet wavelet.                                                      
%     morlfreq   - Compute Morlet wavelet carrier frequency given peak frequency.   
%  
%   See also jRidges, jSpectral, jEllipse.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015--2016 J.M. Lilly --- type 'help jlab_license' for details
 
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
newline=[ '  __________________________________________________________________' char(13) char(13)];
for i=1:length(dirlist)
    if dirlist(i).isdir
        if ~strcmpi(dirlist(i).name(1),'.')
            if ~strcmpi(dirlist(i).name,'doc')&&~strcmpi(dirlist(i).name,'html')&&~strcmpi(dirlist(i).name,'figures')
                 eval(['xx=[xx help(''' dirlist(i).name ''') newline];'])
            end
        end
    end
end
xx
varargout{1}=xx;