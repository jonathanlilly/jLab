function[]=jlab_index(varargin)
%JLAB_INDEX  Alphabetical index into JLAB and JDATA contents.
%
%JLAB_INDEX  Alphabetical index into JLAB and JDATA contents.
%
% JDATA index
%   about_drifters - Global surface drifter dataset from the Global Drifter Program. 
%   about_floats   - Historical dataset of eddy-resolving subsurface floats.          
%   about_ibcao    - International Bathymetric Chart of the Arctic Ocean topography.  
%   about_sandwell - One minute resolution topography data from Smith and Sandwell.  
%   about_tpjaos   - Sea surface height anomalies from the Beckley merged dataset.    
%   readtopo       - Read one-minute topography data from Smith and Sandwell.         
%   topo_copyright - Copyright statement for the Smith and Sandwell topography.      
%
% JLAB index
%   ab2kl         - Converts A and B to ellipse parameters Kappa and Lambda.                              
%   about_jtopo   - One-sixth degree global topography, from Smith and Sandwell + IBCAO.                  
%   allall        - ALLALL(X)=ALL(X(:))                                                                   
%   anatrans      - Analytic part of signal.                                                              
%   anyany        - ANYANY(X)=ANY(X(:))                                                                   
%   aresame       - Test whether two N-D arrays are the same.                                             
%   arrayify      - Converts a set of scalars or arrays into column arrays.                               
%   axeshandles   - Returns handles to all axes children.                                                 
%   bellpoly      - Complete Bell polynomials.                                                            
%   bindata       - Rapidly sort data into adjacent bins.                                                 
%   blocklen      - Counts the lengths of 'blocks' in an array.                                           
%   blocknum      - Numbers the contiguous blocks of an array.                                            
%   blurspec      - Returns the blurred and aliased spectrum given the autocovariance.                    
%   boxon         - Sets 'box' property to 'off'.                                                         
%   boxon         - Sets 'box' property to 'on'.                                                          
%   cell2col      - Converts cell arrays of numeric arrays into 'column-appended' form.                   
%   cellabs       - Absolute value of each element in a cell array.                                       
%   celladd       - Addition acting on each element in a cell array.                                      
%   cellchunk     - Converts cell array data into uniform length 'chunks'.                                
%   celldiv       - Division acting on each element in a cell array.                                      
%   cellength     - Length of each element in a cell array.                                               
%   cellfill      - Fills missing data marked by NaNs in a cell array.                                    
%   cellfirst     - Returns the first (or last) element of each entry in a cell array.                    
%   cellget       - Indexes a cell array of numerical arrays by ID number.                                
%   cellgrid      - Interpolate a cell array of numeric arrays onto a regular grid.                       
%   cellimag      - Imaginary part of each element in a cell array.                                       
%   cellimit      - Limits the ranges of times in a cell array of numerical arrays.                       
%   cellindex     - Applies a cell array of indices to a cell array of column vectors.                    
%   celllog10     - Base ten logarithm of each element in a cell array.                                   
%   cellmax       - Maximum of each element in a cell array.                                              
%   cellmean      - Mean value of each element a cell array.                                              
%   cellmin       - Minimum of each element in a cell array.                                              
%   cellmult      - Multiplication acting on each element in a cell array.                                
%   cellplot      - Rapidly plot all elements of a cell array of numeric arrays.                          
%   cellprune     - Removes all empty cells, or cells less than a specified length.                       
%   cellreal      - Real part of each element in a cell array.                                            
%   cellsize      - Size of each element in a cell array along specified dimension.                       
%   cellsplit     - Splits cell arrays of numeric arrays at data gaps.                                    
%   cellstd       - Standard deviation of each element a cell array.                                      
%   cellstrip     - Strips NaN values from the beginnings or ends of cell arrays.                         
%   choose        - Binomial coefficient: CHOOSE(N,K) = N!K!/(N-K)!                                       
%   closedcurves  - Locate and interpolate closed curves in a possibly periodic domain.                   
%   cms2kmd       - Converts centimeters per second to kilometers per day.                                
%   col2cell      - Converts 'column-appended' data into cell arrays of numeric arrays.                   
%   col2mat       - Expands 'column-appended' data into a matrix.                                         
%   colbreaks     - Insert NANs into discontinuties in a vector.                                          
%   commentlines  - Returns the comment lines from m-files.                                               
%   corfreq       - Coriolis frequency in radians per hour.                                               
%   crop          - Gets rid of whitespace around an image. [by A. Bliss]                                 
%   crop_license  - License statement for CROP by Andrew Bliss                                            
%   cum2mom       - Convert cumulants to moments.                                                         
%   curveinterp   - Interpolate a field or its gradient onto a set of curves.                             
%   curvemoments  - Centroid, area, and many other moments of a closed curve.                             
%   dawson        - The Dawson function and its derivatives. [With P.J. Acklam]                           
%   deg180        - Converts degrees to the range [-180,180].                                             
%   deg360        - Converts degrees to the range [0, 360].                                               
%   degunwrap     - Unwraps arrays given in degrees.                                                      
%   discretecolorbar - Plots a colorbar with discrete variation.                                          
%   divgeom       - Geometric decomposition of eddy vorticity flux divergence.                            
%   dlines        - Add diagonal lines to a plot.                                                         
%   doublen       - Interpolates a time series to double its length.                                      
%   ecconv        - Converts between eccentricity measures.                                               
%   ellband       - Bandwidth of modulated elliptical signals in two or three dimensions.                 
%   ellcurves     - Returns curves corresponding to specified ellipse properties.                         
%   elldiff       - Differentiation of modulated elliptical signals.                                      
%   ellipseplot   - Plot ellipses.                                                                        
%   ellparams     - Ellipse parameters of a modulated bivariate or trivariate oscillation.                
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
%   fourier       - The one-sided Fourier frequencies for a given length time series.                     
%   frac          - Fraction: FRAC(A,B)=A./B                                                              
%   hermfun       - Orthonormal Hermite functions. [with F. Rekibi]                                       
%   hermpoly      - Hermite polynomials. [with F. Rekibi]                                                 
%   hlines        - Add horizontal lines to a plot.                                                       
%   hodograph     - Generate hodograph plots (simple and fancy).                                          
%   imlog         - Imaginary part of log: IMLOG(X)=UNWRAP(IMAG(LOG(X)))                                  
%   inellipse     - Locates points on the interior of ellipses.                                           
%   inregion      - Tests whether lat/lon points lie within a specified box.                              
%   instmom       - Univariate and multivariate instantaneous moments.                                    
%   inticks       - Sets the 'tickdir' property of the current axis to 'in'.                              
%   iseven        - True for even integer values; false otherwise.                                        
%   isodd         - True for odd integer values; false otherwise.                                         
%   isridgepoint  - Finds wavelet ridge points using one of several criterion.                            
%   jdeg2rad      - Converts degrees to radians.                                                          
%   jhelp         - Opens linked JLAB help files in Matlab's internal web browser.                        
%   jlab_addpath  - Adds JLAB and JDATA subdirectories to your Matlab search path.                        
%   jlab_allhelp  - Displays the help comments for all JLAB modules.                                      
%   jlab_changes  - Changes to JLAB in each release.                                                      
%   jlab_highlights - Introduction to some of the most useful routines.                                   
%   jlab_index    - Alphabetical index into JLAB and JDATA contents.                                      
%   jlab_install  - Instructions for installing JLAB.                                                     
%   jlab_license  - License statement and permissions for JLAB package.                                   
%   jlab_makefigs - Makes figures for papers by J.M. Lilly.                                               
%   jlab_runtests - Runs a test suite for JLAB package.                                                   
%   jlab_thanks   - Sources of support for developing this software.                                      
%   jmat2         - 2x2 rotation matrix through specified angle.                                          
%   jmat3         - 3x3 rotation matrix through specified angle.                                          
%   jpcolor       - Modified version of PCOLOR appropriate for cell-centered grids.                       
%   jrad2deg      - Converts radians to degrees.                                                          
%   kl2ab         - Converts ellipse parameters Kappa and Lambda to A and B.                              
%   lansey        - The Lansey modification of Cynthia Brewer's "Spectral" colormap.                      
%   lansey_copyright - Copyright statement for the Lansey colormap.                                       
%   latlon2uv     - Converts latitude and longitude to horizontal velocity.                               
%   latlon2xy     - Converts latitude and longitude into local Cartesian coordinates.                     
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
%   makefigs_analytic - has been moved to JLAB_MAKEFIGS.                                                  
%   makefigs_asilomar - has been moved to JLAB_MAKEFIGS.                                                  
%   makefigs_bandwidth - has been moved to JLAB_MAKEFIGS.                                                 
%   makefigs_closedcurves - Makes a sample figure for CLOSEDCURVES.                                       
%   makefigs_curvemoments - Makes a sample figure for CURVEMOMENTS.                                       
%   makefigs_dawson - Makes a sample figure for DAWSON.                                                   
%   makefigs_dlines - Makes a sample figure for DLINES.                                                   
%   makefigs_doublen - Makes a sample figure for DOUBLEN.                                                 
%   makefigs_ecconv - Makes a sample figure for ECCONV.                                                   
%   makefigs_ellband - Makes a sample figure for ELLBAND.                                                 
%   makefigs_ellipseplot - Makes a sample figure for ELLIPSEPLOT.                                         
%   makefigs_ellvel - Makes a sample figure for ELLVEL.                                                   
%   makefigs_hermfun - Makes a sample figure for HERMFUN.                                                 
%   makefigs_inellipse - Makes a sample figure for INELLIPSE.                                             
%   makefigs_instmom - Makes some sample figures for INSTMOM.                                             
%   makefigs_jpcolor - Makes a sample figure for JPCOLOR.                                                 
%   makefigs_jtopo - Makes a sample figure for ABOUT_JTOPO.                                               
%   makefigs_lansey - Makes a sample figure for LANSEY.                                                   
%   makefigs_maternoise - Makes a sample figure for MATERNOISE.                                           
%   makefigs_maternspec - Makes some sample figures for MATERNSPEC.                                       
%   makefigs_morlfreq - Makes a sample figure for MORLFREQ.                                               
%   makefigs_morlwave - Makes a sample figure for MORLWAVE.                                               
%   makefigs_morsearea - Makes a sample figure for MORSEAREA.                                             
%   makefigs_morsebox - Makes a sample figure for MORSEBOX.                                               
%   makefigs_morsefreq - Makes a sample figure for MORSEFREQ.                                             
%   makefigs_morsespace - Makes some sample figures for MORSESPACE.                                       
%   makefigs_morsewave - Makes some sample figures for MORSEWAVE.                                         
%   makefigs_morsies - has been moved to JLAB_MAKEFIGS.                                                   
%   makefigs_mspec - Makes some sample figures for MSPEC.                                                 
%   makefigs_multivariate - has been moved to JLAB_MAKEFIGS.                                              
%   makefigs_patchcontourf - Makes a sample figure for PATCHCONTOURF.                                     
%   makefigs_periodize - Makes a sample figure for PERIODIZE.                                             
%   makefigs_polysmooth - Makes some sample figures for POLYSMOOTH.                                       
%   makefigs_polysmooth2 - Makes a sample figure for POLYSMOOTH using TPJAOS.MAT.                         
%   makefigs_psi2fields - Makes a sample figure for PSI2FIELDS.                                           
%   makefigs_regionplot - Makes some sample figures for REGIONPLOT.                                       
%   makefigs_ridges - has been moved to JLAB_MAKEFIGS.                                                    
%   makefigs_ridgewalk - Makes a sample figure for RIDGEWALK.                                             
%   makefigs_simplepdf - Makes a sample figure for SIMPLEPDF.                                             
%   makefigs_slidetrans - Makes a sample figure for SLIDETRANS.                                           
%   makefigs_stickvect - Makes a sample figure for STICKVECT.                                             
%   makefigs_superfamily - has been moved to JLAB_MAKEFIGS.                                               
%   makefigs_topoplot - Makes a sample figure for TOPOPLOT.                                               
%   makefigs_trackextract - Makes a sample figure for TRACKEXTRACT.                                       
%   makefigs_trajextract - Makes a sample figure for TRAJEXTRACT.                                         
%   makefigs_trajwrap - Makes a sample figure for TRAJWRAP.                                               
%   makefigs_trivariate - has been moved to JLAB_MAKEFIGS.                                                
%   makefigs_twodhist - Makes a sample figure for TWODHIST.                                               
%   makefigs_twodmed - Makes a sample figure for TWODMED.                                                 
%   makefigs_twodstats - Makes a sample figure for TWODSTATS.                                             
%   makefigs_vortex - has been moved to JLAB_MAKEFIGS.                                                    
%   makefigs_wavespecplot - Makes a sample figure for WAVESPECPLOT.                                       
%   makefigs_wavetrans - Makes a sample figure for WAVETRANS.                                             
%   makefigs_widgist - Makes some sample figures for WIGDIST.                                             
%   mat2col       - Compress NAN-padded matrix data into long columns.                                    
%   materncfun    - Returns the normalization function C_ALPHA for a Matern process.                      
%   maternchol    - Cholesky decomposition of Matern and fBm covariances.                                 
%   materncov     - Autocovariance of the Matern random process and variations.                           
%   maternedge    - Long-time cutoff edge for the Matern impulse response function.                       
%   maternimp     - Impulse response function for the Matern random process.                              
%   maternoise    - Realizations of the Matern process and variations, including fBm.  [with A. Sykulski] 
%   maternspec    - Fourier spectrum of the Matern random process and variations.                         
%   matinv        - Fast inversion of arrays of small matrices.                                           
%   matmult       - Matrix multiplication for arrays of matrices.                                         
%   matsave       - Create and save structure of variables as a mat-file.                                 
%   maxmax        - MAXMAX(X)=MAX(X(ISFINITE(X)))                                                         
%   minmin        - MINMIN(X)=MIN(X(ISFINITE(X)))                                                         
%   mom2cum       - Convert moments to cumulants.                                                         
%   morlfreq      - Compute Morlet wavelet carrier frequency given peak frequency.                        
%   morlwave      - Morlet wavelet.                                                                       
%   morseafun     - Returns the generalized Morse wavelet amplitude "a".                                  
%   morsearea     - Time-frequency concentration area of Morse wavelets. [with F. Rekibi]                 
%   morsebox      - Heisenberg time-frequency box for generalized Morse wavelets.                         
%   morsederiv    - Frequency-domain derivatives of generalized Morse wavelets.                           
%   morsefreq     - Frequency measures for generalized Morse wavelets. [with F. Rekibi]                   
%   morsehigh     - High-frequency cutoff of the generalized Morse wavelets.                              
%   morsemom      - Frequency-domain moments of generalized Morse wavelets.                               
%   morseproj     - Projection coefficient for two generalized Morse wavelets.                            
%   morseprops    - Properties of the demodulated generalized Morse wavelets.                             
%   morsespace    - Logarithmically-spaced frequencies for generalized Morse wavelets.                    
%   morsewave     - Generalized Morse wavelets of Olhede and Walden (2002).                               
%   morsexpand    - Generalized Morse wavelets via time-domain Taylor series.                             
%   mspec         - Multitaper power and cross spectra.                                                   
%   msvd          - Singular value decomposition for polarization analysis.                               
%   ncinterp      - Interpolate field from NetCDF file onto specified positions.                          
%   nocontours    - Removes contours from a CONTOURF plot.                                                
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
%   polysmooth    - Smoothing scattered 2D data with local polynomial fitting.                            
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
%   ridgewalk     - Extract wavelet transform ridges, including bias estimates.                           
%   rot           - Complex-valued rotation:  ROT(X)=EXP(SQRT(-1)*X)                                      
%   sampletimes   - Computes mean sampling intervals and their statistics.                                
%   sig2latlon    - Converts an oscillatory signal to lat/lon displacements.                              
%   simplepdf     - Gaussian, uniform, Cauchy, and exponential pdfs.                                      
%   sleptap       - Calculate Slepian tapers.                                                             
%   slidetrans    - Sliding-window ('moving-window') Fourier transform.                                   
%   specdiag      - Diagonalize a 2 x 2 spectral matrix.                                                  
%   sphere2uvw    - Converts a 3D spherical vector to a 3D Cartesian vector.                              
%   spherecurl    - Curl of a vector field on the surface of a sphere.                                    
%   spheredist    - Computes great circle distances on a sphere.                                          
%   spherediv     - Divergence of a vector field on the surface of a sphere.                              
%   spheregrad    - Gradient of a field on the surface of a sphere.                                       
%   spherelap     - Laplacian of a field on the surface of a sphere.                                      
%   spheresort    - Sorted great circle distances to nearby points on the earth.                          
%   squared       - Squares the modulus of its argument:  SQUARED(X)=ABS(X).^2                            
%   standalone    - Create stand-alone version of an m-file, including dependencies.                      
%   stickvect     - Plots "stick vectors" for multicomponent velocity time series.                        
%   tidefreq      - Frequencies of the eight major tidal components.                                      
%   timeseries_boundary - Apply boundary conditions to data before transform.                             
%   tmat          - 2x2 complex grouping matrix.  TMAT = [1  i; 1 -i] / SQRT(2)                           
%   to_grab_from_caller - Returns a string to grab variable values from caller.                           
%   to_overwrite  - Returns a string to overwrite original arguments.                                     
%   topo_copyright - Copyright statement for the Smith and Sandwell topography.                           
%   topoplot      - Plot regional or global topography at one-sixth degree resolution.                    
%   trackextract  - Extracts alongtrack altimetry segments within given region.                           
%   trajchunk     - Converts Lagrangian trajectories into chunks based on the Coriolis period.            
%   trajextract   - Extracts Lagrangian trajectory segments within given region.                          
%   trajunwrap    - Unwraps Lagrangian trajectories from a periodic domain.                               
%   trajwrap      - Wraps Lagrangian trajectories to fit within a periodic domain.                        
%   twodhist      - Two-dimensional histogram.                                                            
%   twodmed       - Median value of a function of two variables.                                          
%   twodsort      - Distances from data points to nearby grid points.                                     
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
%   xlin          - Sets x-axis scale to linear.                                                          
%   xlog          - Sets x-axis scale to logarithmic.                                                     
%   xoffset       - Offsets lines in the x-direction after plotting.                                      
%   xtick         - Sets locations of x-axis tick marks.                                                  
%   xy2latlon     - Converts local Cartesian coordinates into latitude and longitude.                     
%   xyz2latlon    - Converts 3D Cartesian coordinates into latitude and longitude.                        
%   yearfrac      - Converts a DATENUM into 'year.fraction' and 'month.fraction'.                         
%   ylin          - Sets y-axis scale to linear.                                                          
%   ylog          - Sets y-axis scale to log.                                                             
%   yoffset       - Offsets lines in the y-direction after plotting.                                      
%   ytick         - Sets locations of y-axis tick marks.                                                  
%   ztick         - Sets locations of z-axis tick marks.                                                  
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015--2016 J.M. Lilly --- type 'help jlab_license' for details

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
jdata_dir=[jlab_dir(1:end-4) 'jdata'];

for k=1:2
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

disp('%JLAB_INDEX  Alphabetical index into JLAB and JDATA contents.')
disp('%')
disp('% JDATA index')
disp(comments2)
disp('%')
disp('% JLAB index')
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
