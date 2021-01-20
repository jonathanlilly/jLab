%JLAB_CHANGES   Changes to JLAB in each release.
%
%   Changes new in version 1.6.9
%
%   Figure making for two new publications, see JLAB_MAKEFIGS:
%      
%   Lilly and Perez-Brunius (2021b). Extracting statistically significant 
%       eddy signals from large Lagrangian datasets, with application to 
%       the Gulf of Mexico. In review at Nonlinear Processes in Geophysics.
%
%   Lilly and Perez-Brunius (2021a). A gridded surface current product for
%       the Gulf of Mexico from consolidated  drifter measurements.  
%       Accepted to Earth Science System Data.
%
%   Major new functions for Lagrangian eddy analysis and statistics:
%
%   eddyridges    - Coherent eddy ridges from Lagrangian trajectories.   
%   noisedrifters - Create a noise Lagrangian dataset matching mean and variance.
%   eddylevels    - Eddy ridge significance levels using the survival function.
%   griddrifters  - Average drifter velocities onto a space/time 3D grid.
%
%   Other new functions:
%
%   spheretrans   - Wavelet transform for oscillations on the surface of a sphere.
%   ridgemult     - Ridge multiplicity, the number of simultaneous ridges present.
%   catstruct     - Concatenates the array elements of a cell array of structures.
%   structindex   - Applies an index to all array-valued fields in a structure.
%   mconf         - Confidence intervals for the multitaper spectral estimate.
%   chisquared    - The chi-squared distribution.
%   eddyfit2d     - Least squares fit of 2D velocity data to an eddy profile.
%   simpleddy     - Streamfunction, velocity, and vorticity for various eddy profiles.
%   monthstats    - Mean month and standard deviation using circular statistics.
%   jhanning      - Hanning window.
%   celldot       - Dot product for arrays of column vectors.
%   haxby         - The Haxby colormap.  
%   seminfhaxby   - The seminf-Haxby colormap.
%
%   Updated datasets:  
%   
%   JTOPO global topography is now at 1/12th resolution.
%   
%   Changes and improvements:
%
%   SLEPTAP vast speed improvements for handling multiple time series.
%   LINECOLOR bugfixes.
%   PRINTALL organization improvement and bugfix for R2020a.
%   MATERNFIT now includes an option to average over ensemble members.
%   MATERNFIT, MATERNSPEC, MATERNCOV, and MATERNOISE refactoring, support 
%        for new extensions of the Matern process, support real-valued as
%        well as complex-valued processes, and input argument changes. 
%   MSPEC now uses two-sided, not one-sided, normalization for real data.
%   POLYSMOOTH considerable speed improvements. 
%   SPHERESORT and TWODSORT improvements and major change to output format.
%   SPHERELAP now computes all the terms in the Hessian matrix.
%   AB2KL bugfix for B=0.
%   READTOPO bugfixes when using interpolation option.
%   PERIODINDEX new, faster algorithm based on ridge age.
%   VMOMENT bugfix affecting all odd moments.
%   LETTERLABELS bugfix to avoid legend handles, if present. 
%   CELLPACK functionality change to key off of the first input argument.
%   VSTD redefined for complex-valued arguments.
%   YEARFRAC bugfix for NUM consisting of a ND array.
%   ELLROSSBY new definition of Rossby number based on inferred vorticity.
%   TRAJEXTRACT and POLYSMOOTH bugfixes for sample figures. 
%   TWODSTATS and TWODMED fixed sample figures.
%   LINECOLOR can now take either a colormap matrix or a name as input.
%   NCLOAD can now handle filenames including hyphens '-'.
%   NCLOAD now supports cell conversion for the 'trajectory' feature type.
%   NCLOAD now converts all non-double numeric variables into doubles.
%   NCLOAD is now self-contained.
%   -----------------------------------------------------------------------
%
%   Changes new in version 1.6.6
%
%   Figure making for a new publications:
%      
%   Lilly (2018).  Kinematics of a fluid ellipse in a linear flow. Fluids, 
%        3 (1) 16: 1--29.
%
%   New functions:
% 
%   jfig         - Shorthand for tweaking figures.
%   jprint       - Print to a specified directory and crop the resulting file.
%   ncload       - Load all variables from a NetCDF file and unpack any columns.
%   ridgetrim    - Trim edge effect regions from wavelet ridges.
%   colorquant   - Set the color axis according to quantiles of the data.
%   ellpol       - Polarization parameters of an elliptical signal.
%   printall     - Print and close all open figures.
%   cellmed      - Median value of each element a cell array.
%   cellsum      - Sum of each element a cell array, possibly weighted.
%   interplatlon - Interpolation for working with latitude and longitude.
%   jhanning      - Hanning window.
%
%   Updated datasets:  All datasets now include both a mat-file version and
%        a NetCDF version.  The NetCDF versions are larger but load faster.
%        Use the jLab routine NCLOAD to load a whole NetCDF file at once.
%
%   floats.mat and floats.nc     - Major update to the historical dataset of 
%        eddy-resolving subsurface floats.  See about_floats.
%   drifters.mat and drifters.nc - The global surface drifter dataset, 
%       updated through April 1, 2018.  See about_drifters.
%   tpjaos.mat and tpjaos.nc     - Alongtrack SSH anomalies from the 
%       Beckley dataset, updated through May 30, 2018.  See about_tpjaos.
%   sandwell.mat and sandwell.nc - One minute resolution topography data 
%       from Smith and Sandwell. See about_sandwell.
%   ibcao.mat and ibcao.nc       -  International Bathymetric Chart of the 
%      Arctic Ocean topography. See about_ibcao.
%
%   Major bugfix
%
%   SPHERELAP for computing the Laplacian in spherical geometry 
%   unfortunately had a major bug in earlier versions that rendered its
%   output incorrect.  This has been corrected in version 1.6.6.
%
%   Changes and improvements:
% 
%   VSIZE changed to return 1's rather than 0's for singular dimensions.
%   NCINTERP bugfix to correctly handle [-180,180] or [0,360] longitudes.
%   NCINTERP buxfix affecting processing chunks in non-periodic domains.
%   ELLIPSEPLOT bugfix for interacting with M_MAP.
%   JLAB_ADDPATH bugfix to add root JDATA directory.
%   MSPEC now computes the spectrum along an arbitrary dimension.
%   MSPEC now removes the mean by default, but no longer detrends.
%   ANATRANS improved to correctly handle the Nyquist frequency.
%   DIVGEOM now handles different DX and DY; note input argument change.
%   FINDFILES bugfix for 'include' and 'exclude' flags.  
%   TOPOPLOT no longer plots continental shelves by default.
%   CELLGRID now accepts array-valued DT arguments.
%   CELLSPLIT now accepts array-valued TOL arguments.
%   TWODSTATS bugfix for complex-valued input using HISTCOUNTS2 algorithm.
%   MATINV now works for matrices up to dimension 12. 
%   LATLON2XY added new correction for small angles and additional tests.
%   SPHERESORT and TWODSORT now have options for truncating size of output.
%   POLYSMOOTH has major changes, including speed and memory improvements.
%   MATERNFIT now supports using NLopt, https://nlopt.readthedocs.io. 
%   MATERNFIT now works using Nelder-Mead without the Optimization Toolbox. 
%   MATERNFIT corrections to AICC output field.
%   RIDGEWALK now supports trimming edge effects using RIDGETRIM.
%   RIDGEWALK simplifications and speed improvements.
%   RIDGEWALK bugfix to sample figure generation.
%   PERIODINDEX simplifications and bugfix.
%   FOURIER now outputs one- or two-sided frequency arrays.
%   FOURIER outputs a cell array of frequency arrays given vector input.
%   CELLFILL, -STRIP, and -PACK now use INFs (not NaNs) as the data flag.
%   CELLMEAN and -STD now support weighted means and standard deviations.
%   ELLPARAMS input format simplifications.
%   READTOPO bugfixes and additional testing.
%   -----------------------------------------------------------------------
%
%   Changes new in version 1.6.5
%         
%   Support for a new publication, see module jMatern:
%
%   Lilly, Sykulski, Early, and Olhede (2017).  Fractional Brownian motion,
%      the Matern process, and stochastic modeling of turbulent dispersion.
%      Nonlinear Processes in Geophysics, 24: 481--514.
%
%   New functions:
%
%   maternfit  - Parametric spectral fit to the Matern form. [with A. Sykulski]
%   
%   Changes and improvements:
%
%   XTICK, YTICK, ZTICK bugfix for specifying empty tickmark sets.  
%   FINDFILES now has an option to ignore all directories named 'private'.
%   STANDALONE bugfix for recursive searches and support for files lists.
%   STANDALONE now uses 'private' sub-directory for dependency files.
%
%   -----------------------------------------------------------------------
%
%   Changes new in version 1.6.4
%
%   PROVEC now expecting sampling interval as first argument, in hours.
%   STICKVECT bugfix for cutting off the end of the time series.
%   STICKVECT improved commenting for how to choose the scale factor SCALE.
%
%   -----------------------------------------------------------------------
%
%   Changes new in version 1.6.3
%
%   This is a major new release, which includes the following changes:
%
%   Support for a new publication, see module jWavelet:
%
%   Lilly (2017).  Element analysis: a wavelet-based method for analyzing
%      time-localized events in noisy time series. Submitted.  Available at
%      http://jmlilly.net/jmlpubs.html.
%
%   New functions for element analysis using generalized Morse wavelets:
%
%   transmax     - Locates modulus maximum points of wavelet transform.
%   transmaxdist - Distributions of wavelet transform maxima in noise.
%   morseregion  - Generalized Morse wavelet time-frequency concentration region.
%   isomax       - Returns those transform maxima that are isolated from others.
%   maxprops     - Returns properties of wavelet transform maxima.
%   max2eddy     - Converts transform maxima into oceanic coherent eddy properties.
%
%   New function for fast, accurate remapping of non-plaid grids on the sphere:
%
%   sphereinterp  - Fast linear interpolation on the sphere from non-plaid grids.
%
%   Updated datasets
%
%   drifters.mat  - The global surface drifter dataset from NOAA's Global
%        Drifter Program, updated through Sept. 2016.  See about_drifters.
%   tpjaos.mat    - Alongtrack sea surface height anomalies from the 
%        Beckley dataset, updated through Oct. 12, 2016.  See about_tpjaos.
% 
%   Changes and improvements
%
%   JPCOLOR now SQUEEZEs input arguments by default. 
%   POLYSMOOTH improved organization of output arguments, speed improvments,
%        and improved documentation.
%   SPHERESORT fixed a bug in the parallelized version that led to a shift
%        in output arguments if the output included empty latitude bands.  
%   LATLON2UV improvemented argument handling for matrix input.
%   NCINTERP now supports input longitudes as [-180,180] or [0, 360].
%   YEARFRAC now passes Inf values through, as well as NaNs.
%   ELLIPSEPLOT and LINESTYLE now support the new Matlab colors.
%   RIDGEWALK now supports cell array input.
%   READTOPO bugfix for longitude regions spanning 360 degrees.
%   MORSEHIGH simplifications.
%   MORSESPACE minor re-definition of low-frequency cutoff parameter.
%   VSUM modified for backwards compatibility with earlier Matlab versions.
%   DAWSON has been renamed JDAWSON, to avoid confusion with Matlab's
%        version in the Symbolic Math Toolbox. 
%   SLEPTAP added sign check for non-interpolated data tapers.
%   -----------------------------------------------------------------------
%
%   Changes new in version 1.6.2 
%
%   This is a major new release, which includes the following changes:
%         
%   Support for a new publication, see module jMatern:
%
%   Lilly, Sykulski, Early, and Olhede (2016).  Fractional Brownian motion,
%      the Matern process, and stochastic modeling of turbulent dispersion.
%      Available at http://jmlilly.net/jmlpubs.html.   
%
%   New functions
%
%   ncinterp      - Interpolate field from NetCDF file onto specified positions.
%   cellget       - Indexes a cell array of numerical arrays by ID number.
%   cellimit      - Limits the ranges of times in a cell array of numerical arrays.
%   cellstd       - Standard deviation of each element a cell array.
%   cellpack      - Removes all NaN values from cell arrays of numeric arrays.
%   degunwrap     - Unwraps arrays given in degrees.
%   trackextract  - Extracts alongtrack altimetry segments within given region.
%
%   Changes and improvements
%
%   DRIFTERS.MAT is now updated with data through December 2015.
%   ABOUT_DRIFTERS contains updated and simplified internal processing.
%   TPJAOS.MAT updated through November 2015.
%
%   TWOSPECPLOT bugfix for omitting tidal lines.
%   TWOSPECPLOT bugfix swapping Cartesian spectra.
%   POLYSMOOTH corrected to account for missing data values in sorted input.
%   TRACKEXTRACT now strips out any empty tracks in the region.
%   JPCOLOR now works for non-uniformly spaced axis arrays.
%   NCINTERP is now greatly accelerated by loading large chunks into memory.
%   PERIODINDEX updated to handle new format output but RIDGEWALK.
%   CELLPLOT bugfix to handle empty cells.
%   YEARFRAC modified to handle potential empty cells.
%   RIDGEWALK output formats simplified.
%   CELL2COL and COL2CELL now work with arrays having multiple columns.
%   LETTERLABELS updated for recent changes to handle graphics conventions.
%   TOPOPLOT now supports linestyles in LINESTYLE format.
%   TOPOPLOT bugfix for x-axis limits exceeding 360 degrees.
%   TOPOPLOT bugfix when contours are input. 
%   CELLFILL bugfix to prevent transposition of cell arrays.
%   CELLGRID, CELLMEAN, CELLSTD, and SAMPLETIMES now support parallelization.
%   CELLGRID bugfix for data containing some or all NANs.  
%   CELLGRID now uses the 'pchip' method by default, as does CELLFILL.
%   CELLPACK bugfix to address inaccurate report of points filled.
%   LATLON2UV now optionally returns the acceleration.
%   LATLON2UV bugfix for empty or length one datasets.
%   LATLON2UV now longer has the option to return U and V separately. 
%   LANSEY bugfix for returning specified size of colormap.
%   TRAJFILL has been removed.  Use CELLFILL instead. 
%   TRAJEXTRACT can now handle input numerical arrays such as ID number.
%   TRAJEXTRACT input argument order change.
%   TRAJEXTRACT now supports variable overwriting.
%   TRAJCHUNK bugfix for length one or zero cell elements. 
%   TRAJCHUNK can now handle input numerical arrays such as ID number. 
%   TWOSPECPLOT changes to make it easier to overlay spectral plots.
%   MATERNCHOL bugfix for numerical noise in imaginary part of matrix.
%   VMEAN, VSUM, VMOMENT, and VSTD speed improvements.
%   VMEAN, VSUM, VMOMENT, and VSTD now no longer exclude Infs, only NaNs.
%
%   -----------------------------------------------------------------------
%
%   Changes new in version 1.6.1
%
%   jLab is now on GitHub at https://github.com/jonathanlilly/jLab.
%
%   Changes and improvements
%
%   WAVETRANS now supports parallelization.
%   LINECOLOR bugfix for colormaps of lengths different from 64.
%   CELLFIRST now also supports a 'last' flag to return the last element.
%   YEARFRAC now also returns 'month.fraction' in additon to 'year.fraction.'
%   CELLGRID and CELLFILL bugfixes for time series of length one or zero.
%   SAMPLETIMES bugfix for time series of length one.
%   SPHEREDIST bugfix for row vectors input.
%   TOPOPLOT bugfix for vector of depth contours input.
%
%   -----------------------------------------------------------------------
%
%   Changes new in version 1.6
%
%   New datasets
%
%   about_ibcao - International Bathymetric Chart of the Arctic Ocean topography.
%
%   New functions
%
%   jlab_index - Alphabetical index into JLAB and JDATA contents.
%   standalone - Create stand-alone version of an m-file, including dependencies.
%   yearfrac   - Convert date from 'datenum' format to 'year.fraction'.
%   lansey     - The Lansey modification of Cynthia Brewer's "Spectral" colormap.
%  
%   celldiv    - Division acting on each element in a cell array.
%   cellstrip  - Strips NaN values from the beginnings or ends of cell arrays.
%   cellgrid   - Interpolate a cell array of numeric arrays onto a regular grid.
%   cellsplit  - Splits cell arrays of numeric arrays at data gaps.
%   cellprune  - Removes all empty cells, or cells less than a specified length.
%   cellfill   - Fills missing data marked by NaNs in a cell array.
%   cellfirst  - Returns the first element of each entry in a cell array.
%
%   Changes and improvements
%
%   All example figures now echo code to the Matlab command window.
%
%   MATERNCOV, MATERNSPEC, and MATERNOISE extensions and simplifications.
%   ABOUT_JTOPO dataset now includes high-latitude topography data from IBCAO.  
%   MSPEC, SLEPTAP, and MATERNSPEC now support parallelization.  
%   TWODHIST and TWODSTATS are now vastly faster thanks to Matlab's 
%       HISTCOUNTS2, available as of release 2015b.
%   TWODHIST and TWODSTATS now support parallelization.
%   LATLON2UV now supports parallelization.
%   CELLCHUNK now supports parallelization.
%   CELLCHUNK bugfix for 50% overlap mode. 
%   ARESAME now works with cell arrays of numerical arrays.
%   CORFREQ is now defined to be positive in the Northern Hemisphere, 
%       and negative in the Southern Hemisphere, following convention.
%   MSPEC bugfix for periodogram applied to multi-component datasets.
%   TOPOPLOT now plots topography underneath other items on the plot.
%   MATERNCHOL improved handling of numerical non-positive-definiteness.
%   VSWAP now works with cell arrays.  
%   ELLIPSEPLOT bugfix for empty index input.
%   PATCHCONTOURF bugfix for use with M_MAP.
%   _________________________________________________________________
%
%   Changes new in version 1.5 
%
%   This is a major new release, which includes the following changes:
%       
%     --- Compatible with the graphics changes introduced in Matlab R2014b 
%     --- Completely new organization with a more modular approach
%     --- Many obsolete functions removed and redundancies eliminated
%
%   New functions
%
%   divgeom        - Geometric decomposition of eddy vorticity flux divergence.
%   blurspec       - Returns the blurred and aliased spectrum given the autocovariance.
%   findpath       - Returns the full pathname of a directory on the Matlab search path.
%   patchcontourf  - Generate filled contours using patches, with specified colors.
%   cellchunk      - Converts cell array data into uniform length 'chunks'.
%   trajchunk      -  Converts Lagrangian trajectories into chunks based on the Coriolis period.
%   cellmean       - Mean value of each element a cell array or set of cell arrays.
%   uv2latlon      - Integrates horizontal velocity to give latitude and longitude.
%   sig2latlon     - Converts an oscillatory signal to lat/lon displacements.
%   fminsearchbnd  - fminsearch, but with bound constraints by transformation. [By J. D'Errico]
%   whichdir       - Returns directory name containing file in search path.
%   jlab_allhelp   - Displays the help comments for all JLAB modules.
%
%   Improved datasets
%
%   about_drifters - Modified to include data through June 2014.  Please
%      note that the earlier verison of drifters.mat was missing half the 
%      data on account of a read error.  The new version corrects this.
%
%   Minor changes and improvements
%
%   SPHERESORT speed improvements, bug fix, and support for parallel computing.
%   READTOPO now works with Smith and Sandwell v. 18.1
%   FINDFILES now optionally searches directories recurvisely. 
%   PACKBOTH, PACKROWS, and PACKCOLS are now combined into PACKFIG.
%   JLAB_MAKEFIGS is now one file that makes figures for recent published papers. 
%   MATERNSPEC bugfix for input parameter arrays of length greater than one.
%   FLOATREGION now rejects individual trajectory segements outside of region.
%   CELLPLOT fixed minor bug in which CELLPLOT did not plot the first point.
%   TWOSPECPLOT no longer opens a new figure window by default.
%   RIDGEWALK now outputs ridges in a cell array, with one ridge per cell.  
%
%   ELLIPSEPLOT, MORSEBOX, FONTSIZE, TOPOPLOT, NOCONTOURS, TOPOPLOT, 
%   DISCRETECOLORBAR, JLAB_MAKEFIGS, and sample figures modified to work 
%   with Matlab version R2014b.
%
%   ELLVEL, ELLROSSBY, ELLBAND, ELLPARAMS, ELLRAD, and RIDGELEN now all 
%   correctly return empty arrays, or empty cell arrays, given empty input.
%
%   MATERNSPEC, MATERNCOV, MATERNIMP, and MATERNOISE improvements.  All 
%   extended to accommodate non-unit sample rates.
%
%   Removed a large number of obsolete functions; consolidated many others.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002--2019 J.M. Lilly --- type 'help jlab_license' for details

help jlab_changes

%   These were held back from version 0.95.
%
%   gaussprofile - Wavelet transform profile of a Gaussian with a Gaussian wavelet.
%   eddycensus   - Census of coherent eddies from streamfunction snapshots.
%   curveflow  Transport into and out of a region bounded by a closed curve.
%   eddyslice    - Slice through an eddy due to a turning background flow.
%   eddycurve    - Trajectory of an eddy center in a turning background flow.
%   eddyguess    - Eddy properties in a current meter record from a wavelet transform.
%   eddyfit      - Best fit parameters for an eddy advected past a current meter mooring.
%   wavetransderiv  - Rapidly calculate the wavelet transform of a signal's derivative.