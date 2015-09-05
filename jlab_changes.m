%JLAB_CHANGES   Changes to JLAB in each release.
%
%   Changes new in version 1.6.1
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
%   (C) 2002--2015 J.M. Lilly --- type 'help jlab_license' for details

help jlab_changes

%   These were help back from version 0.95.
%
%   gaussprofile - Wavelet transform profile of a Gaussian with a Gaussian wavelet.
%   eddycensus   - Census of coherent eddies from streamfunction snapshots.
%   curveflow  Transport into and out of a region bounded by a closed curve.
%   eddyslice    - Slice through an eddy due to a turning background flow.
%   eddycurve    - Trajectory of an eddy center in a turning background flow.
%   eddyguess    - Eddy properties in a current meter record from a wavelet transform.
%   eddyfit      - Best fit parameters for an eddy advected past a current meter mooring.
%   wavetransderiv  - Rapidly calculate the wavelet transform of a signal's derivative.