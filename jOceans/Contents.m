% jOceans:  Oceanography-specific data and model analysis tools
%
% Conversions for Lagrangian trajectories
%   latlon2uv   - Converts latitude and longitude to horizontal velocity.  
%   uv2latlon   - Integrates horizontal velocity to give latitude and longitude.  
%
% Manipulating Lagrangian trajectories [see also jCell]
%   trajextract  - Extracts Lagrangian trajectory segments within given region.        
%   trajunwrap   - Unwraps Lagrangian trajectories from a periodic domain.              
%   trajwrap     - Wraps Lagrangian trajectories to fit within a periodic domain.       
%   trajchunk    - Chunks Lagrangian trajectories based on the Coriolis period.
%   griddrifters - Average drifter velocities onto a space/time 3D grid.
%
% Idealized numerical model tools
%   psi2fields   - Velocity and other fields from the streamfunction. [with P.E. Isachsen]    
%   periodize    - Returns a doubly periodic version of an input array.   
%
% Lagrangian eddy identification and analysis
%   eddyridges    - Coherent eddy ridges from Lagrangian trajectories.   
%   noisedrifters - Create a noise Lagrangian dataset matching mean and variance.
%   eddylevels    - Eddy ridge significance levels using the survival function.
%
% Eulerian eddy identification and analysis
%   inellipse    - Locates points on the interior of ellipses.                          
%   closedcurves - Locate and interpolate closed curves in a possibly periodic domain.
%   curvemoments - Centroid, area, and many other moments of a closed curve.          
%   divgeom      - Geometric decomposition of eddy vorticity flux divergence.     
%   eddyfit2d    - Least squares fit of 2D velocity data to an eddy profile.
%   simpleddy     - Streamfunction, velocity, and vorticity for various eddy profiles.
%
% Plotting tools for mooring data
%   hodograph  - Generate hodograph plots (simple and fancy).                                   
%   provec     - Generate progressive vector diagrams (simple and fancy).             
%   stickvect  - Plots "stick vectors" for multicomponent velocity time series. 
%
% Alongtrack altimetry tools
%   trackextract  - Extracts alongtrack altimetry segments within given region.
%
% Wind-driven surface currents
%    windtrans  - Ekman-like transfer-functions for the wind-driven response.
%
% NetCDF tools
%   ncinterp    - One-line interpolation from 3D lat/lon/time field in NetCDF file.
%
% Date and time
%   monthstats  - Mean month and standard deviation using circular statistics.
%
% Topography tools and data
%   jtopo.mat   - One-sixth degree global topography, from Smith and Sandwell + IBCAO.              
%   topoplot    - Plot regional or global topography at one-sixth degree resolution.
% 
% Low-level functions
%   besselktilde  - K-type Bessel function after factoring off exponential decay.
%   besselitilde  - I-type Bessel function after factoring off exponential growth.
%   curveinterp   - Interpolate a field or its gradient onto a set of curves.           
%   orbitbreaks   - Separate orbit into passes based on turning points. 
%   trianglepath  - Moving instrument path composed of adjacent triangles.
%   topo_copyright - Copyright statement for the Smith and Sandwell topography.
%
%  See also jData.


help joceans
