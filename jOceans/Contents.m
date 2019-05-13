% jOceans:  Oceanography-specific data and model analysis tools
%
% Conversions for Lagrangian trajectories
%   latlon2uv   - Converts latitude and longitude to horizontal velocity.  
%   uv2latlon   - Integrates horizontal velocity to give latitude and longitude.  
%
% Manipulating Lagrangian trajectories [see also jCell]
%   trajextract - Extracts Lagrangian trajectory segments within given region.        
%   trajunwrap  - Unwraps Lagrangian trajectories from a periodic domain.              
%   trajwrap    - Wraps Lagrangian trajectories to fit within a periodic domain.       
%   trajchunk   - Converts cell array data into chunks based on the Coriolis period.
%
% Idealized numerical model tools
%   psi2fields   - Velocity and other fields from the streamfunction. [with P.E. Isachsen]    
%   periodize    - Returns a doubly periodic version of an input array.   
%
% Eulerian eddy identification and analysis
%   inellipse    - Locates points on the interior of ellipses.                          
%   closedcurves - Locate and interpolate closed curves in a possibly periodic domain.
%   curvemoments - Centroid, area, and many other moments of a closed curve.          
%   divgeom      - Geometric decomposition of eddy vorticity flux divergence.     
%
% Plotting tools for mooring data
%   hodograph  - Generate hodograph plots (simple and fancy).                                   
%   provec     - Generate progressive vector diagrams (simple and fancy).             
%   stickvect  - Plots "stick vectors" for multicomponent velocity time series. 
%
% Alongtrack altimetry tools
%   trackextract  - Extracts alongtrack altimetry segments within given region.
%
% NetCDF tools
%   ncinterp    - One-line interpolation from 3D lat/lon/time field in NetCDF file.
%
% Topography tools and data
%   jtopo.mat   - One-sixth degree global topography, from Smith and Sandwell + IBCAO.              
%   topoplot    - Plot regional or global topography at one-sixth degree resolution.
% 
%  See also jData.

%  Low-level functions
%  topo_copyright - Copyright statement for the Smith and Sandwell topography.
%  curveinterp  - Interpolate a field or its gradient onto a set of curves.           
%  orbitbreaks  - Separate orbit into passes based on turning points.          

help joceans
