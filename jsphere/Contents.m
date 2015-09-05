% jSphere:  Tools for latitude, longitude, and spherical geometry
%
% Simple lat / lon conversions
%   deg180     - Converts degrees to the range [-180,180].                        
%   deg360     - Converts degrees to the range [0, 360].                          
%   jdeg2rad   - Converts degrees to radians.                                     
%   jrad2deg   - Converts radians to degrees.   
%
% Distances and regions 
%   inregion   - Tests whether lat/lon points lie within a specified box.         
%   spheredist - Computes great circle distances on a sphere.            
%
% Earth constants and frequencies
%   radearth   - The radius of the earth in kilometers.   
%   corfreq    - Coriolis frequency in radians per hour.   
%   tidefreq   - Frequencies of the eight major tidal components.
%
% Coordinate transformations 
%   latlon2xy  - Converts latitude and longitude into local Cartesian coordinates.
%   latlon2xyz - Converts latitude and longitude into 3D Cartesian coordinates.   
%   xy2latlon  - Converts local Cartesian coordinates into latitude and longitude.
%   xyz2latlon - Converts 3D Cartesian coordinates into latitude and longitude.   
%   sphere2uvw - Converts a 3D spherical vector to a 3D Cartesian vector. 
%
% Differential operators
%   spherecurl - Curl of a vector field on the surface of a sphere.               
%   spherediv  - Divergence of a vector field on the surface of a sphere.         
%   spheregrad - Gradient of a field on the surface of a sphere.                  
%   spherelap  - Laplacian of a field on the surface of a sphere.   
%
% Plotting tools
%   latratio   - Set plot aspect ratio for latitude / longitude plot.    
%   lonshift   - Shifts longitude origin for plotting purposes.                   
%   regionplot - Plots a simple box indicating a latitude / longitude region.    

help jSphere