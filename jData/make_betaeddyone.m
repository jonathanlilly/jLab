function[varargout]=make_betaeddyone(varargin)
%MAKE_BETAEDDYONE  Create the BetaEddyOne nonlinear QG eddy simulation. 
%
%   MAKE_BETAEDDYONE creates the BetaEddyOne simulation:
%
%   Early, J. J., and J. M. Lilly (2023). BetaEddyOne: A long-lived 1.5 
%      layer quasigeostrophic eddy on a beta plane, with Lagrangian 
%      particles (1.0.2) [Data set]. Zenodo. 
%      https://doi.org/10.5281/zenodo.8200055
%
%   To run it, you'll need to install the GLOceanKit from GitHub
%
%       https://github.com/Energy-Pathways-Group/GLOceanKit
%
%   specifically the netCDF-conventions branch.  You'll also need to
%   recursively set the Matlab path to include GLOceanKit/Matlab and all 
%   subdirectories.  
%
%   'make_betaeddyone --create' generates BetaEddyOne.nc.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2024 J.M. Lilly and J.J. Early --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--create')
    make_betaeddyone_integrate;
    make_betaeddyone_metadata;
elseif strcmp(varargin{1}, '--integrate')
    make_betaeddyone_integrate;
elseif strcmp(varargin{1}, '--metadata')
    make_betaeddyone_metadata;
end

end

%--------------------------------------------------------------------------
function[]=make_betaeddyone_integrate
%   This script makes the simulation BetaEddyOne.  

Lx = 2000e3;
Ly = 1000e3;

Nx = 512;
Ny = 256;

latitude = 24;

wvt = WVTransformSingleMode([Lx, Ly], [Nx, Ny], h=0.8, latitude=latitude);

x0 = 3*Lx/4;
y0 = Ly/2;
A = 0.15;
L = 80e3;
wvt.setSSH(@(x,y) A*exp( - ((x-x0).^2 + (y-y0).^2)/L^2) );

%figure, pcolor(wvt.x/1000,wvt.y/1000,wvt.ssh.'), shading interp, axis equal
%colorbar('eastoutside')
%xlim([min(wvt.x) max(wvt.x)]/1000); ylim([min(wvt.y) max(wvt.y)]/1000);

[K,L,~] = wvt.kljGrid;
outputVar = WVVariableAnnotation('zeta',{'x','y','z'},'1/s', 'vertical component of relative vorticity');
outputVar.attributes('short_name') = 'ocean_relative_vorticity';
fs = @(wvt) wvt.diffX(wvt.v) - wvt.diffY(wvt.u); % simple definition, but computationally inefficient
f = @(wvt) wvt.transformToSpatialDomainWithF(-(wvt.g/wvt.f) * (K.^2 + L.^2) .* wvt.A0t);
wvt.addOperation(WVOperation('zeta',outputVar,f),overwriteExisting=1);

outputVar = WVVariableAnnotation('nu',{'x','y','z'},'1/s', 'normal strain');
fs = @(wvt) wvt.diffX(wvt.u) - wvt.diffY(wvt.v); % simple definition, but computationally inefficient
f = @(wvt) wvt.transformToSpatialDomainWithF( (wvt.g/wvt.f) * (2*K.*L) .* wvt.A0t);
wvt.addOperation(WVOperation('nu',outputVar,f));

outputVar = WVVariableAnnotation('sigma',{'x','y','z'},'1/s', 'shear strain');
fs = @(wvt) wvt.diffX(wvt.v) + wvt.diffY(wvt.u); % simple definition, but computationally inefficient
f = @(wvt) wvt.transformToSpatialDomainWithF( -(wvt.g/wvt.f) * (K.^2 - L.^2) .* wvt.A0t);
wvt.addOperation(WVOperation('sigma',outputVar,f));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the integrator with the model
model = WVModel(wvt,nonlinearFlux=QGPVE(wvt,shouldUseBeta=1,u_damp=wvt.uMax));

% set initial positions for a bunch of floats
[xFloat,yFloat] = ndgrid(wvt.x,wvt.y);
xFloat = reshape(xFloat,1,[]);
yFloat = reshape(yFloat,1,[]);
nTrajectories = length(xFloat);
model.setDrifterPositions(xFloat,yFloat,[],'ssh','qgpv','u','v','zeta','nu','sigma');

model.setupIntegrator(timeStepConstraint="advective", outputInterval=86400);

model.createNetCDFFileForModelOutput('BetaEddyOne.nc',shouldOverwriteExisting=1);
model.setNetCDFOutputVariables('A0','ssh','qgpv','u','v','zeta','nu','sigma');
model.integrateToTime(365*86400);

ncfile = model.ncfile;
% [x,y] = ncfile.readVariables('drifter_x','drifter_y');
% qgpv = ncfile.readVariables('drifter_qgpv');
% 
% figure, plot(x.',y.')
end
%--------------------------------------------------------------------------


function[]=make_betaeddyone_metadata
%   This script adds NetCDF-compliant metadata to BetaEddyOne.

ncfile = NetCDFFile('BetaEddyOne.nc');

%% Highly recommended ACDD conventions
% Descriptions of these four attributes are found here:
% https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3#Highly_Recommended
ncfile.addAttribute('title','BetaEddyOne'); 
ncfile.addAttribute('summary','A long-lived 1.5 layer quasigeostrophic eddy on a beta plane, with Lagrangian particles'); 
ncfile.addAttribute('Conventions','CF-1.10, ACDD-1.3'); 

% Using this vocabulary for keywords as recommended by ACDD:
% https://gcmd.earthdata.nasa.gov/KeywordViewer/scheme/all?gtm_search=wave&gtm_scheme=all
ncfile.addAttribute('keywords','ocean currents, component process models, mesoscale eddies, potential vorticity, ocean waves, turbulence')


%% Recommended ACDD conventions
% https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3#Recommended
ncfile.addAttribute('id','10.5281/zenodo.8200055');
ncfile.addAttribute('naming_authority','Zenodo');
ncfile.addAttribute('product_version','1.0.2');

%%
newHistory = sprintf('%s: adding ACDD annotations to file.',datetime('now'));
ncfile.addAttribute('history','14-July-2024 version 1.0.2 created updating references to dataset and to paper | 06-Nov-2023 version 1.0.1 created updating some typographic errors in units, changing spelling of name, and adding some additional metadata | 31-Jul-2023 original verison 1.0.0 created');%suggested

%This is the only CF convention that is not also an ACDD convention...seems redundant with creator_institution below but whatever 
ncfile.addAttribute('institution','NorthWest Research Associates');
ncfile.addAttribute('processing_level','Unprocessed original model output'); %textual description 
ncfile.addAttribute('acknowledgment','J. Early''s work on this model run was supported by grant number 2049521 from the Physical Oceanography program of the United States National Science Foundation.  J. M. Lilly''s contribution to the NetCDF data product was carried out with the support of NSF award 2220291.');
ncfile.addAttribute('license','Creative Commons Attribution 4.0 International, https://creativecommons.org/licenses/by/4.0/');
ncfile.addAttribute('standard_name_vocabulary','CF Standard Name Table');

%in case you change in the future: 
%ncfile.addAttribute('date_modified',datestr(now,1));%suggested
ncfile.addAttribute('creator_name','Jeffrey J. Early');
ncfile.addAttribute('creator_email','jearly@nwra.com');
ncfile.addAttribute('creator_type','person');
ncfile.addAttribute('creator_url','https://jeffreyearly.com');
ncfile.addAttribute('creator_institution','NorthWest Research Associates'); %suggested
ncfile.addAttribute('project','Global Eddy-Driven Transport Estimated from in situ Lagrangian Observations (NSF award 2048552);  A Coordinate-Free Framework for Improving Eddy Parameterizations (NSF award 2220291)');
ncfile.addAttribute('publisher_name','Jonathan M. Lilly') %The name of the person (or other entity specified by the publisher_type attribute) responsible for publishing the data file or product to users, with its current metadata and format.
ncfile.addAttribute('publisher_email','jmlilly@psi.edu');
ncfile.addAttribute('publisher_url','http://www.jmlilly.net')
ncfile.addAttribute('publisher_institution','Planetary Science Institute') %suggested

%skipping these unless you want to acknowledge another contributor 
%ncfile.addAttribute('contributor_name',''); %suggested
%ncfile.addAttribute('contributor_role',''); %suggested

addGeospatialAnnotations(ncfile);

ncfile.addAttribute('platform','Models'); %suggested
ncfile.addAttribute('platform_vocabulary','GCMD, https://gcmd.earthdata.nasa.gov/static/kms/'); %suggested
%From https://docs.unidata.ucar.edu/netcdf-java/4.6/userguide/metadata/DataDiscoveryAttConvention.html
%The "cdm_data_type" attribute gives the THREDDS data type appropriate for this dataset. E.g., "Grid", "Image", "Station", "Trajectory", "Radial".
ncfile.addAttribute('cdm_data_type','Grid'); 

newReference = "Lilly, Feske, Fox-Kemper, and Early (2024).  Integral theorems for the gradient of a vector field, with a fluid dynamical application.  Proceedings of the Royal Society of London, Series A.  480 (2293): 20230550, 1–30.  doi 10.1098/rspa.2023.0550.";
appendStringAttribute(ncfile,'references',newReference);

ncfile.close();
end


function appendStringAttribute(ncfile,attname,string)
if isKey(ncfile.attributes,attname)
    existingStrings = ncfile.attributes(attname);
    existingStrings =[existingStrings,string];
else
    existingStrings = string;
end
ncfile.addAttribute(attname,existingStrings);
end

function addGeospatialAnnotations(ncfile)
% Let's say we have some model simulation with domain size (Lx, Ly) and we
% are going to pretend the box is located at its midpoint at lat0, lon0.
lat0 = ncfile.readVariables('latitude');
lon0 = 0;
Lx = ncfile.readVariables('Lx');
Ly = ncfile.readVariables('Ly');

[x0,y0] = LatitudeLongitudeToTransverseMercator(lat0,lon0,lon0=lon0);
xll = x0 - Lx/2;
yll = y0 - Ly/2;
xur = x0 + Lx/2;
yur = y0 + Ly/2;

[latMin,lonMin] = TransverseMercatorToLatitudeLongitude(xll,yll,lon0=lon0);
[latMax,lonMax] = TransverseMercatorToLatitudeLongitude(xur,yur,lon0=lon0);

%geospatial_bounds, geospatial_bounds_crs, geospatial_bounds_vertical_crs
ncfile.addAttribute('geospatial_lat_min',latMin);%XXX
ncfile.addAttribute('geospatial_lat_max',latMax);%XXX
ncfile.addAttribute('geospatial_lat_units','degree_north');
ncfile.addAttribute('geospatial_lon_min',lonMin);%XXX
ncfile.addAttribute('geospatial_lon_max',lonMax);%XXX
ncfile.addAttribute('geospatial_lon_units','degree_east');

t = ncfile.readVariables('t');


ncfile.addAttribute('time_coverage_start',string(datetime(t(1),ConvertFrom='posixtime',Format='yyyy-MM-dd HH:mm:ss.SSS')));%XXX
ncfile.addAttribute('time_coverage_end',string(datetime(t(end),ConvertFrom='posixtime',Format='yyyy-MM-dd HH:mm:ss.SSS')));%XXX
if length(t)>1
    ncfile.addAttribute('time_coverage_resolution',sprintf('PT%dS',t(2)-t(1)));%XXX
end
end

function [x,y] = LatitudeLongitudeToTransverseMercator(lat, lon, options)
arguments
    lat (:,1) double {mustBeNumeric,mustBeReal}
    lon (:,1) double {mustBeNumeric,mustBeReal}
    options.lon0 (1,1) double {mustBeNumeric,mustBeReal}
    options.k0 (1,1) double {mustBeNumeric,mustBeReal} = 0.9996;
end
k0 = options.k0;
if ~isfield(options,'lon0')
    lon0 = min(lon) + (max(lon)-min(lon))/2;
else
    lon0 = options.lon0;
end

% These are the *defined* values for WGS84
WGS84a=6378137;
WGS84invf=298.257223563;

% Convert from degrees to radians
phi = lat*pi/180;

% Compute a few trig functions that we'll be using a lot
s = sin(phi);
c = cos(phi);
t = tan(phi);
s2 = s.*s;
c2 = c.*c;
t2 = t.*t;

% Compute v and e2 -- pieces of this could be precomputed and #define
f = 1 / WGS84invf;
e2 = f*(2 - f);
v = WGS84a ./ sqrt( 1 - e2*s2);
e2 = e2 / (1 - e2); % From this point forward e2 will actually be (e^prime)^2
e2c2 = e2*c2;

deltaLambda = (lon - lon0)*pi/180;
d2c2 = deltaLambda.*deltaLambda.*c2;

% Terms to compute x.
T7 = 1 - t2 + e2c2;
T8 = 5 +t2.*(t2- 18) + e2c2.*(14 - 58*t2); % + 13.*e4*c4 + 4.*e6*c6 - 64.*t2*e4*c4 - 24.*t2*e6*c6;
T9 = 61 - t2.*(479 - t2.*(179 - t2));

x = k0 .* v .* c .* deltaLambda .* (1 + (d2c2/6).*( T7 + (d2c2/20).*( T8 + (d2c2/42).*T9 )));

% Terms to compute y.
T3 = 5 - t2 + e2c2.*(9 + 4*e2c2);
T4 = 61 - t2.*(58 - t2) + 270*e2c2 - 330*t2.*e2c2; % + 445.*e4*c4 + 324.*e6*c6 - 680.*t2*e4*c4 + 88.*e8*c8 - 600.*t2*e6*c6 - 192.*t2*e8*c8;
T5 = 1385 - t2.*(3111 - t2.*(543 - t2));

y = k0 * MeridionalArcPROJ4(phi) + (k0 * v .* s .* c / 2) .* deltaLambda.*deltaLambda .* (1 + (d2c2/12).*( T3 + (d2c2/30).*( T4 + (d2c2/56.).*T5)));
end

function [lat,lon] = TransverseMercatorToLatitudeLongitude(x, y, options)
arguments
    x (:,1) double {mustBeNumeric,mustBeReal}
    y (:,1) double {mustBeNumeric,mustBeReal}
    options.lon0 (1,1) double {mustBeNumeric,mustBeReal}
    options.k0 (1,1) double {mustBeNumeric,mustBeReal} = 0.9996;
end
lon0 = options.lon0;
k0 = options.k0;

% These are the *defined* values for WGS84
WGS84a=6378137;
WGS84invf=298.257223563;

phi = InverseMeridionalArcPROJ4(y/k0);

if ( abs( phi ) >= 2*pi )
    if y<0
        lat = -90;
    else
        lat = 90;
    end
    lon = 0;
    return;
end

% Compute a few trig functions that we'll be using a lot
s = sin(phi); s2 = s.*s;
c = cos(phi); c2 = c.*c;
t = tan(phi); t2 = t.*t;

% Compute v and e2 -- pieces of this could be precomputed and #define
f = 1 / WGS84invf;
e2 = f*(2 - f);
v = WGS84a ./ sqrt( 1 - e2*s2);
e2 = e2 / (1 - e2); % From this point forward e2 will actually be (e^prime)^2
e2c2 = e2*c2;

T11 = 5 + 3*t2 + e2c2.*(1 - 4*e2c2 - 9*t2);
T12 = 61 + t2.*(90 - 25.*e2c2 + 45*t2) + 46*e2c2;
T13 = 1385 + t2.*(3633. + t2.*(4095 + t2.*1575.));

d = x./(v*k0);
d2 = d.*d;

phi = phi - 0.5*(1 + e2c2).*t.*d2.*(1 - (d2/12).*( T11 - (d2/30).*( T12 - (d2/56).*T13)));
lat = phi*180./pi;

T15 = 1 + 2*t2 + e2c2;
T16 = 5 + t2.*(28 + 24*t2 + 8*e2c2) + 6*e2c2;
T17 = 61 + t2 .* (662 + t2 .* (1320 + 720 * t2));

lon = lon0 + 180/pi*(d./c).*(1 - (d2/6).*( T15 - (d2/20).*( T16 - (d2/42).*T17 )));
end
