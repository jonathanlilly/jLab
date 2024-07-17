function[varargout]=make_goldocean(varargin)
%MAKE_GOLDOCEAN  Create quarter-degree gridded files from GOLD simulations.
%
%   MAKE_GOLDOCEAN creates NetCDF files of hourly sea surface height from 
%   two different yearlong simulations using the Generalized Ocean Layer 
%   Dynamics (GOLD) model run by Harper Simmons, one that includes tidal 
%   forcing and one that does not.
% 
%   The resulting files, GoldOcean.nc and GoldOceanNT.nc, are available on
%   Zenodo at https://zenodo.org/doi/10.5281/zenodo.11396564 and 
%   https://zenodo.org/doi/10.5281/zenodo.11396570 respectively.
%
%   The SSH fields are interpolated from the model grid of nominal 1/8 
%   resolution to a uniform 1/4 degree grid, the same grid used for the 
%   AVISO/CMEMS gridded altimetry products.
%
%   This file is provided for completeness, although it will not be 
%   possible to run it without access to the original model files.
%
%   make_goldocean --f generates two sample figures. 
%
%   Usage: make_goldocean --create
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2024 J.M. Lilly --- type 'help jlab_license' for details
 
if nargin==0
    help make_goldocean
elseif nargin>0
    if strcmpi(varargin{1}, '--create')
        make_gold_gridded('notides')
        make_gold_gridded('tides')
        return
    elseif strcmpi(varargin{1}, '--f')
        make_gold_gridded_figures
        return
    end
end


function[]=make_gold_gridded(str)

%This is the AVISO quarter degree grid
lon=[-180+1/8:1/4:180-1/8]';
lat=[-90+1/8:1/4:90-1/8]';
%Hours since January 1, 2007 according to Harper
time=[datenum(2007,1,1):1/24:(datenum(2008,1,1)-1/24)]'-datenum(1950,1,1);
%datestr(time(1)+datenum(1950,1,1)  %ok

%Gathering the mapping information
geofile='/Volumes/Thrills/Gold/NetCDF/ocean_geometry.nc';
lato=ncread(geofile,'geolat')';%Geographic longitudes of h-points
lono=ncread(geofile,'geolon')';%Geographic latitudes of h-points
tic;[dx,dy,index,bool,~]=sphereinterp(lato,lono,lat,lon,'parallel','periodic');toc;
  
writedir='/Users/lilly/Desktop/Dropbox/NetCDF/';

if strcmpi(str(1:3),'not')
    infile='/Volumes/Thrills/Gold/Model Archive/gold-notides/SSH_surface__0005_nt.nc';
    name='GoldOceanNT.nc';
    ncfile = NetCDFFile([writedir name]);
    ncfile.addAttribute('title','GoldOceanNT: Hourly sea surface height from a yearlong global ocean simulation without tidal forcing');
    ncfile.addAttribute('summary','Sea surface height from a yearlong simulation with the GOLD general circulation model without tidal forcing, performed at a nominal 1/8 degree resolution and interpolated from the model grid to a uniform quarter-degree grid, the same grid used by the AVISO/CMEMS gridded altimetry products.');
    ncfile.addAttribute('id','10.5281/zenodo.11396570');
elseif strcmpi(str(1:3),'tid')
    infile='/Volumes/Thrills/Gold/Model Archive/gold-withtides/SSH_surface__0005_wt.nc';
    name='GoldOcean.nc';
    ncfile = NetCDFFile([writedir name]);
    ncfile.addAttribute('title','GoldOcean: Hourly sea surface height from a yearlong global ocean simulation');
    ncfile.addAttribute('summary','Sea surface height from a yearlong simulation with the GOLD general circulation model, performed at a nominal 1/8 degree resolution and interpolated from the model grid to a uniform quarter-degree grid, the same grid used by the AVISO/CMEMS gridded altimetry products.')
    ncfile.addAttribute('id','10.5281/zenodo.11396564');
end

%--------------------------------------------------------------------------
% Define the global attributes

%Some useful links
%https://foundations.projectpythia.org/core/data-formats/netcdf-cf.html
%http://cfconventions.org/cf-conventions/cf-conventions.html
%https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3
%https://www.unidata.ucar.edu/software/udunits/udunits-2.2.28/udunits2.html#Database
%https://github.com/Unidata/netcdf4-python/issues/442  %no reference year 0

%coverage_content_type	An ISO 19115-1 code to indicate the source of the data 
%(image, thematicClassification, physicalMeasurement, auxiliaryInformation, 
%qualityInformation, referenceInformation, modelResult, or coordinate).
%https://ngdc.noaa.gov/wiki/index.php/ISO_19115_and_19115-2_CodeList_Dictionaries#MD_CoverageContentTypeCode

% Using this vocabulary for keywords as recommended by ACDD:
% https://gcmd.earthdata.nasa.gov/KeywordViewer/scheme/all?gtm_search=wave&gtm_scheme=all
ncfile.addAttribute('keywords','ocean general circulation models, ocean circulation, sea surface height');
ncfile.addAttribute('Conventions','CF-1.10, ACDD-1.3');

% Recommended ACDD conventions
% https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3#Recommended
ncfile.addAttribute('naming_authority','Zenodo');
ncfile.addAttribute('product_version','1.0.0');
ncfile.addAttribute('history','30-May-2024 v. 1.0.0 initial version');

%This is the only CF convention that is not also an ACDD convention...seems redundant with creator_institution below but whatever 
ncfile.addAttribute('institution','University of Washington Applied Physics Laboratory');
ncfile.addAttribute('source','GOLD model as described in Simmons and Alford (2012)');
ncfile.addAttribute('processing_level','Interpolated from model grid to uniform quarter-degree grid by J. M. Lilly using the jLab function SPHEREINTERP'); %textual description 
ncfile.addAttribute('acknowledgment','The simulation on which this product is based was performed by H. Simmons with the support of ONR grant N00014-09-1-0399 and NSF grant OCE-0968838.  The creation of this NetCDF data product from the GOLD model output was carried out by J. M. Lilly with the support of NASA grant 80NSSC21K1823.');
ncfile.addAttribute('license','Creative Commons Attribution 4.0 International, https://creativecommons.org/licenses/by/4.0/');
ncfile.addAttribute('standard_name_vocabulary','CF Standard Name Table v79');

ncfile.addAttribute('date_created','2024-05-29T');
% %in case you change in the future: 
% %ncfile.addAttribute('date_modified',datestr(now,1));%suggested
ncfile.addAttribute('creator_name','Harper Simmons');
ncfile.addAttribute('creator_email','hsimmons@apl.washington.edu');
ncfile.addAttribute('creator_type','person');
ncfile.addAttribute('creator_url','https://www.apl.washington.edu/people/profile.php?last_name=Simmons&first_name=Harper');
ncfile.addAttribute('creator_institution','University of Washington Applied Physics Laboratory'); %suggested
ncfile.addAttribute('project','Near-inertial wave studies using historical mooring records and a high-resolution general circulation model (ONR grant N00014-09-1-0399), Collaborative Research: Representing internal-wave driven mixing in global ocean models (NSF grant OCE-0968838), Eddy dynamics from along-track altimetry (NASA grant 80NSSC21K1823)');
ncfile.addAttribute('publisher_name','Jonathan M. Lilly') %The name of the person (or other entity specified by the publisher_type attribute) responsible for publishing the data file or product to users, with its current metadata and format.
ncfile.addAttribute('publisher_email','jmlilly@psi.edu');
ncfile.addAttribute('publisher_url','http://www.jmlilly.net')
ncfile.addAttribute('publisher_institution','Planetary Science Institute') %suggested
ncfile.addAttribute('references','Simmons, H.L., and M.H. Alford (2012). Simulating the long-range swell of internal waves generated by ocean storms. Oceanography 25(2):30–41, http://dx.doi.org/10.5670/oceanog.2012.39.');%suggested

ncfile.addAttribute('platform','Models'); %suggested
ncfile.addAttribute('platform_vocabulary','GCMD, https://gcmd.earthdata.nasa.gov/static/kms/'); %suggested
%From https://docs.unidata.ucar.edu/netcdf-java/4.6/userguide/metadata/DataDiscoveryAttConvention.html
%The "cdm_data_type" attribute gives the THREDDS data type appropriate for this dataset. E.g., "Grid", "Image", "Station", "Trajectory", "Radial".
ncfile.addAttribute('cdm_data_type','Grid'); 

ncfile.addAttribute('geospatial_lat_min',min(lat));
ncfile.addAttribute('geospatial_lat_max',max(lat));
ncfile.addAttribute('geospatial_lat_resolution','0.25 degree'); %suggested
ncfile.addAttribute('geospatial_lon_min',min(lon));
ncfile.addAttribute('geospatial_lon_max',max(lon));
ncfile.addAttribute('geospatial_lon_resolution','0.25 degree');  %suggested

ncfile.addAttribute('time_coverage_start','2007-01-01T00:00:00Z');
ncfile.addAttribute('time_coverage_end','2007-12-31T23:00:00Z');
ncfile.addAttribute('time_coverage_duration','P0001-00-00T00:00:00');%one year
ncfile.addAttribute('time_coverage_resolution','P0000-00-00T01:00:00');%one hour

%Notes to self: these are not needed because they are assumed by default
%netcdf.putAtt(ncid,varid,'geospatial_lat_units','degree_north');
%netcdf.putAtt(ncid,varid,'geospatial_lon_units','degree_east');

ncid = ncfile.ncid;
% Define the dimensions ... in order of time, depth, lat, lon
dimitime = netcdf.defDim(ncid,'time',length(time));
dimilat = netcdf.defDim(ncid,'lat',length(lat));
dimilon = netcdf.defDim(ncid,'lon',length(lon));

varid = netcdf.defVar(ncid,'time','NC_DOUBLE',dimitime);
netcdf.putAtt(ncid,varid,'standard_name','time');
netcdf.putAtt(ncid,varid,'long_name','time of model output');
netcdf.putAtt(ncid,varid,'units','days since 1950-01-01 00:00:00');
netcdf.putAtt(ncid,varid,'time_zone','UTC') ;
netcdf.putAtt(ncid,varid,'axis','T')
netcdf.putAtt(ncid,varid,'calendar','standard');
netcdf.putAtt(ncid,varid,'coverage_content_type','coordinate');

varid = netcdf.defVar(ncid,'lat','NC_DOUBLE',dimilat);
netcdf.putAtt(ncid,varid,'standard_name','latitude');
netcdf.putAtt(ncid,varid,'long_name','latitude');
netcdf.putAtt(ncid,varid,'units','degrees_north');
netcdf.putAtt(ncid,varid,'axis','Y');
netcdf.putAtt(ncid,varid,'coverage_content_type','coordinate');

varid = netcdf.defVar(ncid,'lon','NC_DOUBLE',dimilon);
netcdf.putAtt(ncid,varid,'standard_name','longitude');
netcdf.putAtt(ncid,varid,'long_name','longitude');
netcdf.putAtt(ncid,varid,'units','degrees_east');
netcdf.putAtt(ncid,varid,'axis','X');
netcdf.putAtt(ncid,varid,'coverage_content_type','coordinate');

varid = netcdf.defVar(ncid,'ssh','NC_FLOAT',[dimilon dimilat dimitime]);
netcdf.putAtt(ncid,varid,'coordinates','lon lat');
netcdf.putAtt(ncid,varid,'standard_name','sea_surface_height_above_geoid');
netcdf.putAtt(ncid,varid,'long_name','model sea surface height');
netcdf.putAtt(ncid,varid,'units','m');
netcdf.putAtt(ncid,varid,'coverage_content_type','modelResult');
netcdf.defVarFill(ncid,varid,false,nan);

varid = netcdf.defVar(ncid,'mss','NC_FLOAT',[dimilon dimilat]);
netcdf.putAtt(ncid,varid,'coordinates','lon lat');
netcdf.putAtt(ncid,varid,'standard_name','sea_surface_height_above_geoid');
netcdf.putAtt(ncid,varid,'long_name','model sea surface height averaged over all hourly time slices');
netcdf.putAtt(ncid,varid,'units','m');
netcdf.putAtt(ncid,varid,'cell_methods','time: mean');
netcdf.putAtt(ncid,varid,'coverage_content_type','modelResult');
netcdf.defVarFill(ncid,varid,false,nan);

%deflate seems to slow things down a lot
%netcdf.defVarDeflate(ncid,ssh_ID,true,true,1);%shuffle and deflate level 1
%netcdf.defVarDeflate(ncid,ssh_ID,true,true,2);%shuffle and deflate level 2

%recommended deflate level is 5 according to 
%http://climate-cms.wikis.unsw.edu.au/NetCDF_Compression_Tools
%https://www.unidata.ucar.edu/blogs/developer/entry/netcdf_compression

netcdf.endDef(ncid);
%--------------------------------------------------------------------------
%writing the variables
ncwrite([writedir name],'lat',lat);
ncwrite([writedir name],'lon',lon);
ncwrite([writedir name],'time',time);

%--------------------------------------------------------------------------
%sea surface height
tic
for i=1:length(time)
    i
    ssho=ncread(infile,'SSH',[1 1 i],[inf inf 1])';
    ssho(abs(ssho)>1e3)=nan;
    ssh=sphereinterp(dx,dy,index,bool,ssho)';
    ncwrite([writedir name],'ssh',ssh,[1 1 i]);
end
toc%That took 1.2 hours
%--------------------------------------------------------------------------
%ssh without tides
%after-the-fact hack for notides #456, which has a write error in the original file
if strcmpi(str(1:3),'not')
    ssh1=ncread([writedir name],'ssh',[1 1 455],[inf inf 1]);
    ssh2=ncread([writedir name],'ssh',[1 1 457],[inf inf 1]);
    ncwrite([writedir name],'ssh',ssh1/2+ssh2/2,[1 1 456]);
end

%--------------------------------------------------------------------------
%compute mean sea surface
mss=zeros(length(lon),length(lat));
for i=1:length(time)
    i
    mss=mss+ncread([writedir name],'ssh',[1 1 i],[inf inf 1]);
end
mss=mss./length(time);
ncwrite([writedir name],'mss',mss);
%--------------------------------------------------------------------------
ncfile.close();
% 
% ncid = netcdf.open('GoldOcean.nc','NC_WRITE')
% netcdf.redef(ncid)
% [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,1)
% netcdf.renameVar(ncid,1,'lat')
% [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,0)
% 
% [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,2)
% netcdf.renameVar(ncid,2,'lon')
% [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,2)
% 
% 
% ncid = netcdf.open('GoldOceanNT.nc','NC_WRITE')
% netcdf.redef(ncid)
% [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,1);
% netcdf.renameVar(ncid,1,'lat')
% [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,2);
% netcdf.renameVar(ncid,2,'lon')

function[]=make_gold_gridded_figures

filename={'GoldOcean','GoldOceanNT'};
for i=1:2
    lat = ncread([filename{i} '.nc'],'lat');
    lon = ncread([filename{i} '.nc'],'lon');
    ssh = ncread([filename{i} '.nc'],'ssh',[1 1 1],[inf inf 1]);
    mss = ncread([filename{i} '.nc'],'mss');

    figure
    jpcolor(lon,lat,100*(ssh-mss)')
    topoplot continents
    latratio(30)
    ylim([-77.5 90])
    h=colorbar('eastoutside');
    colormap(flipud(crameri('broc')))
    caxis([-25 25])
    title(['Snapshot of Sea Surface Height Anomaly from ' filename{i}])
    h.Label.String="Sea Surface Height Anomaly (cm)";

    set(gcf,'paperposition',[1/2 1 10 6.5])
    fontsize 16 12 12 12
    jprint(pwd,[filename{i} '_snapshot'])
end