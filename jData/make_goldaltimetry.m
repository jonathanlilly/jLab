function[varargout]=make_goldaltimetry(varargin)
%MAKE_GOLDALTIMETRY  Create a synthetic altimeter dataset and gridded statistics.
%
%   MAKE_GOLDALTIMETRY creates two files for used in observing system
%   simulation experiments of Jason-class satellite altimetry. 
% 
%   GoldAlongTrack is a synthetic Jason-class alongtrack altimeter dataset 
%   from GOLD model output, specifically a yearlong no-tides simulation.
%
%   GoldStatistics contains statistics of the GOLD no-tides simulation 
%   within each Jason altimeter cycle.
%
%   This file is provided for completeness, although it will not be 
%   possible for you to run it without access to the original model files.
%
%   Usage: make_goldaltimetry --create
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2024 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--create')
    make_goldalongtrack;
    %make_goldstatistics;
    return
elseif strcmpi(varargin{1}, '--f')
    make_goldalongtrack_figure
    return
end
 
function[]=make_goldalongtrack_figure


%/*************************************************************************
%reading in goldalongtrack data
filename='GoldAlongTrack.nc';
lat=ncread(filename,'lat');
lon=ncread(filename,'lon');
cycle_time=ncread(filename,'cycle_time')+datenum(1950,1,1);
%mss=ncread(filename,'mss');
%time=ncread(filename,'time')+datenum(1950,1,1);
sla=ncread(filename,'sla_inst')*100;
mss=ncread(filename,'mss')*100;
mss=mss-mean(mss,"all","omitmissing");
atd=ncread(filename,'atd');
%\*************************************************************************


%/*************************************************************************
sla_std=vstd(sla,3);

figure
%--------------------------------------------------------------------------
ax(1)=subplot(2,1,1);
%mat=log10(100*(1-sla_count./size(sla,3)));mat(mat==2)=nan;
jpcolor(1:size(sla,2),vmean(atd,2),mss)
caxis([-185 185]),vlines(127.5,'0.3k:'),axis tight,colormap(ax(1),flipud(crameri('broc')))
title('The GoldAlongTrack Dataset')
%--------------------------------------------------------------------------
ax(2)=subplot(2,1,2);
jpcolor(1:size(sla,2),vmean(atd,2),sla_std)
colorquant(0.5),vlines(127.5,'0.3k:'),axis tight,colormap(ax(2),flipud(crameri('davos')))
xlabel(['(Descending Tracks 1--127)\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,' ...
    'Track Number\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,(Ascending Tracks 128--254)'])
%--------------------------------------------------------------------------
ax=packfig(2,1,'rows');
for i=1:2
   axes(ax(i))
   %text(245,18,['(' char( real('a')+i-1) ')'])
   hc=colorbar('EastOutside');
   ylabel('Along-Track Distance (Mm)')
   if i==2
       hc.Label.String='SSH Standard Deviation (cm)';
       hc.Ticks=[0:5:50];
   elseif i==1
       hc.Label.String='Mean SSH (cm)';
   end
   %plot([92+0*1i 92+2*1i],'w','linewidth',2)
   %plot([92+0*1i 92+2*1i],'k','linewidth',1.5)
   %plot([92+maxmax(atd)*1i 92+(maxmax(atd)-2)*1i],'w','linewidth',2)
   %plot([92+maxmax(atd)*1i 92+(maxmax(atd)-2)*1i],'k','linewidth',1.5)
   set(gcf,'color','w');set(gca,'color',0.6*[1 1 1]);set(gcf, 'InvertHardcopy', 'off')
   vlines(127.5,'0.3k')
   xticks([25:25:300]),yticks([2:2:18])%,fixlabels([0,-1])
   %yticks([0:.5:3.5]),
   posax=get(gca,'Position');
   pos=hc.Position;
   hc.Position=[pos(1) pos(2)+0.5*pos(3) pos(3)/2 pos(4)-pos(3)];
   set(gca,'Position',posax)%thin x colorbar
   set(gca,'TickLen',get(gca,'TickLen')/2)
end
%--------------------------------------------------------------------------
fontsize 8 8 8 8
%set(gcf,'paperposition',[1 1 8.5 9]),
set(gcf,'paperposition',[1 1 7 5]),
jprint(pwd,'goldalongtrack','png','-r300')
%\*************************************************************************

function[]=make_goldalongtrack
%gold ssh interpolated onto tpjaos tracks

%/*************************************************************************
dir='/Users/lilly/Desktop/Dropbox/NetCDF/';

lat=ncread([dir 'JasonAlongTrack.nc'],'lat');
lon=ncread([dir 'JasonAlongTrack.nc'],'lon');
time_offset=ncread([dir 'JasonAlongTrack.nc'],'time_offset');
cycle_time=ncread([dir 'JasonAlongTrack.nc'],'cycle_time');

%figure,plot(lono,'.'),vlines(64.98)

%-------------------------------------------------------------------------
%the model lat,lon, and num
geofile='/Volumes/Thrills/Gold/NetCDF/ocean_geometry.nc';
lato=ncread(geofile,'geolat')';%Geographic longitudes of h-points
lono=ncread(geofile,'geolon')';%Geographic latitudes of h-points
numo=(datenum(2007,1,1):1/24:(datenum(2008,1,1)-1/24))';

%We will use a 12th degree grid for intermediate interpolation above 64.98 
%where the grid becomes non-plaid 
%lon1=[-180+1/24:1/12:180-1/24]';
lon1=(0+1/24:1/12:360-1/24)';
lat1=(65-1/24:1/12:67)';
tic;[dx,dy,index,bool,~]=sphereinterp(lato,lono,lat1,lon1,'parallel','periodic');toc;
%-------------------------------------------------------------------------
%the tpjaos lat, lon, and num
dir=['/Users/lilly/Desktop/Dropbox/NetCDF/'];
lat=ncread([dir 'JasonAlongTrack.nc'],'lat');
lon=ncread([dir 'JasonAlongTrack.nc'],'lon');
dnum=ncread([dir 'JasonAlongTrack.nc'],'time_offset');
num=ncread([dir 'JasonAlongTrack.nc'],'cycle_time')+datenum(1950,1,1);
[y,~,~]=datevec(num);
num=num(y==2007);
num=num(1:end-1);%omit the last point because the whole ten-day window doesn't fit
%\*************************************************************************

%/*************************************************************************
%Make synoptic version of GOLDALONGTRACK, interpolating within original fields
%I'm interpolating the model lato,lono onto jason tracks lat,lon
infile='/Volumes/Thrills/Gold/Model Archive/gold-notides/SSH_surface__0005_nt.nc';
%--------------------------------------------------------------------------
%local mean version
tic
mssho=zeros(size(lono,1),size(lono,2),length(num));
for i=1:length(num)
    i
    [~,a]=min(abs(num(i)+minmin(dnum)-numo));%first point
    numi=numo(a:a+238);%238=9.916*24;
    %(num(i)+minmin(dnum)-numi(1))*24
    %(num(i)+maxmax(dnum)-numi(end))*24

    ssho=ncread(infile,'SSH',[1 1 a],[inf inf 238]);%9.92 day window    
    ssho(abs(ssho)>1e3)=nan;
    ssho=permute(ssho,[2 1 3]);
    if i==2%fix a write error affecting a single slice in the second cycle
        ssho(:,:,151)=ssho(:,:,150)/2+ssho(:,:,152)/2;
    end
    %size(ssho)
    mssho(:,:,i)=vmean(ssho,3);
end    
toc

ssh=zeros(size(lon,1),size(lon,2),length(num));
for i=1:length(num)
    i
    latobool=lato(:,1)<64.98;
    latbool=lat>=65-1/24;

    %below lat 64.97, the model grid is plaid and I can use interplatlon 
    msshoi=mssho(:,:,i);
    sshi=interplatlon(lato(latobool,1),lono(1,:),[],msshoi(latobool,:),lat,lon);
    
    %above lat 64.97, the model grid is non-plaid and I use an intermediary
    ssh1=sphereinterp(dx,dy,index,bool,msshoi);
    sshi(latbool)=interplatlon(lat1,lon1,[],ssh1,lat(latbool),lon(latbool));
    ssh(:,:,i)=sshi;
end
ssh_mean=ssh;
toc %42 minutes
%--------------------------------------------------------------------------
%instantaneous version
tic
ssh=zeros(size(lon,1),size(lon,2),length(num));
num1=zeros(size(num));
for i=1:length(num)
    i
    [~,ii]=min(abs(numo-num(i)));%midpoint
    numi=numo([ii-5*24:ii+5*24-1]);
    num1(i)=mean(numi);

    ssho=ncread(infile,'SSH',[1 1 ii-5*24],[inf inf 10*24]);%10 day window    
    ssho(abs(ssho)>1e3)=nan;
    ssho=permute(ssho,[2 1 3]);
    
    %below lat 64.97, the model grid is plaid and I can use interplatlon 
    latobool=lato(:,1)<64.98;
    %anyany(num(i)+dnum>max(numi)|num(i)+dnum<min(numi))
    sshi=interplatlon(lato(latobool,1),lono(1,:),numi,ssho(latobool,:,:),lat,lon,num(i)+dnum);

    %above lat 64.97, the model grid is non-plaid and I use an intermediary
    ssh1=zeros(size(dx,1),size(dx,2),size(ssho,3));
    for k=1:size(ssho,3)
        ssh1(:,:,k)=sphereinterp(dx,dy,index,bool,ssho(:,:,k));
    end
    latbool=lat>=65-1/24;
    sshi(latbool)=interplatlon(lat1,lon1,numi,ssh1,lat(latbool),lon(latbool),num(i)+dnum(latbool));
    ssh(:,:,i)=sshi;
end
ssh_inst=ssh;
toc

%fill bad data points with along-track interpolation
for i=1:length(num)
    ssh_mean(:,:,i)=fillbad(ssh_mean(:,:,i),nan,inf);
    ssh_inst(:,:,i)=fillbad(ssh_inst(:,:,i),nan,inf);
end

%set to NaNs all points that do not correspond to oceans in JasonAlongTrack
stf=ncread([dir 'JasonAlongTrack.nc'],'stf');
for i=1:length(num)
    sshi=ssh_mean(:,:,i);
    sshi(stf~=0)=nan;
    ssh_mean(:,:,i)=sshi;

    sshi=ssh_inst(:,:,i);
    sshi(stf~=0)=nan;
    ssh_inst(:,:,i)=sshi;
end

%--------------------------------------------------------------------------
%mean sea surface ... interpolate from 1/4 degree model fields
msso=ncread([dir 'GoldOceanNT.nc'],'mss');%in version 1.0.0 the "NT" was incorrectly omitted
latmss=ncread([dir 'GoldOceanNT.nc'],'lat');
lonmss=ncread([dir 'GoldOceanNT.nc'],'lon');
mss=interplatlon(latmss,lonmss,[],msso',lat,lon);
%\*************************************************************************

%/*************************************************************************
%start saving the NetCDF file
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

writedir='/Users/lilly/Desktop/Dropbox/NetCDF/';
name_nc='GoldAlongTrack.nc';
ncfile = NetCDFFile([writedir name_nc]);

% Highly recommended ACDD conventions
% Descriptions of these four attributes are found here:
% https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3#Highly_Recommended
ncfile.addAttribute('title','GoldAlongTrack: A synthetic Jason-class alongtrack altimeter dataset from GOLD model output');
ncfile.addAttribute('summary','Model sea level anomalies and mean sea surface height, from a yearlong simulation with the GOLD general circulation model without tidal forcing performed at a nominal 1/8 degree resolution, sampled in the same way as JasonAlongTrack, a reformatted version of the Integrated Multi-Mission Ocean Altimeter Data for Climate Research Version 5.1.  Model sea level anomalies consist of a 3D array with dimensions of along-track direction by geographically sorted track number by cycle.');
ncfile.addAttribute('Conventions','CF-1.10, ACDD-1.3'); 

% Using this vocabulary for keywords as recommended by ACDD:
% https://gcmd.earthdata.nasa.gov/KeywordViewer/scheme/all?gtm_search=wave&gtm_scheme=all
ncfile.addAttribute('keywords','ocean general circulation models, ocean circulation, sea surface height, Jason-class altimeter');

% Recommended ACDD conventions
% https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3#Recommended
ncfile.addAttribute('id','10.5281/zenodo.11410611/GoldAlongTrack');
ncfile.addAttribute('naming_authority','Zenodo');
ncfile.addAttribute('product_version','1.0.1');
ncfile.addAttribute('history','19-June-2024 v. 1.0.1 fixed incorrect mean sea surface used in computing sea level anomaly | 01-June-2024 v. 1.0.0 initial version');

%This is the only CF convention that is not also an ACDD convention...seems redundant with creator_institution below but whatever 
ncfile.addAttribute('institution','University of Washington Applied Physics Laboratory');
ncfile.addAttribute('source','GOLD model as described in Simmons and Alford (2012)');
ncfile.addAttribute('processing_level','Interpolated from model grid to JasonAlongTrack tracks by J. M. Lilly using the jLab function INTERPLATLON');
ncfile.addAttribute('acknowledgment','The simulation on which this product is based was performed by H. Simmons with the support of ONR grant N00014-09-1-0399 and NSF grant OCE-0968838.  The creation of this NetCDF data product from the GOLD model output was carried out by J. M. Lilly with the support of NASA grant 80NSSC21K1823.');
ncfile.addAttribute('license','Creative Commons Attribution 4.0 International, https://creativecommons.org/licenses/by/4.0/');
ncfile.addAttribute('standard_name_vocabulary','CF Standard Name Table v79');
ncfile.addAttribute('comment','The organization of this dataset follows that of JasonAlongTrack, a reformatted version of the Brian Beckley and Richard Ray''s TPJAOS dataset, available from https://zenodo.org/doi/10.5281/zenodo.10055670.  Note that time is approximated by time_offset+cycle_time; see the comment under the time_offset variable for details.');

ncfile.addAttribute('date_created','2024-06-01T');
ncfile.addAttribute('date_modified','2024-06-19T');%suggested
ncfile.addAttribute('creator_name','Harper Simmons');
ncfile.addAttribute('creator_email','hsimmons@apl.washington.edu');
ncfile.addAttribute('creator_type','person');
ncfile.addAttribute('creator_url','https://www.apl.washington.edu/people/profile.php?last_name=Simmons&first_name=Harper');
ncfile.addAttribute('creator_institution','University of Washington Applied Physics Laboratory'); %suggested
ncfile.addAttribute('project','Near-inertial wave studies using historical mooring records and a high-resolution general circulation model (ONR grant N00014-09-1-0399), Collaborative Research: Representing internal-wave driven mixing in global ocean models (NSF grant OCE-0968838), Eddy dynamics from along-track altimetry (NASA grant 80NSSC21K1823)');
ncfile.addAttribute('publisher_name','J. M. Lilly');
ncfile.addAttribute('publisher_email','jmlilly@psi.edu');
ncfile.addAttribute('publisher_url','http://www.jmlilly.net');
ncfile.addAttribute('publisher_institution','Planetary Science Institute'); %suggested
ncfile.addAttribute('references','Simmons, H.L., and M.H. Alford (2012). Simulating the long-range swell of internal waves generated by ocean storms. Oceanography 25(2):30–41, http://dx.doi.org/10.5670/oceanog.2012.39.');%suggested

ncfile.addAttribute('platform','Models'); %suggested
ncfile.addAttribute('platform_vocabulary','GCMD, https://gcmd.earthdata.nasa.gov/static/kms/'); %suggested
%From https://docs.unidata.ucar.edu/netcdf-java/4.6/userguide/metadata/DataDiscoveryAttConvention.html
%The "cdm_data_type" attribute gives the THREDDS data type appropriate for this dataset. E.g., "Grid", "Image", "Station", "Trajectory", "Radial".
ncfile.addAttribute('cdm_data_type','Trajectory');
ncfile.addAttribute('featureType','trajectory');

ncfile.addAttribute('geospatial_lat_min',minmin(lat));
ncfile.addAttribute('geospatial_lat_max',maxmax(lat));
ncfile.addAttribute('geospatial_lon_min',minmin(lon));
ncfile.addAttribute('geospatial_lon_max',maxmax(lon));

ncfile.addAttribute('time_coverage_start','2007-01-08T18:46:04');%datestr(num(1))
ncfile.addAttribute('time_coverage_end','2007-12-21T19:54:35');%datestr(num(end))
ncfile.addAttribute('time_coverage_duration','P0001-00-00T00:00:00');%one year
ncfile.addAttribute('time_coverage_resolution','P0000-00-00T00:00:01');%one second

ncid = ncfile.ncid;
% Define the dimensions
dimitime = netcdf.defDim(ncid,'cycle_time',length(num1));
dimialong = netcdf.defDim(ncid,'along_track',size(lat,1));
dimitrack = netcdf.defDim(ncid,'track_number',size(lat,2));

cycle_time_ID = netcdf.defVar(ncid,'cycle_time','double',dimitime);
netcdf.putAtt(ncid,cycle_time_ID,'standard_name','time');
netcdf.putAtt(ncid,cycle_time_ID,'long_name','the mean of all valid observations times within each cycle');
netcdf.putAtt(ncid,cycle_time_ID,'units','days since 1950-01-01 00:00:00');
netcdf.putAtt(ncid,cycle_time_ID,'axis','T')
netcdf.putAtt(ncid,cycle_time_ID,'time_zone','UTC') ;
netcdf.putAtt(ncid,cycle_time_ID,'calendar','standard');
netcdf.putAtt(ncid,cycle_time_ID,'cell_methods','along_track: track_number: mean');%see 7.3.1 CF 
netcdf.putAtt(ncid,cycle_time_ID,'coverage_content_type','coordinate');

time_offset_ID = netcdf.defVar(ncid,'time_offset','double',[dimialong dimitrack]);
netcdf.putAtt(ncid,time_offset_ID,'standard_name','time');
netcdf.putAtt(ncid,time_offset_ID,'long_name','the mean difference, averaged over all cycles, between the time at each observation point and cycle_time');
netcdf.putAtt(ncid,time_offset_ID,'comment','For the ith along-track location, jth track, and kth cycle, time_offset(i,j)+cycle_time(k) approximates time(i,j,k) very closely, with an RMS error of 0.30 seconds and a maximum error of 3.8 seconds.');
netcdf.putAtt(ncid,time_offset_ID,'coordinates','lat lon');
netcdf.putAtt(ncid,time_offset_ID,'units','days');
netcdf.putAtt(ncid,time_offset_ID,'time_zone','UTC') ;
netcdf.putAtt(ncid,time_offset_ID,'calendar','standard');
netcdf.putAtt(ncid,time_offset_ID,'cell_methods','cycle_time: mean');%see 7.3.1 CF 
netcdf.putAtt(ncid,time_offset_ID,'coverage_content_type','referenceInformation');
netcdf.defVarFill(ncid,time_offset_ID,false,nan);

lat_ID = netcdf.defVar(ncid,'lat','double',[dimialong dimitrack]);
netcdf.putAtt(ncid,lat_ID,'standard_name','latitude');
netcdf.putAtt(ncid,lat_ID,'long_name','latitude');
netcdf.putAtt(ncid,lat_ID,'units','degrees_north');
netcdf.putAtt(ncid,lat_ID,'coverage_content_type','referenceInformation');

lon_ID = netcdf.defVar(ncid,'lon','double',[dimialong dimitrack]);
netcdf.putAtt(ncid,lon_ID,'standard_name','longitude');
netcdf.putAtt(ncid,lon_ID,'long_name','longitude');
netcdf.putAtt(ncid,lon_ID,'units','degrees_east');
netcdf.putAtt(ncid,lon_ID,'coverage_content_type','referenceInformation');

rev_ID = netcdf.defVar(ncid,'rev','int',[dimialong dimitrack]);
netcdf.putAtt(ncid,rev_ID,'long_name','reference orbit revolution number');
netcdf.putAtt(ncid,rev_ID,'coordinates','lat lon');
netcdf.putAtt(ncid,rev_ID,'valid_min',1);
netcdf.putAtt(ncid,rev_ID,'valid_max',127);
netcdf.putAtt(ncid,rev_ID,'comment','specifies the number of the revolution within a near 10-day repeat orbit')
netcdf.putAtt(ncid,rev_ID,'coverage_content_type','referenceInformation');

atd_ID = netcdf.defVar(ncid,'atd','double',[dimialong dimitrack]);
netcdf.putAtt(ncid,atd_ID,'long_name','great circle along-track distance from southernmost point');
netcdf.putAtt(ncid,atd_ID,'coordinates','lat lon');
netcdf.putAtt(ncid,atd_ID,'units','Mm')
netcdf.putAtt(ncid,atd_ID,'coverage_content_type','referenceInformation');
netcdf.putAtt(ncid,atd_ID,'comment','1 Mm = 1000 km')

dfc_ID = netcdf.defVar(ncid,'dfc','double',[dimialong dimitrack]);
netcdf.putAtt(ncid,dfc_ID,'long_name','along-track distance from coast');
netcdf.putAtt(ncid,dfc_ID,'coordinates','lat lon');
netcdf.putAtt(ncid,dfc_ID,'units','km')
netcdf.putAtt(ncid,dfc_ID,'coverage_content_type','auxillaryInformation');

stf_ID = netcdf.defVar(ncid,'stf','int',[dimialong dimitrack]);
netcdf.putAtt(ncid,stf_ID,'long_name','surface type flag');
netcdf.putAtt(ncid,stf_ID,'comment','0 - ocean; 1 and 3 -land; 2 - inland sea or lake');
netcdf.putAtt(ncid,stf_ID,'coordinates','lat lon');
netcdf.putAtt(ncid,stf_ID,'coverage_content_type','auxiliaryInformation');

depth_ID = netcdf.defVar(ncid,'depth','double',[dimialong dimitrack]);
netcdf.putAtt(ncid,depth_ID,'standard_name','sea_floor_depth_below_geoid');
netcdf.putAtt(ncid,depth_ID,'long_name','sea floor depth');
netcdf.putAtt(ncid,depth_ID,'coordinates','lat lon');
netcdf.putAtt(ncid,depth_ID,'source','1 arc-minute bathymetry/topography grid ETOPO1 (Amante and Eakins, 2009');
netcdf.putAtt(ncid,depth_ID,'units','km')
netcdf.putAtt(ncid,depth_ID,'positive','down')
netcdf.putAtt(ncid,depth_ID,'coverage_content_type','auxillaryInformation');

%--------------------------------------------------------------------------
mss_ID = netcdf.defVar(ncid,'mss','float',[dimialong dimitrack]);
netcdf.putAtt(ncid,mss_ID,'units','m');
netcdf.putAtt(ncid,mss_ID,'coordinates','lat lon');
netcdf.putAtt(ncid,mss_ID,'standard_name','sea_surface_height_above_geoid');
netcdf.putAtt(ncid,mss_ID,'long_name','GOLD mss sampled along JasonAlongTrack tracks');
netcdf.putAtt(ncid,mss_ID,'comment','model sea surface height interpolated onto a quarter degree grid, averaged over all hourly model snapshots within the yearlong simulation, and interpolated onto the JasonAlongTrack altimeter tracks');
netcdf.putAtt(ncid,mss_ID,'coverage_content_type','modelResult');
netcdf.putAtt(ncid,mss_ID,'cell_methods','time: mean');%see 7.3.4 CF 
netcdf.defVarFill(ncid,mss_ID,false,nan);
%--------------------------------------------------------------------------
ssh_mean_ID = netcdf.defVar(ncid,'sla_mean','float',[dimialong dimitrack dimitime]);
netcdf.putAtt(ncid,ssh_mean_ID,'units','m');
netcdf.putAtt(ncid,ssh_mean_ID,'coordinates','cycle_time lat lon');
netcdf.putAtt(ncid,ssh_mean_ID,'standard_name','sea_surface_height_above_geoid');
netcdf.putAtt(ncid,ssh_mean_ID,'long_name','cycle-averaged GOLD sea surface height anomaly sampled along JasonAlongTrack tracks');
netcdf.putAtt(ncid,ssh_mean_ID,'comment','sea surface height anomaly relative to the mean sea surface given in mss averaged over all hourly model snapshots within a 9.92 day window centered at each JasonAlongTrack cycle_time and interpolated onto the JasonAlongTrack altimeter tracks');
netcdf.putAtt(ncid,ssh_mean_ID,'coverage_content_type','modelResult');
netcdf.putAtt(ncid,ssh_mean_ID,'cell_methods','time: mean');%see 7.3.4 CF 
netcdf.defVarFill(ncid,ssh_mean_ID,false,nan);
%--------------------------------------------------------------------------
ssh_inst_ID = netcdf.defVar(ncid,'sla_inst','float',[dimialong dimitrack dimitime]);
netcdf.putAtt(ncid,ssh_inst_ID,'units','m');
netcdf.putAtt(ncid,ssh_inst_ID,'coordinates','cycle_time lat lon');
netcdf.putAtt(ncid,ssh_inst_ID,'standard_name','sea_surface_height_above_geoid');
netcdf.putAtt(ncid,ssh_inst_ID,'long_name','instantaneous GOLD sea surface height anomaly sampled along JasonAlongTrack tracks');
netcdf.putAtt(ncid,ssh_inst_ID,'comment','sea surface height anomaly relative to the mean sea surface given in mss interpolated onto the JasonAlongTrack altimeter tracks at the approximate JasonAlongTrack sampling times time_offset+cycle_time');
netcdf.putAtt(ncid,ssh_inst_ID,'coverage_content_type','modelResult');
netcdf.putAtt(ncid,ssh_inst_ID,'cell_methods','time: mean');%see 7.3.4 CF 
netcdf.defVarFill(ncid,ssh_inst_ID,false,nan);
%--------------------------------------------------------------------------
netcdf.endDef(ncid);
%--------------------------------------------------------------------------
netcdf.putVar(ncid,cycle_time_ID,num1-datenum(1950,1,1));
netcdf.putVar(ncid,time_offset_ID,dnum);
netcdf.putVar(ncid,lat_ID,lat);
netcdf.putVar(ncid,lon_ID,lon);
netcdf.putVar(ncid,rev_ID,ncread([dir 'JasonAlongTrack.nc'],'rev'));
netcdf.putVar(ncid,atd_ID,ncread([dir 'JasonAlongTrack.nc'],'atd'));
netcdf.putVar(ncid,dfc_ID,ncread([dir 'JasonAlongTrack.nc'],'dfc'));
netcdf.putVar(ncid,stf_ID,ncread([dir 'JasonAlongTrack.nc'],'stf'));
netcdf.putVar(ncid,depth_ID,ncread([dir 'JasonAlongTrack.nc'],'depth'));
netcdf.putVar(ncid,mss_ID,mss);
netcdf.putVar(ncid,ssh_mean_ID,ssh_mean-mss);%subtract the mss
netcdf.putVar(ncid,ssh_inst_ID,ssh_inst-mss);%subtract the mss
netcdf.close(ncid)
%\*************************************************************************


function[]=make_goldstatistics
%gold ssh on quarter-degree grid time-averaged within JasonAlongTrack cycles ...
%this only takes about 10 minutes

%/*************************************************************************
%model time
numo=[datenum(2007,1,1):1/24:(datenum(2008,1,1)-1/24)]';
%the tpjaos lat, lon, and num
dir=['/Users/lilly/Desktop/Dropbox/NetCDF/'];
dnum=ncread([dir 'JasonAlongTrack.nc'],'time_offset');
num=ncread([dir 'JasonAlongTrack.nc'],'cycle_time')+datenum(1950,1,1);
[y,~,~]=datevec(num);
num=num(y==2007);
num=num(1:end-1);%omit the last point because the whole ten-day window doesn't fit
%-------------------------------------------------------------------------

%median(diff(num),'omitnan') %9.915643933054525

writedir=['/Users/lilly/Desktop/Dropbox/NetCDF/'];
name_nc='GoldOceanNT.nc';

%--------------------------------------------------------------------------
%compute local means and standard deviations
tic
lat=ncread([writedir name_nc],'latitude');
lon=ncread([writedir name_nc],'longitude');
numo=ncread([writedir name_nc],'time')+datenum(1950,1,1);

sla_mean=zeros(length(lon),length(lat),length(num));
sla_std=zeros(length(lon),length(lat),length(num));
sla_snap=zeros(length(lon),length(lat),length(num));
mss=ncread([writedir name_nc],'mss');

meannum=zeros(size(num));
tic
for i=1:length(num)
    i
    [~,a]=min(abs(num(i)+minmin(dnum)-numo));%first point
    meannum(i)=mean(numo(a:a+238));%238=9.916*24;
    %numa(i)=numo(a);
    %(num(i)+minmin(dnum)-numi(1))*24
    %(num(i)+maxmax(dnum)-numi(end))*24

    sla_snap(:,:,i)=ncread([writedir name_nc],'ssh',[1 1 a+238/2],[inf inf 1])-mss;%subtract the mss 

    %meannum(i)==numo(a+238/2) %mean time is the same as center time
    ssho=ncread([writedir name_nc],'ssh',[1 1 a],[inf inf 238]);%9.92 day window    
    sla_mean(:,:,i)=vmean(ssho-mss,3);%subtract the mss
    sla_std(:,:,i)=vstd(ssho-mss,3);%subtract the mss
end    
toc

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

writedir='/Users/lilly/Desktop/Dropbox/NetCDF/';
name='GoldStatistics.nc';
ncfile = NetCDFFile([writedir name]);

% Highly recommended ACDD conventions
% Descriptions of these four attributes are found here:
% https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3#Highly_Recommended
ncfile.addAttribute('title','GoldStatistics: Statistics of the GOLD no-tides ocean simulation within each Jason altimeter cycle');
ncfile.addAttribute('summary','Sea level anomalies from a yearlong simulation with the GOLD general circulation model without tidal forcing, performed at a nominal 1/8 degree resolution and interpolated from the model grid to a uniform quarter-degree grid, then sampled and averaged over each Jason-class altimeter cycle occuring within the time span of the simulation');
ncfile.addAttribute('Conventions','CF-1.10, ACDD-1.3'); 

% Using this vocabulary for keywords as recommended by ACDD:
% https://gcmd.earthdata.nasa.gov/KeywordViewer/scheme/all?gtm_search=wave&gtm_scheme=all
ncfile.addAttribute('keywords','ocean general circulation models, ocean circulation, sea surface height, Jason-class altimeter');

% Recommended ACDD conventions
% https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3#Recommended
ncfile.addAttribute('id','10.5281/zenodo.11410611/GoldStatistics');
ncfile.addAttribute('naming_authority','Zenodo');
ncfile.addAttribute('product_version','1.0.0');
ncfile.addAttribute('history','01-June-2024 v. 1.0.0 initial version');

%This is the only CF convention that is not also an ACDD convention...seems redundant with creator_institution below but whatever 
ncfile.addAttribute('institution','University of Washington Applied Physics Laboratory');
ncfile.addAttribute('source','GOLD model as described in Simmons and Alford (2012)');
ncfile.addAttribute('processing_level','Time slices, averages, and standard deviations from the GoldOceanNT dataset over time intervals matching the 9.92 day repeat cycles within the GoldAlongTrack synthetic along-track altimeter dataset');
ncfile.addAttribute('acknowledgment','The simulation on which this product is based was performed by H. Simmons with the support of ONR grant N00014-09-1-0399 and NSF grant OCE-0968838.  The creation of this NetCDF data product from the GOLD model output was carried out by J. M. Lilly with the support of NASA grant 80NSSC21K1823.');
ncfile.addAttribute('license','Creative Commons Attribution 4.0 International, https://creativecommons.org/licenses/by/4.0/');
ncfile.addAttribute('standard_name_vocabulary','CF Standard Name Table v79');
ncfile.addAttribute('comment','This dataset is intended to be used in conjunction with GoldAlongTrack, in which the GOLD no-tides simulation is sampled in the same way as the JasonAlongTrack along-track altimeter dataset.');

ncfile.addAttribute('date_created','2024-06-01T');
% %in case you change in the future: 
% %ncfile.addAttribute('date_modified',datestr(now,1));%suggested
ncfile.addAttribute('creator_name','Harper Simmons');
ncfile.addAttribute('creator_email','hsimmons@apl.washington.edu');
ncfile.addAttribute('creator_type','person');
ncfile.addAttribute('creator_url','https://www.apl.washington.edu/people/profile.php?last_name=Simmons&first_name=Harper');
ncfile.addAttribute('creator_institution','University of Washington Applied Physics Laboratory'); %suggested
ncfile.addAttribute('project','Near-inertial wave studies using historical mooring records and a high-resolution general circulation model (ONR grant N00014-09-1-0399), Collaborative Research: Representing internal-wave driven mixing in global ocean models (NSF grant OCE-0968838), Eddy dynamics from along-track altimetry (NASA grant 80NSSC21K1823)');
ncfile.addAttribute('publisher_name','J. M. Lilly');
ncfile.addAttribute('publisher_email','jmlilly@psi.edu');
ncfile.addAttribute('publisher_url','http://www.jmlilly.net');
ncfile.addAttribute('publisher_institution','Planetary Science Institute'); %suggested
ncfile.addAttribute('references','Simmons, H.L., and M.H. Alford (2012). Simulating the long-range swell of internal waves generated by ocean storms. Oceanography 25(2):30–41, http://dx.doi.org/10.5670/oceanog.2012.39.');%suggested

ncfile.addAttribute('platform','Models'); %suggested
ncfile.addAttribute('platform_vocabulary','GCMD, https://gcmd.earthdata.nasa.gov/static/kms/'); %suggested
%From https://docs.unidata.ucar.edu/netcdf-java/4.6/userguide/metadata/DataDiscoveryAttConvention.html
%The "cdm_data_type" attribute gives the THREDDS data type appropriate for this dataset. E.g., "Grid", "Image", "Station", "Trajectory", "Radial".
ncfile.addAttribute('cdm_data_type','Grid');
ncfile.addAttribute('featureType','grid');

ncfile.addAttribute('geospatial_lat_min',min(lat));
ncfile.addAttribute('geospatial_lat_max',max(lat));
ncfile.addAttribute('geospatial_lat_resolution','0.25 degree'); %suggested
ncfile.addAttribute('geospatial_lon_min',min(lon));
ncfile.addAttribute('geospatial_lon_max',max(lon));
ncfile.addAttribute('geospatial_lon_resolution','0.25 degree');  %suggested

ncfile.addAttribute('time_coverage_start','2007-01-08T18:46:04');%datestr(num(1))
ncfile.addAttribute('time_coverage_end','2007-12-21T19:54:35');%datestr(num(end))
ncfile.addAttribute('time_coverage_duration','P0001-00-00T00:00:00');%one year
ncfile.addAttribute('time_coverage_resolution','P0000-00-00T00:00:01');%one second

%Notes to self: these are not needed because they are assumed by default
%netcdf.putAtt(ncid,varid,'geospatial_lat_units','degree_north');
%netcdf.putAtt(ncid,varid,'geospatial_lon_units','degree_east');

ncid = ncfile.ncid;
% Define the dimensions ... in order of time, depth, lat, lon
dimitime = netcdf.defDim(ncid,'time',length(meannum));
dimilat = netcdf.defDim(ncid,'lat',length(lat));
dimilon = netcdf.defDim(ncid,'lon',length(lon));

varid = netcdf.defVar(ncid,'time','NC_DOUBLE',dimitime);
netcdf.putAtt(ncid,varid,'standard_name','time');
netcdf.putAtt(ncid,varid,'long_name','time of model output');
netcdf.putAtt(ncid,varid,'comment','these times match the center times of the Jason-class altimetric cycles, as documented in the JasonAlongTrack dataset, during the year 2007');
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

varid = netcdf.defVar(ncid,'sla_mid','NC_FLOAT',[dimilon dimilat dimitime]);
netcdf.putAtt(ncid,varid,'coordinates','lon lat');
netcdf.putAtt(ncid,varid,'standard_name','sea_surface_height_above_geoid');
netcdf.putAtt(ncid,varid,'long_name','instantaneous model sea surface height anomaly relative to mss');
netcdf.putAtt(ncid,varid,'comment','snapshots from the hourly model slice closest to each time, representing the middle of each altimeter cycle');
netcdf.putAtt(ncid,varid,'long_name','model sea surface height');
netcdf.putAtt(ncid,varid,'units','m');
netcdf.putAtt(ncid,varid,'coverage_content_type','modelResult');
netcdf.defVarFill(ncid,varid,false,nan);

varid = netcdf.defVar(ncid,'sla_mean','NC_FLOAT',[dimilon dimilat dimitime]);
netcdf.putAtt(ncid,varid,'coordinates','lon lat');
netcdf.putAtt(ncid,varid,'standard_name','sea_surface_height_above_geoid');
netcdf.putAtt(ncid,varid,'long_name','time-averaged model sea surface height anomaly relative to mss');
netcdf.putAtt(ncid,varid,'comment','the time average is taken over all hourly slices within a 9.92 day (238 hour) cycle centered on each time');
netcdf.putAtt(ncid,varid,'units','m');
netcdf.putAtt(ncid,varid,'cell_methods','time: mean');
netcdf.putAtt(ncid,varid,'coverage_content_type','modelResult');
netcdf.defVarFill(ncid,varid,false,nan);

varid = netcdf.defVar(ncid,'sla_std','NC_FLOAT',[dimilon dimilat dimitime]);
netcdf.putAtt(ncid,varid,'coordinates','lon lat');
netcdf.putAtt(ncid,varid,'standard_name','sea_surface_height_above_geoid');
netcdf.putAtt(ncid,varid,'long_name','standard deviation of model sea surface height anomaly relative to mss');
netcdf.putAtt(ncid,varid,'comment','the standard deviation is computed over all hourly slices within a 9.92 day (238 hour) cycle centered on each time');
netcdf.putAtt(ncid,varid,'units','m');
netcdf.putAtt(ncid,varid,'cell_methods','time: mean');
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

netcdf.endDef(ncid);
%--------------------------------------------------------------------------
%writing the variables
ncwrite([writedir name],'lat',lat);
ncwrite([writedir name],'lon',lon);
ncwrite([writedir name],'time',meannum-datenum(1950,1,1));

ncwrite([writedir name],'sla_mid',sla_snap);
ncwrite([writedir name],'sla_mean',sla_mean);
ncwrite([writedir name],'sla_std',sla_std);
ncwrite([writedir name],'mss',mss);

netcdf.close(ncid)
%\*************************************************************************

   
