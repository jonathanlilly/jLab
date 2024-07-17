function[]=make_gulfflow(varargin)
%MAKE_GULFFLOW  Create a gridded velocity dataset for the Gulf of Mexico.
%
%   GULFFLOW is a dataset containing space/time gridded sea surface 
%   velocities from all available surface drifter data from the Gulf of
%   Mexico, as described in 
%
%        Lilly, J. M. and P. Perez-Brunius (2021a). A gridded surface
%            current product for the Gulf of Mexico from consolidated  
%            drifter measurements.  Submitted to Earth Science System Data.
%
%   The dataset is available at 
%
%        Lilly, J. M. and P. Perez-Brunius (2021). GulfFlow: A gridded 
%            surface current product for the Gulf of Mexico from 
%            consolidated drifter measurements. [Data set]. Zenodo. 
%            https://doi.org/10.5281/zenodo.3978793
%
%   This file creates the gridded drifter dataset from the consolidated 
%   surface drifter datasets.  It is provided for completeness, although it 
%   will not be possible for you to recreate it because you do not have
%   access to the proprietary GulfDriftersAll.nc dataset is it based on.
%
%   The code for conversion from Matlab files to NetCDF was written by
%   Paula Garcia Carrillo of CICESE with minor changes by J. Lilly.
%
%   Usage: make_gulfflow --create
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2020--2021 J.M. Lilly --- type 'help jlab_license' for details

if nargin==0
    help make_gulfflow
elseif nargin>0
    if strcmpi(varargin{1}, '--create')
        gulfflow_matfiles
        gulfflow_onequarter_netcdf
        gulfflow_onetwelfth_netcdf
    elseif strcmpi(varargin{1}, '--matfiles')
        gulfflow_matfiles
    elseif strcmpi(varargin{1}, '--netcdf')
        gulfflow_onequarter_netcdf
        gulfflow_onetwelfth_netcdf
    elseif strcmpi(varargin{1}, '--models')
        gulfflow_models
    end
end


function[]=gulfflow_matfiles

%/*************************************************************************
basedir='/Users/lilly/Desktop/Dropbox/NetCDF/';

%ncload([basedir '/former/GulfDriftersAll_100.nc']);
%use GulfDriftersAll_100
%length(find(cell2col(time)>datenum(2018,8,1)&cell2col(time)<datenum(2019,8,1)))
%length(find(cell2col(time)>datenum(2017,8,1)&cell2col(time)<datenum(2018,8,1)))

%Monthly velocity with 50% overlap
ncload([basedir 'GulfDriftersAll.nc']);
use GulfDriftersAll

%--------------------------------------------------------------------------
%modify to separately account for drifters with different drogue statuses
n=length(time);
for i=1:length(time)
    i
    [ids{2*n+i},time{2*n+i},lat{2*n+i},lon{2*n+i},u{2*n+i},v{2*n+i},source(2*n+i)]=...
        vindex(ids{i},time{i},lat{i},lon{i},u{i},v{i},source(i)+30,isnan(drogue{i}),1);
    [ids{n+i},time{n+i},lat{n+i},lon{n+i},u{n+i},v{n+i},source(n+i)]=...
        vindex(ids{i},time{i},lat{i},lon{i},u{i},v{i},source(i)+15,drogue{i}==0,1);
    [ids{i},time{i},lat{i},lon{i},u{i},v{i},source(i)]=...
        vindex(ids{i},time{i},lat{i},lon{i},u{i},v{i},source(i),drogue{i}==1,1);
end

ExpName=GulfDriftersAll.exp_names;
ExpName(16:30,:)=GulfDriftersAll.exp_names;
ExpName(31:45,:)=GulfDriftersAll.exp_names;

for i=1:15
    ExpName(i,9:18)=   '   drogued';
    ExpName(i+15,9:18)=' undrogued';
    ExpName(i+30,9:18)='   unknown';
end
%--------------------------------------------------------------------------
%length(find(cell2col(time)>datenum(2018,8,1)&cell2col(time)<datenum(2019,8,1)))
%length(find(cell2col(time)>datenum(2017,8,1)&cell2col(time)<datenum(2018,8,1)))
%new version has about 1000 more measurements during the last year

%length(find(cell2col(time)>datenum(2019,6,1)))  %only 2587 past this date
%length(find(cell2col(time)>datenum(2018,6,1)))  %versus 67576 past this date
%length(find(cell2col(time)<datenum(1993,6,1)))  %versus 12300 in first year
%--------------------------------------------------------------------------
%648 is 24 x  27 years
%672 is 24 x  28 years

lono=[-99:1/4:-80.5];
lato=[18:1/4:31];
gulfflow=griddrifters(time,lat,lon,u,v,filled,source,lato,lono,1992,7,672);
gulfflow.ExpName=ExpName;
matsave gulfflow

lono=[-99-1/24:1/12:-80.5];
lato=[18-1/24:1/12:31];
gulfflow_onetwelfth=griddrifters(time,lat,lon,u,v,filled,source,lato,lono,1992,7,672);
gulfflow_onetwelfth.ExpName=ExpName;
matsave gulfflow_onetwelfth
%\*************************************************************************


function[]=gulfflow_onequarter_netcdf
%by Paula Garcia Carrillo with modifications by J.M. Lilly 
writedir=['/Users/lilly/Desktop/Dropbox/NetCDF/'];

load gulfflow
time=gulfflow.num; %jml
%timei=gulfflow.num; num=timei;
%time=(timei*86400)-(datenum(1900,1,1,0,0,0)*86400);%%seconds
lon=gulfflow.lon; lat=gulfflow.lat;u=gulfflow.u; v=gulfflow.v; 
um=vmean(gulfflow.u,3); vm=vmean(gulfflow.v,3);
epsuu=gulfflow.epsuu; epsvv=gulfflow.epsvv; epsuv=gulfflow.epsuv; 
eps2=vmean(gulfflow.epsuu+gulfflow.epsvv,3);%jml
count=gulfflow.counts; total_count=gulfflow.count;%jml
ExpName=gulfflow.ExpName; %ExpName=char(ExpName); ExpName=ExpName';

%name_nc=['GulfFlow_GriddedFields_onequarter.nc'];
name_nc=['GulfFlow_OneQuarter.nc'];
% netCDF
% Define the global attributes
ncid = netcdf.create([writedir name_nc],'NETCDF4');
varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid,varid,'site_code','Gulf of Mexico');
netcdf.putAtt(ncid,varid,'platform_code','drifters_gridded_velocity_fields');
netcdf.putAtt(ncid,varid,'title','GulfFlow gridded velocity product, 1/4 deg');
netcdf.putAtt(ncid,varid,'product_version','1.1.0');
netcdf.putAtt(ncid,varid,'summary','gridded surface velocity product for the Gulf of Mexico derived from surface drifters, 1/4 deg resolution');
netcdf.putAtt(ncid,varid,'project','Implementacion de redes de observacion oceanograficas (fisicas geoquimicas, ecologicas) para la generacion de escenarios ante posibles contingencias relacionadas a la exploracion y produccion de hidrocarburos en aguas profundas del Golfo de Mexico Consorcio de Investigacion del Golfo de Mexico (CIGoM)');
netcdf.putAtt(ncid,varid,'naming_authority','cigom.org');
netcdf.putAtt(ncid,varid,'id','gridded_surface_velocity_fields');
netcdf.putAtt(ncid,varid,'source','drifters');
netcdf.putAtt(ncid,varid,'principal_investigator','Paula Perez-Brunius');
netcdf.putAtt(ncid,varid,'principal_investigator_email','brunius@cicese.mx');
netcdf.putAtt(ncid,varid,'principal_investigator_url','https://giola.cicese.mx');
netcdf.putAtt(ncid,varid,'institution','Centro de Investigación Científica y de Educación Superior de Ensenada, CICESE');
netcdf.putAtt(ncid,varid,'area','Gulf of Mexico');
netcdf.putAtt(ncid,varid,'geoespatial_lat_min',min(lat));
netcdf.putAtt(ncid,varid,'geoespatial_lat_max',max(lat));
netcdf.putAtt(ncid,varid,'geoespatial_lat_units','degree_north');
netcdf.putAtt(ncid,varid,'geoespatial_lon_min',min(lon));
netcdf.putAtt(ncid,varid,'geoespatial_lon_max',max(lon));
netcdf.putAtt(ncid,varid,'geoespatial_lon_units','degree_east');
netcdf.putAtt(ncid,varid,'time_coverage_start',datestr(min(time)));
netcdf.putAtt(ncid,varid,'time_coverage_end',datestr(max(time)));
netcdf.putAtt(ncid,varid,'time_coverage_resolution','monthly');
netcdf.putAtt(ncid,varid,'instrument','drifters');
netcdf.putAtt(ncid,varid,'cdm_data_type','grid');
netcdf.putAtt(ncid,varid,'data_type','grid data');
netcdf.putAtt(ncid,varid,'Conventions','CF-1.6, ACDD-1.3');
netcdf.putAtt(ncid,varid,'netcdf_version','4.3');
netcdf.putAtt(ncid,varid,'standard_name_vocabulary','CF Standard Name Table');
netcdf.putAtt(ncid,varid,'publisher_name','Sistema de Manejo Integral de Datos CIGOM');
netcdf.putAtt(ncid,varid,'publisher_email','smid@cigom.org');
netcdf.putAtt(ncid,varid,'publisher_url','http://smid.cigom.org');
netcdf.putAtt(ncid,varid,'acknowledgment','La base de datos GulfFlow (DOI:10.5281/zenodo.3978794) es producto del Consorcio de Investigación del Golfo de México (CIGoM) y fue parcialmente financiado por el Fondo Sectorial CONACYT-SENER-Hidrocarburos, México, proyecto 201441. The GulfFlow database (DOI:10.5281/zenodo.3978794) is the product of the Gulf of Mexico Research Consortium (CIGoM) and was partially funded by the CONACYT-SENER-Hydrocarbons Sector Fund, Mexico, project 201441.');
netcdf.putAtt(ncid,varid,'keywords','drifters, velocity, gulfflow, gridded fields, Gulf of Mexico, surface circulation');
netcdf.putAtt(ncid,varid,'date_created','05-May-2020');
netcdf.putAtt(ncid,varid,'date_modified',datestr(now,1));
netcdf.putAtt(ncid,varid,'history',[datestr(now,1) ' v. 1.1.0 Updated to span an additional year and to include additional Global Drifter Program data in the last two years; changed time to be in days; corrected fill value to be NaN | ' ...
                                    '05-May-2020 v. 1.0.0 initial version']);
netcdf.putAtt(ncid,varid,'processing_level','processed quality controlled by J. M. Lilly');
netcdf.putAtt(ncid,varid,'creator_type','group');
netcdf.putAtt(ncid,varid,'contributor_name','J. M. Lilly');
netcdf.putAtt(ncid,varid,'contributor_email','j.m.lilly@theissresearch.org');


% Define the dimensions
dimitime = netcdf.defDim(ncid,'time',length(time));
dimilon = netcdf.defDim(ncid,'lon',length(lon));
dimilat = netcdf.defDim(ncid,'lat',length(lat));
%dimiN = netcdf.defDim(ncid,'N',size(count,4));%jml removed this
dimiexp = netcdf.defDim(ncid,'sources',size(ExpName,1)); 
dimistr = netcdf.defDim(ncid,'string_length',size(ExpName,2)); 

%jml removed this
%N_ID = netcdf.defVar(ncid,'N','int',[dimiN]);
%netcdf.putAtt(ncid,N_ID,'standard_name','drogue_id');
%netcdf.putAtt(ncid,N_ID,'long_name','Identifier for drogued and undrogued drifters for each dataset, used for count and total_count');
%netcdf.putAtt(ncid,N_ID,'units','No units');

exp_str_ID = netcdf.defVar(ncid,'source_name','char',[dimistr dimiexp]);
netcdf.putAtt(ncid,exp_str_ID,'standard_name','dataset_ID');
netcdf.putAtt(ncid,exp_str_ID,'long_name','names and drogue statuses of drifter datasets');%jml edited this

%jml
time_ID = netcdf.defVar(ncid,'time','double',[dimitime]);
netcdf.putAtt(ncid,time_ID,'standard_name','serial_date_num');
netcdf.putAtt(ncid,time_ID,'long_name','serial date number, or Matlab''s datenum');
netcdf.putAtt(ncid,time_ID,'units','days since 00-Jan-0000 00:00:00');
netcdf.putAtt(ncid,time_ID,'axis','T')
netcdf.putAtt(ncid,time_ID,'calendar','standard');

%num_ID = netcdf.defVar(ncid,'num','double',[dimitime]);
%netcdf.putAtt(ncid,num_ID,'standard_name','serial_date_num');
%netcdf.putAtt(ncid,num_ID,'long_name','serial date number, Matlab datenum');
%netcdf.putAtt(ncid,num_ID,'units','days');
%netcdf.putAtt(ncid,num_ID,'axis','T')
%netcdf.putAtt(ncid,num_ID,'calendar','standard');

%time_ID = netcdf.defVar(ncid,'time','double',[dimitime]);
%netcdf.putAtt(ncid,time_ID,'standard_name','time');
%netcdf.putAtt(ncid,time_ID,'long_name','time of observations');
%netcdf.putAtt(ncid,time_ID,'units','seconds since 1900-01-01 00:00:00');
%netcdf.putAtt(ncid,time_ID,'axis','T')
%netcdf.putAtt(ncid,time_ID,'calendar','standard');

lat_ID = netcdf.defVar(ncid,'lat','double',[dimilat]);
netcdf.putAtt(ncid,lat_ID,'standard_name','latitude');
netcdf.putAtt(ncid,lat_ID,'long_name','central latitude of the grid');
netcdf.putAtt(ncid,lat_ID,'units','degrees_north');
netcdf.putAtt(ncid,lat_ID,'axis','Y');

lon_ID = netcdf.defVar(ncid,'lon','double',[dimilon]);
netcdf.putAtt(ncid,lon_ID,'standard_name','longitude');
netcdf.putAtt(ncid,lon_ID,'long_name','central longitude of the grid');
netcdf.putAtt(ncid,lon_ID,'units','degrees_east');
netcdf.putAtt(ncid,lon_ID,'axis','X');

u_ID = netcdf.defVar(ncid,'u','float',[dimilat dimilon dimitime]);
netcdf.putAtt(ncid,u_ID,'standard_name','eastward_sea_water_velocity');
netcdf.putAtt(ncid,u_ID,'long_name','monthly mean surface current east component');
netcdf.putAtt(ncid,u_ID,'units','m/s');
netcdf.defVarFill(ncid,u_ID,false,nan);
netcdf.defVarDeflate(ncid,u_ID,true,true,9);

v_ID = netcdf.defVar(ncid,'v','float',[dimilat dimilon dimitime]);
netcdf.putAtt(ncid,v_ID,'standard_name','northward_sea_water_velocity');
netcdf.putAtt(ncid,v_ID,'long_name','monthly mean surface current north component');
netcdf.putAtt(ncid,v_ID,'units','m/s');
netcdf.defVarFill(ncid,v_ID,false,nan)
netcdf.defVarDeflate(ncid,v_ID,true,true,9);

um_ID = netcdf.defVar(ncid,'um','float',[dimilat dimilon]);
netcdf.putAtt(ncid,um_ID,'standard_name','mean_eastward_sea_water_velocity');
netcdf.putAtt(ncid,um_ID,'long_name','global mean surface current east component');
netcdf.putAtt(ncid,um_ID,'units','m/s');
netcdf.defVarFill(ncid,um_ID,false,nan)
netcdf.defVarDeflate(ncid,um_ID,true,true,9);

vm_ID = netcdf.defVar(ncid,'vm','float',[dimilat dimilon]);
netcdf.putAtt(ncid,vm_ID,'standard_name','mean_northward_sea_water_velocity');
netcdf.putAtt(ncid,vm_ID,'long_name','global mean surface current north component');
netcdf.putAtt(ncid,vm_ID,'units','m/s');
netcdf.defVarFill(ncid,vm_ID,false,nan)
netcdf.defVarDeflate(ncid,vm_ID,true,true,9);

epsuu_ID = netcdf.defVar(ncid,'epsuu','float',[dimilat dimilon dimitime]);
netcdf.putAtt(ncid,epsuu_ID,'standard_name','u_autocovariance');
netcdf.putAtt(ncid,epsuu_ID,'long_name','local autocovariance of eastward velocity component');
netcdf.putAtt(ncid,epsuu_ID,'units','m^2/s^2');
netcdf.defVarFill(ncid,epsuu_ID,false,nan)
netcdf.defVarDeflate(ncid,epsuu_ID,true,true,9);

epsvv_ID = netcdf.defVar(ncid,'epsvv','float',[dimilat dimilon dimitime]);
netcdf.putAtt(ncid,epsvv_ID,'standard_name','v_autocovariance');
netcdf.putAtt(ncid,epsvv_ID,'long_name','local autocovariance of northward velocity component');
netcdf.putAtt(ncid,epsvv_ID,'units','m^2/s^2');
netcdf.defVarFill(ncid,epsvv_ID,false,nan)
netcdf.defVarDeflate(ncid,epsvv_ID,true,true,9);

epsuv_ID = netcdf.defVar(ncid,'epsuv','float',[dimilat dimilon dimitime]);
netcdf.putAtt(ncid,epsuv_ID,'standard_name','u_v_covariance');
netcdf.putAtt(ncid,epsuv_ID,'long_name','local covariance of northward and eastward velocity components');
netcdf.putAtt(ncid,epsuv_ID,'units','m^2/s^2');
netcdf.defVarFill(ncid,epsuv_ID,false,nan)
netcdf.defVarDeflate(ncid,epsuv_ID,true,true,9);

eps2_ID = netcdf.defVar(ncid,'eps2','float',[dimilat dimilon]);
netcdf.putAtt(ncid,eps2_ID,'standard_name','subgridscale_variability');
netcdf.putAtt(ncid,eps2_ID,'long_name','trace of the time-mean of the local covariance, representing variability unresolved by the gridded field');
netcdf.putAtt(ncid,eps2_ID,'units','m^2/s^2');
netcdf.defVarFill(ncid,eps2_ID,false,nan)
netcdf.defVarDeflate(ncid,eps2_ID,true,true,9);

count_ID = netcdf.defVar(ncid,'count','int',[dimilat dimilon dimitime dimiexp]);
netcdf.putAtt(ncid,count_ID,'standard_name','number_observations_dataset');
netcdf.putAtt(ncid,count_ID,'long_name','number of hourly observations of each source dataset (see ExpName). Values 115 are the count of velocities from drifters from each of the 15 experiments that have been flagged as drogued, values 15 through 30 for drifters that have been flagged as having lost their drogue, and values 31 through 45 for drifters whose drogue status is unknown.');
netcdf.defVarDeflate(ncid,count_ID,true,true,9);

total_count_ID = netcdf.defVar(ncid,'total_count','int',[dimilat dimilon dimitime]);
netcdf.putAtt(ncid,total_count_ID,'standard_name','number_observations');
netcdf.putAtt(ncid,total_count_ID,'long_name','number of hourly observations from all source datasets')%jml
netcdf.defVarDeflate(ncid,total_count_ID,true,true,9);
netcdf.endDef(ncid);
% Then store my main variable

%netcdf.putVar(ncid,N_ID,N);
netcdf.putVar(ncid,exp_str_ID,ExpName);
%netcdf.putVar(ncid,num_ID,num);
netcdf.putVar(ncid,time_ID,time);
netcdf.putVar(ncid,lat_ID,lat);
netcdf.putVar(ncid,lon_ID,lon);
netcdf.putVar(ncid,u_ID,u);
netcdf.putVar(ncid,v_ID,v);
netcdf.putVar(ncid,um_ID,um);
netcdf.putVar(ncid,vm_ID,vm);
netcdf.putVar(ncid,epsuu_ID,epsuu);
netcdf.putVar(ncid,epsvv_ID,epsvv);
netcdf.putVar(ncid,epsuv_ID,epsuv);
netcdf.putVar(ncid,eps2_ID,eps2);
netcdf.putVar(ncid,count_ID,count);
netcdf.putVar(ncid,total_count_ID,total_count);


netcdf.close(ncid)


function[]=gulfflow_onetwelfth_netcdf
%by Paula Garcia Carrillo with modifications by J.M. Lilly 
writedir=['/Users/lilly/Desktop/Dropbox/NetCDF/'];

load gulfflow_onetwelfth
time=gulfflow_onetwelfth.num; %jml
%timei=gulfflow_onetwelfth.num; num=timei;
%time=(timei*86400)-(datenum(1900,1,1,0,0,0)*86400);%%seconds
lon=gulfflow_onetwelfth.lon; lat=gulfflow_onetwelfth.lat;
u=gulfflow_onetwelfth.u;v=gulfflow_onetwelfth.v;
um=vmean(gulfflow_onetwelfth.u,3); vm=vmean(gulfflow_onetwelfth.v,3);%jml
epsuu=gulfflow_onetwelfth.epsuu;epsvv=gulfflow_onetwelfth.epsvv;epsuv=gulfflow_onetwelfth.epsuv;
eps2=vmean(gulfflow_onetwelfth.epsuu+gulfflow_onetwelfth.epsvv,3);%jml
count=gulfflow_onetwelfth.counts;total_count=gulfflow_onetwelfth.count;%jml
ExpName=gulfflow_onetwelfth.ExpName;%ExpName=char(ExpName);ExpName=ExpName';

%name_nc=['GulfFlow_GriddedFields_onetwelfth.nc'];
name_nc=['GulfFlow_OneTwelfth.nc'];
% netCDF
% Define the global attributes
ncid = netcdf.create([writedir name_nc],'NETCDF4');
varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid,varid,'site_code','Gulf of Mexico');
netcdf.putAtt(ncid,varid,'platform_code','drifters_gridded_velocity_fields');
netcdf.putAtt(ncid,varid,'title','GulfFlow gridded velocity product, 1/12 deg');
netcdf.putAtt(ncid,varid,'product_version','1.1.0');
netcdf.putAtt(ncid,varid,'summary','gridded surface velocity product for the Gulf of Mexico derived from surface drifters, 1/12 deg resolution');
netcdf.putAtt(ncid,varid,'project','Implementacion de redes de observacion oceanograficas (fisicas geoquimicas, ecologicas) para la generacion de escenarios ante posibles contingencias relacionadas a la exploracion y produccion de hidrocarburos en aguas profundas del Golfo de Mexico Consorcio de Investigacion del Golfo de Mexico (CIGoM)');
netcdf.putAtt(ncid,varid,'naming_authority','cigom.org');
netcdf.putAtt(ncid,varid,'id','gridded_surface_velocity_fields');
netcdf.putAtt(ncid,varid,'source','drifters');
netcdf.putAtt(ncid,varid,'principal_investigator','Paula Perez-Brunius');
netcdf.putAtt(ncid,varid,'principal_investigator_email','brunius@cicese.mx');
netcdf.putAtt(ncid,varid,'principal_investigator_url','https://giola.cicese.mx');
netcdf.putAtt(ncid,varid,'institution','Centro de Investigación Científica y de Educación Superior de Ensenada, CICESE');
netcdf.putAtt(ncid,varid,'area','Gulf of Mexico');
netcdf.putAtt(ncid,varid,'geoespatial_lat_min',min(lat));
netcdf.putAtt(ncid,varid,'geoespatial_lat_max',max(lat));
netcdf.putAtt(ncid,varid,'geoespatial_lat_units','degree_north');
netcdf.putAtt(ncid,varid,'geoespatial_lon_min',min(lon));
netcdf.putAtt(ncid,varid,'geoespatial_lon_max',max(lon));
netcdf.putAtt(ncid,varid,'geoespatial_lon_units','degree_east');
netcdf.putAtt(ncid,varid,'time_coverage_start',datestr(min(time)));
netcdf.putAtt(ncid,varid,'time_coverage_end',datestr(max(time)));
netcdf.putAtt(ncid,varid,'time_coverage_resolution','monthly');
netcdf.putAtt(ncid,varid,'instrument','drifters');
netcdf.putAtt(ncid,varid,'cdm_data_type','grid');
netcdf.putAtt(ncid,varid,'data_type','grid data');
netcdf.putAtt(ncid,varid,'Conventions','CF-1.6, ACDD-1.3');
netcdf.putAtt(ncid,varid,'netcdf_version','4.3');
netcdf.putAtt(ncid,varid,'standard_name_vocabulary','CF Standard Name Table');
netcdf.putAtt(ncid,varid,'publisher_name','Sistema de Manejo Integral de Datos CIGOM');
netcdf.putAtt(ncid,varid,'publisher_email','smid@cigom.org');
netcdf.putAtt(ncid,varid,'publisher_url','http://smid.cigom.org');
netcdf.putAtt(ncid,varid,'acknowledgment','La base de datos GulfFlow (DOI:10.5281/zenodo.3978794) es producto del Consorcio de Investigación del Golfo de México (CIGoM) y fue parcialmente financiado por el Fondo Sectorial CONACYT-SENER-Hidrocarburos, México, proyecto 201441. The GulfFlow database (DOI:10.5281/zenodo.3978794) is the product of the Gulf of Mexico Research Consortium (CIGoM) and was partially funded by the CONACYT-SENER-Hydrocarbons Sector Fund, Mexico, project 201441.');
netcdf.putAtt(ncid,varid,'keywords','drifters, velocity, gulfflow, gridded fields, Gulf of Mexico, surface circulation');
netcdf.putAtt(ncid,varid,'date_created','05-May-2020');
netcdf.putAtt(ncid,varid,'date_modified',datestr(now,1));
netcdf.putAtt(ncid,varid,'history',[datestr(now,1) ' v. 1.1.0 Updated to span an additional year and to include additional Global Drifter Program data in the last two years; changed time to be in days  | ' ...
                                    '05-May-2020 v. 1.0.0 initial version']);
netcdf.putAtt(ncid,varid,'processing_level','processed quality controlled by J. M. Lilly');
netcdf.putAtt(ncid,varid,'creator_type','group');
netcdf.putAtt(ncid,varid,'contributor_name','J. M. Lilly');
netcdf.putAtt(ncid,varid,'contributor_email','j.m.lilly@theissresearch.org');

% Define the dimensions
dimitime = netcdf.defDim(ncid,'time',length(time));
dimilon = netcdf.defDim(ncid,'lon',length(lon));
dimilat = netcdf.defDim(ncid,'lat',length(lat));
%dimiN = netcdf.defDim(ncid,'N',size(count,4));%jml removed this
dimiexp = netcdf.defDim(ncid,'sources',size(ExpName,1)); 
dimistr = netcdf.defDim(ncid,'string_length',size(ExpName,2)); 

%jml removed this
%N_ID = netcdf.defVar(ncid,'N','int',[dimiN]);
%netcdf.putAtt(ncid,N_ID,'standard_name','drogue_id');
%netcdf.putAtt(ncid,N_ID,'long_name','Identifier for drogued and undrogued drifters for each dataset, used for count and total_count');
%netcdf.putAtt(ncid,N_ID,'units','No units');

exp_str_ID = netcdf.defVar(ncid,'source_name','char',[dimistr dimiexp]);
netcdf.putAtt(ncid,exp_str_ID,'standard_name','dataset_ID');
netcdf.putAtt(ncid,exp_str_ID,'long_name','names and drogue statuses of drifter datasets');%jml edited this

%jml
time_ID = netcdf.defVar(ncid,'time','double',[dimitime]);
netcdf.putAtt(ncid,time_ID,'standard_name','serial_date_num');
netcdf.putAtt(ncid,time_ID,'long_name','serial date number, or Matlab''s datenum');
netcdf.putAtt(ncid,time_ID,'units','days since 00-Jan-0000 00:00:00');
netcdf.putAtt(ncid,time_ID,'axis','T')
netcdf.putAtt(ncid,time_ID,'calendar','standard');

%num_ID = netcdf.defVar(ncid,'num','double',dimitime);
%netcdf.putAtt(ncid,num_ID,'standard_name','serial_date_num');
%netcdf.putAtt(ncid,num_ID,'long_name','serial date number, Matlab datenum');
%netcdf.putAtt(ncid,num_ID,'units','days');
%netcdf.putAtt(ncid,num_ID,'axis','T')
%netcdf.putAtt(ncid,num_ID,'calendar','standard');

%time_ID = netcdf.defVar(ncid,'time','double',dimitime);
%netcdf.putAtt(ncid,time_ID,'standard_name','time');
%netcdf.putAtt(ncid,time_ID,'long_name','time of observations');
%netcdf.putAtt(ncid,time_ID,'units','seconds since 1900-01-01 00:00:00');
%netcdf.putAtt(ncid,time_ID,'axis','T')
%netcdf.putAtt(ncid,time_ID,'calendar','standard');

lat_ID = netcdf.defVar(ncid,'lat','double',[dimilat]);
netcdf.putAtt(ncid,lat_ID,'standard_name','latitude');
netcdf.putAtt(ncid,lat_ID,'long_name','central latitude of the grid');
netcdf.putAtt(ncid,lat_ID,'units','degrees_north');
netcdf.putAtt(ncid,lat_ID,'axis','Y');

lon_ID = netcdf.defVar(ncid,'lon','double',[dimilon]);
netcdf.putAtt(ncid,lon_ID,'standard_name','longitude');
netcdf.putAtt(ncid,lon_ID,'long_name','central longitude of the grid');
netcdf.putAtt(ncid,lon_ID,'units','degrees_east');
netcdf.putAtt(ncid,lon_ID,'axis','X');

u_ID = netcdf.defVar(ncid,'u','float',[dimilat dimilon dimitime]);
netcdf.putAtt(ncid,u_ID,'standard_name','eastward_sea_water_velocity');
netcdf.putAtt(ncid,u_ID,'long_name','monthly mean surface current east component');
netcdf.putAtt(ncid,u_ID,'units','m/s');
netcdf.defVarFill(ncid,u_ID,false,nan);
netcdf.defVarDeflate(ncid,u_ID,true,true,9);

v_ID = netcdf.defVar(ncid,'v','float',[dimilat dimilon dimitime]);
netcdf.putAtt(ncid,v_ID,'standard_name','northward_sea_water_velocity');
netcdf.putAtt(ncid,v_ID,'long_name','monthly mean surface current north component');
netcdf.putAtt(ncid,v_ID,'units','m/s');
netcdf.defVarFill(ncid,v_ID,false,nan);
netcdf.defVarDeflate(ncid,v_ID,true,true,9);

um_ID = netcdf.defVar(ncid,'um','float',[dimilat dimilon]);
netcdf.putAtt(ncid,um_ID,'standard_name','mean_eastward_sea_water_velocity');
netcdf.putAtt(ncid,um_ID,'long_name','global mean surface current east component');
netcdf.putAtt(ncid,um_ID,'units','m/s');
netcdf.defVarFill(ncid,um_ID,false,nan);
netcdf.defVarDeflate(ncid,um_ID,true,true,9);

vm_ID = netcdf.defVar(ncid,'vm','float',[dimilat dimilon]);
netcdf.putAtt(ncid,vm_ID,'standard_name','mean_northward_sea_water_velocity');
netcdf.putAtt(ncid,vm_ID,'long_name','global mean surface current north component');
netcdf.putAtt(ncid,vm_ID,'units','m/s');
netcdf.defVarFill(ncid,vm_ID,false,nan);
netcdf.defVarDeflate(ncid,vm_ID,true,true,9);

epsuu_ID = netcdf.defVar(ncid,'epsuu','float',[dimilat dimilon dimitime]);
netcdf.putAtt(ncid,epsuu_ID,'standard_name','u_autocovariance');
netcdf.putAtt(ncid,epsuu_ID,'long_name','local autocovariance of eastward velocity component');
netcdf.putAtt(ncid,epsuu_ID,'units','m^2/s^2');
netcdf.defVarFill(ncid,epsuu_ID,false,nan);
netcdf.defVarDeflate(ncid,epsuu_ID,true,true,9);

epsvv_ID = netcdf.defVar(ncid,'epsvv','float',[dimilat dimilon dimitime]);
netcdf.putAtt(ncid,epsvv_ID,'standard_name','v_autocovariance');
netcdf.putAtt(ncid,epsvv_ID,'long_name','local autocovariance of northward velocity component');
netcdf.putAtt(ncid,epsvv_ID,'units','m^2/s^2');
netcdf.defVarFill(ncid,epsvv_ID,false,nan);
netcdf.defVarDeflate(ncid,epsvv_ID,true,true,9);

epsuv_ID = netcdf.defVar(ncid,'epsuv','float',[dimilat dimilon dimitime]);
netcdf.putAtt(ncid,epsuv_ID,'standard_name','u_v_covariance');
netcdf.putAtt(ncid,epsuv_ID,'long_name','local covariance of northward and eastward velocity components');
netcdf.putAtt(ncid,epsuv_ID,'units','m^2/s^2');
netcdf.defVarFill(ncid,epsuv_ID,false,nan);
netcdf.defVarDeflate(ncid,epsuv_ID,true,true,9);

eps2_ID = netcdf.defVar(ncid,'eps2','float',[dimilat dimilon]);
netcdf.putAtt(ncid,eps2_ID,'standard_name','subgridscale_variability');
netcdf.putAtt(ncid,eps2_ID,'long_name','trace of the time-mean of the local covariance, representing variability unresolved by the gridded field');
netcdf.putAtt(ncid,eps2_ID,'units','m^2/s^2');
netcdf.defVarFill(ncid,eps2_ID,false,nan);
netcdf.defVarDeflate(ncid,eps2_ID,true,true,9);

count_ID = netcdf.defVar(ncid,'count','int',[dimilat dimilon dimitime dimiexp]);
netcdf.putAtt(ncid,count_ID,'standard_name','number_observations_dataset');
netcdf.putAtt(ncid,count_ID,'long_name','number of hourly observations of each source dataset (see ExpName). Values 115 are the count of velocities from drifters from each of the 15 experiments that have not been flagged as having lost their drogues, while values 1630 are for observation from drifters that have been flagged as having lost their drogue');
netcdf.defVarDeflate(ncid,count_ID,true,true,9);
 
total_count_ID = netcdf.defVar(ncid,'total_count','int',[dimilat dimilon dimitime]);
netcdf.putAtt(ncid,total_count_ID,'standard_name','number_observations');
netcdf.putAtt(ncid,total_count_ID,'long_name','number of hourly observations from all source datasets')%jml
netcdf.defVarDeflate(ncid,total_count_ID,true,true,9);

netcdf.endDef(ncid);
% Then store my main variable

%netcdf.putVar(ncid,N_ID,N);
netcdf.putVar(ncid,exp_str_ID,ExpName);
%snetcdf.putVar(ncid,num_ID,num);
netcdf.putVar(ncid,time_ID,time);
netcdf.putVar(ncid,lat_ID,lat);
netcdf.putVar(ncid,lon_ID,lon);
netcdf.putVar(ncid,u_ID,u);
netcdf.putVar(ncid,v_ID,v);
netcdf.putVar(ncid,um_ID,um);
netcdf.putVar(ncid,vm_ID,vm);
netcdf.putVar(ncid,epsuu_ID,epsuu);
netcdf.putVar(ncid,epsvv_ID,epsvv);
netcdf.putVar(ncid,epsuv_ID,epsuv);
netcdf.putVar(ncid,eps2_ID,eps2);
netcdf.putVar(ncid,count_ID,count);
netcdf.putVar(ncid,total_count_ID,total_count);


netcdf.close(ncid)


function[]=gulfflow_models

%See make_gulfdrifters for creating these files
%/*************************************************************************
%Monthly velocity with 50% overlap
lono=[-99:1/4:-80.5];
lato=[18:1/4:31];
%--------------------------------------------------------------------------
%Same for Aviso 
ncload('GulfDriftersAllAugmented.nc','ids','time','lat','lon','filled','depth','source','cmems_u','cmems_v');
use GulfDriftersAllAugmented
gulfflow_cmems=griddrifters(time,lat,lon,cmems_u,cmems_v,filled,source,lato,lono,1992,7,672);
matsave gulfflow_cmems
%--------------------------------------------------------------------------
%Same for Hycom
ncload('GulfDriftersAllAugmented.nc','ids','time','lat','lon','filled','depth','source','hycom_u','hycom_v');
use GulfDriftersAllAugmented
gulfflow_hycom=griddrifters(time,lat,lon,hycom_u,hycom_v,filled,source,lato,lono,1992,7,672);
matsave gulfflow_hycom
%--------------------------------------------------------------------------
%Same for nemo
ncload('GulfDriftersAllAugmented.nc','ids','time','lat','lon','filled','depth','source','nemo_u','nemo_v');
use GulfDriftersAllAugmented
gulfflow_nemo=griddrifters(time,lat,lon,nemo_u,nemo_v,filled,source,lato,lono,1992,7,672);
matsave gulfflow_nemo
%--------------------------------------------------------------------------
%Same for roms
ncload('GulfDriftersAllAugmented.nc','ids','time','lat','lon','filled','depth','source','roms_u','roms_v');
use GulfDriftersAllAugmented
gulfflow_roms=griddrifters(time,lat,lon,roms_u,roms_v,filled,source,lato,lono,1992,7,672);
matsave gulfflow_roms
%\*************************************************************************




