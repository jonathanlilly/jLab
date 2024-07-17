function[varargout]=make_gomed(varargin)
%MAKE_GOMED   Create a drifter-derived eddy database for the Gulf of Mexico.
%
%   This function creates the Gulf of Mexico Eddy Dataset, or GOMED, which
%   contains estimated properties of coherent eddy-like features in surface
%   drifter data.  For details, see
%
%       Lilly, J. M. and P. Perez-Brunius (2021). Extracting statistically
%           significant eddy signals from large Lagrangian datasets using
%           wavelet ridge analysis, with application to the Gulf of Mexico.
%           Nonlinear Processes in Geophysics. https://doi.org/10.5194/npg-28-181-2021
%
%   The dataset is available at 
%
%        Lilly, J. M. and P. Perez-Brunius (2021).  The Gulf of Mexico Eddy 
%           Dataset (GOMED), a census of statistically significant eddy-
%           like events from all available surface drifter data (V. 1.1.0). 
%           [Data set]. Zenodo. http://doi.org/10.5281/zenodo.3978803.
%
%   This function contains the code for creating GOMED.NC.  It is
%   provided for completeness, although it will not be possible for you to 
%   recreate it because some of the drifter source files are proprietary. 
%   The code can be run on the publicly available portion of the data by 
%   swapping 'GulfDriftersAll' for 'GulfDriftersOpen' below.  
%
%   Thanks to Paula Garcia Carrillo of CICESE for help with the conversion 
%   to NetCDF.
%
%   Usage: make_gomed --create     
%   _________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2020--2021 J.M. Lilly --- type 'help jlab_license' for details    

if nargin==0
    help about_gomed
elseif nargin>0
    if strcmpi(varargin{1}, '--create')
        gomed_create,return
    end
end

function[]=gomed_create

%This is where your files are kept
basedir='/Users/lilly/Desktop/Dropbox/NetCDF/';
    
%first, creating the eddy dataset

%use whichever version is available to you
ncload([basedir 'GulfDriftersAll.nc']);
gulfdrifters=GulfDriftersAll;clear GulfDriftersAll
%ncload([basedir 'GulfDriftersDWDE.nc']);
%gulfdrifters=GulfDriftersDWDE;clear GulfDriftersDWDE
%ncload([basedir 'GulfDriftersOpen.nc']);
%gulfdrifters=GulfDriftersOpen;clear GulfDriftersOpen

%/*************************************************************************
%Here we create random times series with the same spectra
use gulfdrifters

cv=cellpair(u,v);
cv=cellmult(100,cv);%convert m/s to cm/s
clear f spp snn lambda psi
tic
parfor i=1:length(lat)
    i
    [psi{i,1},lambda{i,1}]=sleptap(size(cv{i},1),4);
    [f{i,1},spp{i,1},snn{i,1}]=mspec(cv{i},psi{i,1},lambda{i,1},'adaptive');
end
toc;  %10 minutes on Talaria
matsave gulfspectra f spp snn lambda psi

use gulfspectra
clear gulfnoisedrifters
parfor j=1:10  %making a dataset 10x the size of my original
    j
    [gulfnoisedrifters{j}.lat,gulfnoisedrifters{j}.lon, gulfnoisedrifters{j}.cv]...
            =noisedrifters(time,lat,lon,cv,spp,snn);
    gulfnoisedrifters{j}.time=time;
end
matsave gulfnoisedrifters
%about 10 minutes on Heimdall

%load gulfnoisedrifters
use gulfnoisedrifters{1}

clear f spp snn
parfor i=1:length(lat)
    i
    [f{i,1},spp{i,1},snn{i,1}]=mspec(cv{i},psi{i},lambda{i},'adaptive');
end
toc;  %10 minutes on Talaria
matsave gulfnoisespectra f spp snn
%/*************************************************************************

%/*************************************************************************
%running the eddy analysis
use gulfdrifters
tic;gulfeddies61=eddyridges(time,lat,lon,2,1/64,sqrt(6),1,0,'parallel');time1=toc  %25 minutes on Heimdall
matsave gulfeddies61
notetoself
%minimum length with be 2MP/pi = sqrt(6)*2/pi=1.56

% use gulfeddies61
% f=cellmult(24,corfreq(latres));
% cell2col(omega,xi,f,Ro);
% omtilde=sign(xi).*omega./f;
% plot([omtilde Ro/2]),ylim([-1 1]*2)
clear gulfeddies

load gulfnoisedrifters
gulfnoisedrifters=catstruct(gulfnoisedrifters);
use gulfnoisedrifters
tic;gulfnoiseeddies61=eddyridges(time,lat,lon,2,1/64,sqrt(6),1,0,'parallel');time2=toc
matsave gulfnoiseeddies61  %about 3.5 hours on Heimdall
%\*************************************************************************

%Computation of the eddy significance levels is at the bottom, after 
%writing to the netcdf files

ncload([basedir 'GulfDriftersAll.nc']);
load gulfeddies61
load gulfnoiseeddies61

for ii=1:2
    %/*************************************************************************
    if ii==1
        use gulfeddies61
    else
        use gulfnoiseeddies61
    end
    %chi=cr;
    time=num;

    row_length=cellength(time);
    ridge_id=[1:length(time)]';
    seg_id=ridge_id;
    drifter_id=ridge_id;
    clear ids
    for i=1:length(time)
        ids{i}=i+0*time{i};
        kri=mod(kr{i}(1)-1,length(GulfDriftersAll.id))+1;%modulo for noise
        %aresame(kri,kr{i}(1))  %true for gulfeddies
        seg_id(i)=GulfDriftersAll.id(kri);
        drifter_id(i)=GulfDriftersAll.drifter_id(kri);
    end

    len=cellfirst(len);
    xi_bar=cellmean(xi,'weight',kappa);
    R_bar=sqrt(cellmean(cellmult(R,R)));
    V_bar=sign(xi_bar).*sqrt(cellmean(cellmult(V,V)));
    
    %find filled and drogued values, same for data and for noise
    filled=ids;
    drogue=ids;
    for i=1:length(ids)
        filled{i}=interp1(GulfDriftersAll.time{seg_id(i)},GulfDriftersAll.filled{seg_id(i)},time{i},'nearest');
        drogue{i}=interp1(GulfDriftersAll.time{seg_id(i)},GulfDriftersAll.drogue{seg_id(i)},time{i},'nearest');
    end
    %example:
    %%gomed.segment_id(131)=29
    %plot(GulfDriftersAll.time{29},GulfDriftersAll.filled{29}),hold on
    %cellplot(time(131),filled(131))

    cell2col(ids,time,lat,lon,latres,lonres,zhat,kappa,xi,theta,phi,omega,R,V,upsilon,chi,filled,drogue);


    f=corfreq(latres)./3600;
    omega_ast=sign(xi).*omega/24/3600./f;
    omega_ast_bar=cellmean(col2cell(omega_ast),'weight',col2cell(kappa));

    %remove nans marking tails
    bool=isnan(time);%length(find(bool))
    vindex(ids,time,lat,lon,latres,lonres,zhat,kappa,xi,theta,phi,omega,R,V,upsilon,chi,filled,drogue,f,omega_ast,~bool,1);

    %Following these convention
    %https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3
    %http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html#trajectory-data

    writedir=['/Users/lilly/Desktop/Dropbox/NetCDF/'];

    if ii==1
        name_nc=['gomed.nc'];
    elseif ii==2
        name_nc=['gomed_noise.nc'];
    end

    % netCDF
    % Define the global attributes
    ncid = netcdf.create([writedir name_nc],'NETCDF4');
    varid = netcdf.getConstant('GLOBAL');

    if ii==1
        netcdf.putAtt(ncid,varid,'title','GOMED');
        netcdf.putAtt(ncid,varid,'summary','GOMED, the Gulf of Mexico Eddy Dataset, is dataset containing identified eddy-like modulated oscillations and their properties within all available surface drifter trajectories from the Gulf of Mexico.')
    elseif ii==2
        netcdf.putAtt(ncid,varid,'title','GOMED_NOISE');
        netcdf.putAtt(ncid,varid,'summary','GOMED_NOISE is a version of the Gulf of Mexico Eddy Dataset that has been constructed from a "noise" dataset 10 times the size as the original drifter dataset, reflecting the null hypothesis of an isotropic background spectrum.')
    end

    netcdf.putAtt(ncid,varid,'product_version','1.1.0');
    netcdf.putAtt(ncid,varid,'keywords','Gulf of Mexico, surface drifters, Lagrangian data, coherent eddies, vortices, wavelet analysis');
    netcdf.putAtt(ncid,varid,'Conventions','CF-1.6, ACDD-1.3');
    netcdf.putAtt(ncid,varid,'references','Lilly and Perez-Prunius (2021b). Extracting statistically significant eddy signals from large Lagrangian datasets using wavelet ridge analysis, with application to the Gulf of Mexico.');

    if ii==1
        netcdf.putAtt(ncid,varid,'id','gomed');
    else
        netcdf.putAtt(ncid,varid,'id','gomed_noise');    
    end

    netcdf.putAtt(ncid,varid,'naming_authority','http://www.jmlilly.net');
    netcdf.putAtt(ncid,varid,'history',[datestr(now,1) ' v. 1.1.0 updated to include Global Drifter Program data through mid-2020, modified definition of statistical significance, added additional signifiance measures, and changed some variable names.' ... 
                        ' | 11-Aug-2020 v. 1.0.0 initial version']);
    netcdf.putAtt(ncid,varid,'source','surface drifter');
    netcdf.putAtt(ncid,varid,'platform','surface drifter');
    netcdf.putAtt(ncid,varid,'instrument','WOCE drifter, CODE drifter, SVP drifter, Far Horizon Drifter (FHD), CARTHE drifter, Microstar drifter, DORIS drifter');

    netcdf.putAtt(ncid,varid,'cdm_data_type','trajectory');
    netcdf.putAtt(ncid,varid,'featureType','trajectory');
    %netcdf.putAtt(ncid,varid,'data_type','trajectory data');

    netcdf.putAtt(ncid,varid,'processing_level','Variable levels of processing by original investigators, followed by uniform processing and quality control as described in Lilly and Perez-Brunius (2021a) and references therein.');
    netcdf.putAtt(ncid,varid,'comment',['Thanks to Paula Garcia Carrillo for help with the NetCDF conversion.  ' ...
       'To convert this dataset to cell array form in Matlab, with the jLab toolbox installed from https://github.com/jonathanlilly/jLab, use ncload(''gomed.nc'')']);

    netcdf.putAtt(ncid,varid,'acknowledgment','La base de datos GOMED es producto del Consorcio de Investigacion del Golfo de Mexico (CIGoM) y fue parcialmente financiado por el Fondo Sectorial CONACYT-SENER-Hidrocarburos, Mexico, proyecto 201441  (http://doi.org/10.5281/zenodo.3978804).  The GOMED database is the product of the Gulf of Mexico Research Consortium (CIGoM) and was partially funded by the CONACYT-SENER-Hydrocarbons Sector Fund, Mexico, project 201441 (http://doi.org/10.5281/zenodo.3978804).  The work of J. M. Lilly was partly supported by award #1658564 from the Physical Oceanography program of the United States National Science Foundation.')
    netcdf.putAtt(ncid,varid,'license','This dataset is licensed for non-commercial academic use only, and it is prohibited to share it with third parties, as well as to profit or sell products derived from it.');

    netcdf.putAtt(ncid,varid,'standard_name_vocabulary','CF Standard Name Table');
    netcdf.putAtt(ncid,varid,'date_created','11-Aug-2020');
    netcdf.putAtt(ncid,varid,'date_modified',datestr(now,1));
    netcdf.putAtt(ncid,varid,'creator_name','Jonathan Lilly');
    netcdf.putAtt(ncid,varid,'creator_email','eponym@jmlilly.net');
    netcdf.putAtt(ncid,varid,'creator_url','http://www.jmlilly.net');
    netcdf.putAtt(ncid,varid,'creator_institution','Theiss Research');
    netcdf.putAtt(ncid,varid,'institution','Centro de Investigacion Cientifica y de Educacion Superior de Ensenada, CICESE');
    netcdf.putAtt(ncid,varid,'contributor_name','Jonathan M. Lilly (80%), Paula Perez Brunius (20%)');
    netcdf.putAtt(ncid,varid,'contributor_email','j.m.lilly@theissresearch.org, brunius@cicese.mx');
    netcdf.putAtt(ncid,varid,'contributor_institution','Theiss Research, CICESE');
    netcdf.putAtt(ncid,varid,'project','Implementacion de redes de observacion oceanograficas (fisicas geoquimicas, ecologicas) para la generacion de escenarios ante posibles contingencias relacionadas a la exploracion y produccion de hidrocarburos en aguas profundas del Golfo de Mexico Consorcio de Investigacion del Golfo de Mexico (CIGoM)');
    netcdf.putAtt(ncid,varid,'publisher_name','Sistema de Manejo Integral de Datos CIGoM')
    netcdf.putAtt(ncid,varid,'publisher_email','smid@cigom.org');
    netcdf.putAtt(ncid,varid,'publisher_url','http://smid.cigom.org')
    netcdf.putAtt(ncid,varid,'principal_investigator','Paula Perez Brunius')
    netcdf.putAtt(ncid,varid,'principal_investigator_email','brunius@cicese.mx');
    netcdf.putAtt(ncid,varid,'principal_investigator_url','https://giola.cicese.mx')
    netcdf.putAtt(ncid,varid,'geoespatial_lat_min',min(lat));
    netcdf.putAtt(ncid,varid,'geoespatial_lat_max',max(lat));
    netcdf.putAtt(ncid,varid,'geoespatial_lat_units','degree_north');
    netcdf.putAtt(ncid,varid,'geoespatial_lon_min',min(lon));
    netcdf.putAtt(ncid,varid,'geoespatial_lon_max',max(lon));
    netcdf.putAtt(ncid,varid,'geoespatial_lon_units','degree_east');
    netcdf.putAtt(ncid,varid,'time_coverage_start',datestr(min(time)));
    netcdf.putAtt(ncid,varid,'time_coverage_end',datestr(max(time)));
    netcdf.putAtt(ncid,varid,'time_coverage_resolution','1h');

    %netcdf.putAtt(ncid,varid,'site_code','Gulf of Mexico');
    %netcdf.putAtt(ncid,varid,'platform_code','drifter');
    %netcdf.putAtt(ncid,varid,'area','Gulf of Mexico');

    % Define the dimensions
    dimiobs   = netcdf.defDim(ncid,'obs',length(time));
    dimifloat = netcdf.defDim(ncid,'traj',length(ridge_id));

    dimisig   = netcdf.defDim(ncid,'rho_cols',8);

    ridge_ID = netcdf.defVar(ncid,'ridge_id','int',dimifloat);
    netcdf.putAtt(ncid,ridge_ID,'long_name','ridge ID');
    netcdf.putAtt(ncid,ridge_ID,'cf_role','trajectory_id');

    seg_ID = netcdf.defVar(ncid,'segment_id','int',dimifloat);
    netcdf.putAtt(ncid,seg_ID,'long_name','ID of trajectory segment from within GulfDriftersAll');

    drifter_ID = netcdf.defVar(ncid,'drifter_id','int',dimifloat);
    netcdf.putAtt(ncid,drifter_ID,'long_name','drifter ID from within GulfDriftersAll');

    row_size_ID = netcdf.defVar(ncid,'row_size','int',dimifloat);
    netcdf.putAtt(ncid,row_size_ID,'long_name','number of observations for each trajectory segment');
    netcdf.putAtt(ncid,row_size_ID,'sample_dimension','obs');

    L_ID = netcdf.defVar(ncid,'L','double',dimifloat);
    netcdf.putAtt(ncid,L_ID,'long_name','ridge length in number of cycles completed');
    netcdf.putAtt(ncid,L_ID,'units','dimensionless');
    netcdf.putAtt(ncid,L_ID,'sample_dimension','obs');

    omega_ast_bar_ID = netcdf.defVar(ncid,'omega_ast_bar','double',dimifloat);
    netcdf.putAtt(ncid,omega_ast_bar_ID,'long_name','ridge-averaged nondimensional instantaneous frequency');
    netcdf.putAtt(ncid,omega_ast_bar_ID,'units','dimensionless');
    netcdf.putAtt(ncid,omega_ast_bar_ID,'sample_dimension','obs');

    xi_bar_ID = netcdf.defVar(ncid,'xi_bar','double',dimifloat);
    netcdf.putAtt(ncid,xi_bar_ID,'long_name','ridge-averaged circularity');
    netcdf.putAtt(ncid,xi_bar_ID,'units','dimensionless');
    netcdf.putAtt(ncid,xi_bar_ID,'sample_dimension','obs');

    R_bar_ID = netcdf.defVar(ncid,'R_bar','double',dimifloat);
    netcdf.putAtt(ncid,R_bar_ID,'long_name','ridge-averaged geometric mean radius');
    netcdf.putAtt(ncid,R_bar_ID,'units','km');
    netcdf.putAtt(ncid,R_bar_ID,'sample_dimension','obs');

    V_bar_ID = netcdf.defVar(ncid,'V_bar','double',dimifloat);
    netcdf.putAtt(ncid,V_bar_ID,'long_name','ridge-averaged kinetic energy velocity');
    netcdf.putAtt(ncid,V_bar_ID,'units','m/s');
    netcdf.putAtt(ncid,V_bar_ID,'sample_dimension','obs');

    rho_ID = netcdf.defVar(ncid,'rho','double',[dimifloat dimisig]);
    netcdf.putAtt(ncid,rho_ID,'long_name','ridge significance level ratio using L (first column), L*|xi_bar| (second column), L*|xi_bar|^2 (third column), ... , L*|xi_bar|^6 (seventh column), and |xi_bar| (eigth column)');
    netcdf.putAtt(ncid,rho_ID,'units','dimensionless');
    netcdf.putAtt(ncid,rho_ID,'sample_dimension','obs');

    % Define the main variables

    ids_ID = netcdf.defVar(ncid,'ids','int',dimiobs);
    netcdf.putAtt(ncid,ids_ID,'long_name','ridge IDs for all data points');
    netcdf.defVarDeflate(ncid,ids_ID,true,true,9);

    time_ID = netcdf.defVar(ncid,'time','double',dimiobs);
    netcdf.putAtt(ncid,time_ID,'standard_name','time');
    netcdf.putAtt(ncid,time_ID,'long_name','time');
    netcdf.putAtt(ncid,time_ID,'units','Days since 00-Jan-0000 00:00:00');
    netcdf.putAtt(ncid,time_ID,'axis','T')
    netcdf.putAtt(ncid,time_ID,'calendar','standard');
    netcdf.defVarDeflate(ncid,time_ID,true,true,9);

    filled_ID = netcdf.defVar(ncid,'filled','double',dimiobs);
    netcdf.putAtt(ncid,filled_ID,'long_name','filled on / off flag');
    netcdf.putAtt(ncid,filled_ID,'valid_range','0, 1');
    netcdf.putAtt(ncid,filled_ID,'flag_values','0, 1');
    netcdf.putAtt(ncid,filled_ID,'flag_meanings','at_least_one_valid_original_datapoint_within_plus_or_minus_three_hours otherwise');
    netcdf.defVarDeflate(ncid,filled_ID,true,true,9);

    drogue_ID = netcdf.defVar(ncid,'drogue','double',dimiobs);
    netcdf.putAtt(ncid,drogue_ID,'long_name','drogue on / off flag');
    netcdf.putAtt(ncid,drogue_ID,'valid_range','0, 1');
    netcdf.putAtt(ncid,drogue_ID,'flag_values','0, 1, NaN');
    netcdf.putAtt(ncid,drogue_ID,'flag_meanings','drogue_off drogue_on unknown');
    netcdf.defVarFill(ncid,drogue_ID,false,nan);%since there are undefined values
    netcdf.defVarDeflate(ncid,drogue_ID,true,true,9);

    Y_ID = netcdf.defVar(ncid,'lat','double',dimiobs);
    netcdf.putAtt(ncid,Y_ID,'standard_name','latitude');
    netcdf.putAtt(ncid,Y_ID,'long_name','latitude');
    netcdf.putAtt(ncid,Y_ID,'units','degrees_north');
    netcdf.putAtt(ncid,Y_ID,'axis','Y');
    netcdf.defVarDeflate(ncid,Y_ID,true,true,9);

    X_ID = netcdf.defVar(ncid,'lon','double',dimiobs);
    netcdf.putAtt(ncid,X_ID,'standard_name','longitude');
    netcdf.putAtt(ncid,X_ID,'long_name','longitude');
    netcdf.putAtt(ncid,X_ID,'units','degrees_east');
    netcdf.putAtt(ncid,X_ID,'axis','X');
    netcdf.defVarDeflate(ncid,X_ID,true,true,9);

    Yres_ID = netcdf.defVar(ncid,'latres','double',dimiobs);
    netcdf.putAtt(ncid,Yres_ID,'long_name','estimated time-varying central latitude of ridge');
    netcdf.putAtt(ncid,Yres_ID,'units','degrees_north');
    netcdf.putAtt(ncid,Yres_ID,'axis','Y');
    netcdf.defVarDeflate(ncid,Yres_ID,true,true,9);

    Xres_ID = netcdf.defVar(ncid,'lonres','double',dimiobs);
    netcdf.putAtt(ncid,Xres_ID,'long_name','estimated time-varying central longitude of ridge');
    netcdf.putAtt(ncid,Xres_ID,'units','degrees_east');
    netcdf.putAtt(ncid,Xres_ID,'axis','X');
    netcdf.defVarDeflate(ncid,Xres_ID,true,true,9);

    Xhat_ID = netcdf.defVar(ncid,'x_star','double',dimiobs);
    netcdf.putAtt(ncid,Xhat_ID,'long_name','eastward displacement associated with eddy ellipse');
    netcdf.putAtt(ncid,Xhat_ID,'units','km');
    netcdf.putAtt(ncid,Xhat_ID,'axis','X');
    netcdf.defVarDeflate(ncid,Xhat_ID,true,true,9);

    Yhat_ID = netcdf.defVar(ncid,'y_star','double',dimiobs);
    netcdf.putAtt(ncid,Yhat_ID,'long_name','northward displacement associated with eddy ellipse');
    netcdf.putAtt(ncid,Yhat_ID,'units','km');
    netcdf.putAtt(ncid,Yhat_ID,'axis','Y');
    netcdf.defVarDeflate(ncid,Yhat_ID,true,true,9);

    kappa_ID = netcdf.defVar(ncid,'kappa','double',dimiobs);
    netcdf.putAtt(ncid,kappa_ID,'long_name','root-mean-square radius of eddy ellipse');
    netcdf.putAtt(ncid,kappa_ID,'units','km');
    netcdf.defVarDeflate(ncid,kappa_ID,true,true,9);

    xi_ID = netcdf.defVar(ncid,'xi','double',dimiobs);
    netcdf.putAtt(ncid,xi_ID,'long_name','eddy ellipse circularity 2ab/(a^2+b^2), positive = counterclockwise motion, negative = clockwise motion');
    netcdf.putAtt(ncid,xi_ID,'units','dimensionless');
    netcdf.defVarDeflate(ncid,xi_ID,true,true,9);

    theta_ID = netcdf.defVar(ncid,'theta','double',dimiobs);
    netcdf.putAtt(ncid,theta_ID,'long_name','angle of eddy ellipse major axis counterclockwise from x-axis');
    netcdf.putAtt(ncid,theta_ID,'units','rad');
    netcdf.defVarDeflate(ncid,theta_ID,true,true,9);

    phi_ID = netcdf.defVar(ncid,'phi','double',dimiobs);
    netcdf.putAtt(ncid,phi_ID,'long_name','phase angle setting location of particle with respect to eddy ellipse major axis');
    netcdf.putAtt(ncid,phi_ID,'units','rad');
    netcdf.defVarDeflate(ncid,phi_ID,true,true,9);

    omega_ID = netcdf.defVar(ncid,'omega','double',dimiobs);
    netcdf.putAtt(ncid,omega_ID,'long_name','instantaneous frequency of eddy ellipse');
    netcdf.putAtt(ncid,omega_ID,'units','rad/s');
    netcdf.defVarDeflate(ncid,omega_ID,true,true,9);

    omega_ast_ID = netcdf.defVar(ncid,'omega_ast','double',dimiobs);
    netcdf.putAtt(ncid,omega_ast_ID,'long_name','nondimensional instantaneous frequency of eddy ellipse, positive = cyclonic motion, negative = anticyclonic motion');
    netcdf.putAtt(ncid,omega_ast_ID,'units','dimensionless');
    netcdf.defVarDeflate(ncid,omega_ast_ID,true,true,9);

    R_ID = netcdf.defVar(ncid,'R','double',dimiobs);
    netcdf.putAtt(ncid,R_ID,'long_name','geometric mean radius of eddy ellipse');
    netcdf.putAtt(ncid,R_ID,'units','km');
    netcdf.defVarDeflate(ncid,R_ID,true,true,9);

    V_ID = netcdf.defVar(ncid,'V','double',dimiobs);
    netcdf.putAtt(ncid,V_ID,'long_name','kinetic energy velocity of eddy ellipse, positive = counterclockwise motion, negative = clockwise motion');
    netcdf.putAtt(ncid,V_ID,'units','m/s');
    netcdf.defVarDeflate(ncid,V_ID,true,true,9);

    chi_ID = netcdf.defVar(ncid,'chi','double',dimiobs);
    netcdf.putAtt(ncid,chi_ID,'long_name','estimated nondimensional bias error magnitude of eddy ellipse');
    netcdf.putAtt(ncid,chi_ID,'units','dimensionless');
    netcdf.defVarDeflate(ncid,chi_ID,true,true,9);

    % We are done defining the NetCdf
    netcdf.endDef(ncid);
    % Then store the dimension variables in
    netcdf.putVar(ncid,ridge_ID,ridge_id);
    netcdf.putVar(ncid,seg_ID,seg_id);
    netcdf.putVar(ncid,drifter_ID,drifter_id);
    netcdf.putVar(ncid,row_size_ID,row_length);
    netcdf.putVar(ncid,L_ID,len);
    netcdf.putVar(ncid,omega_ast_bar_ID,omega_ast_bar);
    netcdf.putVar(ncid,xi_bar_ID,xi_bar);
    netcdf.putVar(ncid,R_bar_ID,R_bar);
    netcdf.putVar(ncid,V_bar_ID,V_bar/100);%convert cm/s to m/s
    netcdf.putVar(ncid,rho_ID,nan*zeros(length(xi_bar),8));
    % Then store my main variable
    netcdf.putVar(ncid,ids_ID,ids);
    netcdf.putVar(ncid,time_ID,time);
    netcdf.putVar(ncid,drogue_ID,drogue);
    netcdf.putVar(ncid,filled_ID,filled);
    netcdf.putVar(ncid,Y_ID,lat);
    netcdf.putVar(ncid,X_ID,lon);
    netcdf.putVar(ncid,Yres_ID,latres);
    netcdf.putVar(ncid,Xres_ID,lonres);
    netcdf.putVar(ncid,Xhat_ID,real(zhat)); 
    netcdf.putVar(ncid,Yhat_ID,imag(zhat));
    netcdf.putVar(ncid,kappa_ID,kappa); 
    netcdf.putVar(ncid,xi_ID,xi);
    netcdf.putVar(ncid,theta_ID,theta);
    netcdf.putVar(ncid,phi_ID,phi);
    netcdf.putVar(ncid,omega_ID,omega/24/3600);%convert rad/day to rad/s
    netcdf.putVar(ncid,omega_ast_ID,omega_ast);
    netcdf.putVar(ncid,R_ID,R);
    netcdf.putVar(ncid,V_ID,V/100);%convert cm/s to m/s
    netcdf.putVar(ncid,chi_ID,chi);

    % We're done, close the netCDF
    netcdf.close(ncid)
    %\*************************************************************************
end

%/*************************************************************************
%eddy significance levels
clear  gulfeddies61  gulfnoiseeddies61

clear gomed gomed_noise
ncload([writedir 'gomed.nc']);
ncload([writedir 'gomed_noise.nc']);

dx=1/4;
xmid=[-2:0.01:2];
N=10;

x=gomed.omega_ast_bar;
xnoise=gomed_noise.omega_ast_bar;

rho=vzeros(length(x),8);
rho_noise=vzeros(length(xnoise),8);
for i=1:8
    i
    if i==8
        y=abs(gomed.xi_bar);
        ynoise=abs(gomed_noise.xi_bar);
        ybins=[0:0.001:1];
    else
        y=(gomed.L).*abs(gomed.xi_bar).^(i-1);
        ynoise=(gomed_noise.L).*abs(gomed_noise.xi_bar).^(i-1);
        ybins=[0:0.005:20];
    end
    [rho(:,i),ymid,~,~,Srat]=eddylevels(dx,xmid,ybins,x,y,xnoise,ynoise,N,'symmetric','sort');
    rho_noise(:,i)=10.^interp2(xmid,ymid,log10(Srat),xnoise,ynoise);
end

name_nc=['gomed.nc'];
ncid = netcdf.open([writedir name_nc],'write');
netcdf.putVar(ncid,netcdf.inqVarID(ncid,'rho'),rho);
netcdf.close(ncid)

name_nc=['gomed_noise.nc'];
ncid = netcdf.open([writedir name_nc],'write');
netcdf.putVar(ncid,netcdf.inqVarID(ncid,'rho'),rho_noise);
netcdf.close(ncid)
%\*************************************************************************
