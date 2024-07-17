function[varargout]=make_gulfdrifters(varargin)
%MAKE_GULFDRIFTERS   Create a drifter dataset for the Gulf of Mexico.
%
%   GULFDRIFTERS creates a dataset containing all publicly available 
%   surface drifter data from the Gulf of Mexico, as described in 
%
%        Lilly, J. M. and P. Perez-Brunius (2021a). A gridded surface
%            current product for the Gulf of Mexico from consolidated  
%            drifter measurements.  Submitted to Earth Science System Data.
%
%   The dataset is available at 
%
%        Lilly, J. M. and P. Perez-Brunius (2021). GulfDrifters: A 
%            consolidated surface drifter dataset for the Gulf of Mexico 
%            [Data set]. Zenodo. https://doi.org/10.5281/zenodo.3985916
%
%   This function contains the code for creating GULFDRIFTERS.NC.  It is
%   provided for completeness, although it will not be possible for you to 
%   run in its entirety because some of the source files are proprietary.
%
%   Thanks to Paula Garcia Carrillo of CICESE for help with the conversion 
%   to NetCDF.
%
%   Usage: make_gulfdrifters --create     
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2020--2021 J.M. Lilly --- type 'help jlab_license' for details    

if nargin==0
    help about_gulfdrifters
elseif nargin>0
    if strcmpi(varargin{1}, '--create')
        gulfdrifters_create
    elseif strcmpi(varargin{1}, '--individual')
        gulfdrifters_individual
    elseif strcmpi(varargin{1}, '--primary')
        gulfdrifters_primary
    elseif strcmpi(varargin{1}, '--secondary')
        gulfdrifters_secondary
    elseif strcmpi(varargin{1}, '--tertiary')
        gulfdrifters_tertiary
    elseif strcmpi(varargin{1}, '--netcdf')
        gulfdrifters_netcdf
    elseif strcmpi(varargin{1}, '--augmented')
        gulfdrifters_augmented
    end
end


function[]=gulfdrifters_create

%You'll need to edit these
aomldir='/Volumes/Myrkheim/aoml/hourly';
sourcedir='/Users/lilly/Desktop/Dropbox/Data/gulfdata/source';
netcdfdir='/Users/lilly/Desktop/Dropbox/NetCDF';

disp('Creating GULFDRIFTERS.MAT.  This may take while...')
gulfdrifters_individual(sourcedir);
gulfdrifters_primary;
gulfdrifters_secondary;
gulfdrifters_tertiary;
gulfdrifters_netcdf(netcdfdir);

function[]=gulfdrifters_individual(sourcedir,aomldir)

disp('Processing invididual datasets ... ')

%region1=[-99 -89 17 30.5]; 
%region2=[-89 -81.5 22 30.5];

region1=[-99 -88.5 17 30.5]; %western region
region2=[-88.5 -84.5 21.5 30.5]; %central region
region3=[-84.5 -81.5 22.3 30.5]; %eastern region
region4=[-81.5 -80.5 22.3 26]; %channel region

%regionplot(region1),hold on,regionplot(region2),regionplot(region3),regionplot(region4),


%/*************************************************************************
% Hourly resolution GPS surface drifters from the Global Drifter Program
disp('Processing hourly GPS dataset... ')

ncload('HourlyDrifters','ids','time','lat','lon','u','v','drogue','gap','filled','gps')
%this is J. Lilly's version of the global drifter program's hourly dataset
use HourlyDrifters
num=time;
id=ids;
vindex(lat,lon,id,num,u,v,gap,drogue,filled,gps==1,1);

[lat,lon,id,num,u,v,gap,drogue,filled]=...
    trajextract(lat,lon,id,num,u,v,gap,drogue,filled,{region1,region2,region3,region4});
%vsize(lat,lon,id,num,u,v,gap,drogue,filled)

cv=cellpair(u,v);
t=num;
for i=1:length(num)
    t{i}=nan*t{i};
end
matsave gulfdrifters_hgps id num lat lon cv t drogue filled gap
clear id num lat lon cv t drogue filled time ids u v 
%\*************************************************************************

% /*************************************************************************
% Hourly resolution Argos surface drifters from the Global Drifter Program
disp('Processing hourly Argos GDP dataset... ')

use HourlyDrifters 
num=time;
id=ids;
vindex(lat,lon,id,num,u,v,gap,drogue,filled,gps==0,1);

[lat,lon,id,num,u,v,gap,drogue,filled]=...
    trajextract(lat,lon,id,num,u,v,gap,drogue,filled,{region1,region2,region3,region4});

cv=cellpair(u,v);
t=num;
for i=1:length(num)
    t{i}=nan*t{i};
end

%don't change this regarrangement!  It's so the cells line up with those
%of an earlier region definition so I don't have to redo my visual QC
vindex(id,num,lat,lon,cv,t,gap,drogue,filled,[1:13 18:length(id) 14:17],1);
matsave gulfdrifters_hargos id num lat lon cv t drogue filled gap
clear id num lat lon cv t drogue filled time ids u v 
%\*************************************************************************

%/*************************************************************************
%Global Drifter Program
disp('Processing GDP dataset... ')

ncload GlobalDrifters
%this is J. Lilly's version of the global drifter program's six-hourly dataset
use GlobalDrifters

t=temp;
num=time;
id=ids;
cv=cellpair(u,v);

[lat,lon,id,num,cv,t,drogue,filled]=...
    trajextract(lat,lon,id,num,cv,t,drogue,filled,{region1,region2,region3,region4});
make gulfdrifters_gdp id num lat lon cv t drogue 
%remove those for which id was already accounted for in hargos and hgps
id=cellfirst(gulfdrifters_gdp.id);
id1=cellfirst(gulfdrifters_hgps.id);
id2=cellfirst(gulfdrifters_hargos.id);
bool=false(size(id));
for i=1:length(id)
    if anyany(id1==id(i))||anyany(id2==id(i))
       bool(i)=true;
    end
end
% %Verify that these are covered by gps and argos
% figure
% use gulfdrifters.gdp
% cellplot(lon(bool),lat(bool),'k'),hold on
% use gulfdrifters.hgps
% cellplot(lon,lat,'r')
% use gulfdrifters.hargos
% cellplot(lon,lat,'g')
% %good
use gulfdrifters_gdp 
vindex(id,num,lat,lon,cv,drogue,t,filled,~bool,1);
matsave gulfdrifters_gdp id num lat lon cv t drogue filled
% %Verify that these are not covered by gps and argos
% figure
% use gulfdrifters_gdp,hold on
% cellplot(lon,lat,'2k')
% use gulfdrifters_hgps,hold on
% cellplot(lon,lat,'0.5r')
% use gulfdrifters_hargos
% cellplot(lon,lat,'0.5g')
%\*************************************************************************

%/*************************************************************************
%Southern Gulf of Mexico drifters
disp('Processing SGOM dataset... ')

%Columns are buoy#, year, month, day, hr(GMT), min(GMT), latitude(N), 
%longitude(E), speed(m/s), direction(degrees from north), flag1, flag2, 
%flag3, flag4,flag5, which are qc flags. Data have been linearly 
%interpolated between gaps. Bad data are set to -9999 if transmitting from 
%land (flag 1), with speeds larger than 3m/s (Flag3), or if  interpolation 
%gaps are larger than 6hrs (Flag4). Flag5 indicates speed peaks larger than
%2 standard deviations of data 12 hrs before and after. Flag2 if drifter is
%at less than 50m bottom depth (drogue is at 45m).

%Working with processed data
dirname=[sourcedir '/sgom/qc_data'];
files=findfiles(dirname,'dat'); 
clear id num lat lon cv flag1 flag2 flag3 flag4 flag5
for i=1:length(files)
    data=load([dirname '/' files{i}]);   
    id{i,1}=data(:,1);
    num{i,1}=datenum(data(:,2),data(:,3),data(:,4),data(:,5),data(:,6),0*data(:,6));
    lat{i,1}=vswap(data(:,7),-9999,nan);  
    lon{i,1}=vswap(data(:,8),-9999,nan);  
    cv{i,1}=100*vswap(data(:,9),-9999,nan).*rot(-2*pi/360*data(:,10)+pi/2);
    flag1{i,1}=data(:,11);
    flag2{i,1}=data(:,12);
    flag3{i,1}=data(:,13);
    flag4{i,1}=data(:,14);
    flag5{i,1}=data(:,15);
end

%note that flag4 says when the interpolation gap is greater than 6 hours,
%but the data is NaN there anyway so we don't need it

%Something has happened with id 80426:  time is out of joint
jj=find(cellfirst(id)==80426);
[num{jj},index]=sort(num{jj});
[lat{jj},lon{jj},cv{jj},flag1{jj},flag2{jj},flag3{jj},flag4{jj},flag5{jj}]=vindex(lat{jj},lon{jj},cv{jj},...
    flag1{jj},flag2{jj},flag3{jj},flag4{jj},flag5{jj},index,1);
filled=flag4;
[lat,lon,id,num,cv,filled]=trajextract(lat,lon,id,num,cv,filled,{region1,region2,region3,region4});
matsave gulfdrifters_sgom id num lat lon cv filled
%\*************************************************************************

% %Remove all nans 
% cellpack(num,id,lat,lon,cv,flag1,flag2,flag3,flag4,flag5);
% 
% vsize(num,id,lat,lon,cv,flag1,flag2,flag3,flag4,flag5)
% cellsplit(num,id,lat,lon,cv,flag1,flag2,flag3,flag4,flag5,2);%Split with longest gap at 2 days
% vsize(num,id,lat,lon,cv,flag1,flag2,flag3,flag4,flag5)
% 
% %Check gap statistics
% [dt,sigdt,meddt,maxdt]=sampletimes(num);
% 
% %Interpolate to regular grid
% cellgrid(num,id,lat,lon,cv,flag1,flag2,flag3,flag4,flag5,1/24);
% 
% %Fill bad data points
% %cellfill(num,id,lat,lon,cv,flag1,flag2,flag3,flag4,flag5);
% 
% [dt,sigdt,meddt,maxdt]=sampletimes(num);
% aresame(max(maxdt),min(dt),1e-4)
% 
% %Number of interpolated data points
% length(nonnan(cell2col(num)))  %879177
% 
% %Limit to the region
% [lat,lon,id,num,cv,flag1,flag2,flag3,flag4,flag5]=trajextract(lat,lon,id,num,cv,flag1,flag2,flag3,flag4,flag5,{region1,region2,region3,region4});
% 
% %Number of interpolated data points
% length(nonnan(cell2col(num)))  %874769
% 
% %Remove segments shorter than one week
% index=find(cellength(num)>7*24);  
% vindex(lat,lon,id,num,cv,flag1,flag2,flag3,flag4,flag5,index,1);
% 
% %Number of interpolated data points
% length(nonnan(cell2col(num)))  %  867019
% 
% %Sort by decreasing length
% [len,sorter]=sort(cellength(num),'descend');
% vindex(lat,lon,id,num,cv,flag1,flag2,flag3,flag4,flag5,sorter,1);
% 
% id=cellfirst(id);
% clear flags
% flags.flag1=flag1;
% flags.flag2=flag2;
% flags.flag3=flag3;
% flags.flag4=flag4;
% flags.flag5=flag5;
% make gulfdrifters.sgom id num lat lon cv flags 
% %\*************************************************************************

%/*************************************************************************
%Northern Gulf of Mexico drifters
disp('Processing NGOM dataset... ')

%From Paula Re: NGOM drifters
%This is exactly the same as the southern Gulf of Mexico drifters same 
%design (parachute-drogue) at 45m, sampling every hour, not possible to
%know if  drogue is still on. The qc is exactly the same with the same 
%flags. Info about the drifters and how they correlate with winds and 
%geostrophic velocities can be found in our article of 2012.

%Same format as campeche drifters
dirname=[sourcedir '/ngom/qc_data'];
files=findfiles(dirname,'dat'); 

clear id num lat lon cv flag1 flag2 flag3 flag4 flag5
for i=1:length(files)
    data=load([dirname '/' files{i}]);   
    id{i,1}=data(:,1);
    num{i,1}=datenum(data(:,2),data(:,3),data(:,4),data(:,5),data(:,6),0*data(:,6));
    lat{i,1}=vswap(data(:,7),-9999,nan);  
    lon{i,1}=vswap(data(:,8),-9999,nan);  
    cv{i,1}=100*vswap(data(:,9),-9999,nan).*rot(-2*pi/360*data(:,10)+pi/2);
    flag1{i,1}=data(:,11);
    flag2{i,1}=data(:,12);
    flag3{i,1}=data(:,13);
    flag4{i,1}=data(:,14);
    flag5{i,1}=data(:,15);
end
%cell2col(lat,flag4);
%length(find(flag4==1&isfinite(lat)))  %zero 

%--------------------------------------------------------------------------
%Fix two-hour time difference at 3785 in #80
jj=80;
ii=3784;
index=[1:ii ii+2:4646];
id{jj}(index)=id{jj};
num{jj}(index)=num{jj};
lat{jj}(index)=lat{jj};
lon{jj}(index)=lon{jj};
cv{jj}(index)=cv{jj};
flag1{jj}(index)=flag1{jj};
flag2{jj}(index)=flag2{jj};
flag3{jj}(index)=flag3{jj};
flag4{jj}(index)=flag4{jj};
flag5{jj}(index)=flag5{jj};

num{jj}(ii+1)=num{jj}(ii)+1/24;
lat{jj}(ii+1)=nan;
lon{jj}(ii+1)=nan;
cv{jj}(ii+1)=nan;
flag1{jj}(ii+1)=nan;
flag2{jj}(ii+1)=nan;
flag3{jj}(ii+1)=nan;
flag4{jj}(ii+1)=nan;
flag5{jj}(ii+1)=nan;
filled=flag4;
%--------------------------------------------------------------------------
[lat,lon,id,num,cv,filled]=trajextract(lat,lon,id,num,cv,filled,{region1,region2,region3,region4});
matsave gulfdrifters_ngom id num lat lon cv filled
%\*************************************************************************

%/*************************************************************************
%GLAD data 
disp('Processing GLAD dataset... ')

filename=[sourcedir '/glad/GLAD_15min_filtered.dat'];
fid=fopen(filename);
%Forward over header
for i=1:5,tline = fgetl(fid);end
data=textscan(fid,'%7s %f %{yyyy-MM-dd}D %{HH:mm:ss.SSS}D %f %f %f %f %f %f');

data=data(2:end);
[id,num,time,lat,lon,err,u,v,cverr]=deal(data{:});
num=datenum(num)+datenum(time)-floor(datenum(time));
cv=100*(u+1i*v);  %Units of u and v are m/s
cverr=cverr*100;

colbreaks(id,num,lat,lon,cv,err,cverr);
col2cell(id,num,lat,lon,cv,err,cverr);

[lat,lon,id,num,cv,err,cverr]=trajextract(lat,lon,id,num,cv,err,cverr,{region1,region2,region3,region4});
matsave gulfdrifters_glad id num lat lon cv err cverr
%\*************************************************************************

%/*************************************************************************
%LASER drifters
disp('Processing LASER dataset... ')

filename=[sourcedir '/laser/laser_spot_drifters_clean_v15.dat'];
fid=fopen(filename);
%Forward over header
for i=1:5,tline = fgetl(fid);end
%data=textscan(fid,'%2s %f %{yyyy-MM-dd}D %{HH:mm:ss.SSSSSS}D %f %f %f %f %f %f',1);
data=textscan(fid,'%2s %f %{yyyy-MM-dd}D %{HH:mm:ss.SSSSSS}D %f %f %f %f %f %f');

%data=data(2:end);
[str,id,num,time,lat,lon,err,u,v,cverr]=deal(data{:});
drogue=true(size(id));
for i=1:length(str)
    if strcmp(str{i}(1),'U')|| strcmp(str{i}(1),'V')  %these mean drogue off
        drogue(i)=false;
    end
end
num=datenum(num)+datenum(time)-floor(datenum(time));
cv=100*(u+1i*v);  %Units of u and v are m/s
cverr=cverr*00;

colbreaks(id,num,lat,lon,cv,err,cverr,drogue);
col2cell(id,num,lat,lon,cv,err,cverr,drogue);
[lat,lon,id,num,cv,err,cverr,drogue]=trajextract(lat,lon,id,num,cv,err,cverr,drogue,{region1,region2,region3,region4});
matsave gulfdrifters_laser id num lat lon cv err cverr drogue
%\*************************************************************************
 
%/*************************************************************************
%Ocean Circulation Group
disp('Processing OCG dataset... ')

%These are available for download from
%http://ocgweb.marine.usf.edu/drifter_cite.php?ref=125
%After you click on "Agree"

clear id lon lat num cv
dirname=[sourcedir '/ocg/RedTide'];
files=findfiles(dirname,'csv'); 
for i=1:length(files)
    fid = fopen([dirname '/' files{i}], 'rt');
    data = textscan(fid, '%s %s %f %f %f %f %f %f %f','Delimiter',',', 'CollectOutput',1, 'HeaderLines',1);
    fclose(fid);
    id{i}=i+0*data{2}(:,1);
    lat{i}=data{2}(:,1);
    lon{i}=data{2}(:,2);
    %t{i}=data{2}(:,3);
    cv{i}=100*(data{2}(:,6)+1i*data{2}(:,7));   %Units are m/s
    for j=1:length(data{1})
        num{i}(j,1)=datenum([char(data{1}(j,1)) ' ' char(data{1}(j,2))]);
    end
end
N=length(num);

%this has a half-hour time step; the other two are hourly ...
dirname=[sourcedir '/ocg/OilSpill/CoastGuard'];
files=findfiles(dirname,'csv'); 
for i=1:length(files)
    fid = fopen([dirname '/' files{i}], 'rt');
    data = textscan(fid, '%s %s %f %f %f %f %f %f','Delimiter',',', 'CollectOutput',1, 'HeaderLines',1);
    fclose(fid);
    %decimate these internally to one hour
    data{1}=data{1}(1:2:end,:);
    data{2}=data{2}(1:2:end,:);
    id{i+N}=i+0*data{2}(:,1)+N;
    lat{i+N}=data{2}(:,1);
    lon{i+N}=data{2}(:,2);
    cv{i+N}=100*(data{2}(:,5)+1i*data{2}(:,6)); %Units are m/s
    for j=1:length(data{1})
        num{i+N}(j,1)=datenum([char(data{1}(j,1)) ' ' char(data{1}(j,2))]);
    end
end
N=length(num);

dirname=[sourcedir '/ocg/OilSpill/OceanCirculationGroup'];
files=findfiles(dirname,'csv'); 
for i=1:length(files)
    fid = fopen([dirname '/' files{i}], 'rt');
    data = textscan(fid, '%s %s %f %f %f %f %f %f %f','Delimiter',',', 'CollectOutput',1, 'HeaderLines',1);
    fclose(fid);
    id{i+N}=i+0*data{2}(:,1);
    lat{i+N}=data{2}(:,1);
    lon{i+N}=data{2}(:,2);
    %t{i+N}=data{2}(:,3);
    cv{i+N}=100*(data{2}(:,6)+1i*data{2}(:,7));  %Units are m/s
    for j=1:length(data{1})
        num{i+N}(j,1)=datenum([char(data{1}(j,1)) ' ' char(data{1}(j,2))]);
    end
end
N=length(num);
[lat,lon,id,num,cv]=trajextract(lat,lon,id,num,cv,{region1,region2,region3,region4});

matsave gulfdrifters_ocg id num lat lon cv 
%\*************************************************************************

%/*************************************************************************
%SCULP drifters
disp('Processing SCULP dataset... ')
%This data is from Carter Ohlmann, who sent me the files on July 18, 2011
%
%His email:
%
%HI Jonathan,
%
%THere were 3 specific experiments.  I recall the first 2 are explained in Ohlmann and Niiler
%(2005; Progress O.).  The 3rd is mentioned in Ohlmann et al (2001; JGR).  The 3rd is "biased"
%in that eddies observed in altimetry were seeded.  The first two data sets are online 
%(http://www.icess.ucsb.edu/~carter/data/).  I would have to do some digging to find the 3rd.  
%
%At least I got ya something quick....
%-Carter

%From the Readme
%For SCULP_1, and SCULP_2, time is referenced to the beginning of 1992 and 1996, respectively.  
%
%But this is incorrect: Reference from SCULP_1 is 1993, not 1992. 

dirname=[sourcedir '/sculp/sculp1'];
files=findfiles(dirname,'*','include','pos'); 
clear id lon lat t num cv
n=0;
for i=1:length(files)
    data=load([dirname '/' files{i}],'-ascii');
    %fid=fopen(filename);
    %data = textscan(fid, '%f %f %f');
    if ~isempty(data) %&& i~=81 && i~=132 && i~=261  %Skip some bad files
        n=n+1;
        id{n,1}=i+0*data(:,1);
        num{n,1}=data(:,1)+datenum(1993,1,1)-1;
        lon{n,1}=data(:,2);
        lat{n,1}=data(:,3);
    end 
end
% files=findfiles(dirname,'*','include','vel'); 
% n=0;
% for i=1:length(files)
%     data=load([dirname '/' files{i}],'-ascii');
%     %fid=fopen(filename);
%     %data = textscan(fid, '%f %f %f');
%     if ~isempty(data) %&& i~=81 && i~=132 && i~=261  %Skip some bad files
%         n=n+1;
%         cv{n,1}=data(:,2)+1i*data(:,1);
%     end 
% end
%Does this data have gaps?  No
[dt,sigdt,meddt,maxdt]=sampletimes(num);%aresame(dt,maxdt)  
cv=latlon2uv(num,lat,lon);
[lat,lon,id,num,cv]=trajextract(lat,lon,id,num,cv,{region1,region2,region3,region4});
matsave gulfdrifters_sculp1 id num lat lon cv 
%--------------------------------------------------------------------------
dirname=[sourcedir '/sculp/sculp2'];
files=findfiles(dirname,'*','include','pos'); 
clear id lon lat t num cv
n=0;
for i=1:length(files)
    data=load([dirname '/' files{i}],'-ascii');
    if ~isempty(data)  %&& i~=110   %Skip some bad files
        n=n+1;
        id{n,1}=i+0*data(:,1)+399;
        num{n,1}=data(:,1)+datenum(1996,1,1)-1;
        lon{n,1}=data(:,2);
        lat{n,1}=data(:,3);
    end 
end
%Verify some dates from Ohlmann and Niiler (2005), PiO
%figure,plot(cell2col(num))
%hlines(datenum(1992,10,14))
%hlines(datenum(1996,2,5))
%Could be off by one day here. 
[dt,sigdt,meddt,maxdt]=sampletimes(num);%aresame(dt,maxdt)  
cv=latlon2uv(num,lat,lon);
[lat,lon,id,num,cv]=trajextract(lat,lon,id,num,cv,{region1,region2,region3,region4});
matsave gulfdrifters_sculp2 id num lat lon cv 
%\************************************************************************

% Different numbers of position and velocity files!
% files=findfiles(dirname,'*','include','vel');
% for i=1:length(files)
%     idvel(i,1)=str2num(files{i}(15:19));
% end
% files=findfiles(dirname,'*','include','pos');
% for i=1:length(files)
%     idpos(i,1)=str2num(files{i}(15:19));
% end

%/************************************************************************
%DWDE
disp('Processing DWDE dataset... ')
% From Paula:
% So you know, the id of each drifter from each DWDE experiment specifies 
% if it is a CODE drifter (*COD*), a Microstar drifter (*MIC*, in which 
% case flag5=1 means it lost its drogue), or a DORIS (homemade UABC
% drifters, *DOR*) which also have drogues hanging from them but no sensor 
% to detect if it is lost.
%--------------------------------------------------------------------------
load gulfdrifters_DWDEs
num=cell(0,1);
lat=cell(0,1);
lon=cell(0,1);
cv=cell(0,1);
filled=cell(0,1);
drogue=cell(0,1);
type=[];

for i=1:4
   N=length(eval(['gulfdrifters_pgc.dwde.DWDE' int2str(i) '.num']));
   eval(['num(end+1:end+N)=gulfdrifters_pgc.dwde.DWDE' int2str(i) '.num;'])
   eval(['lat(end+1:end+N)=gulfdrifters_pgc.dwde.DWDE' int2str(i) '.lat;'])
   eval(['lon(end+1:end+N)=gulfdrifters_pgc.dwde.DWDE' int2str(i) '.lon;'])
   eval(['cv(end+1:end+N)=gulfdrifters_pgc.dwde.DWDE' int2str(i) '.cv;'])
   eval(['filled(end+1:end+N)=gulfdrifters_pgc.dwde.DWDE' int2str(i) '.flags.flag4;'])
   eval(['drogue(end+1:end+N)=gulfdrifters_pgc.dwde.DWDE' int2str(i) '.flags.flag5;'])
   eval(['type=[type;gulfdrifters_pgc.dwde.DWDE' int2str(i) '.id(:,7:9)];'])
%   eval(['depth(end+1:end+N)=gulfdrifters_pgc.dwde.DWDE' int2str(i) '.depth;'])
end

for i=1:length(drogue)
    drogue{i}=double(~drogue{i});
end

id1=10000+str2num(gulfdrifters_pgc.dwde.DWDE1.id(:,end-3:end));
id2=10000+str2num(gulfdrifters_pgc.dwde.DWDE2.id(:,end-3:end));
id3=10000+str2num(gulfdrifters_pgc.dwde.DWDE3.id(:,end-3:end));
id4=10000+str2num(gulfdrifters_pgc.dwde.DWDE4.id(:,end-3:end));

id1=[id1;id2;id3;id4];
id=num;
for i=1:length(id)
    id{i}=id1(i)+0*id{i};
end

matsave gulfdrifters_dwde id num lat lon cv filled drogue type
% %--------------------------------------------------------------------------
% %reading data from the more recent NetCDF file ... 
% %not using this so I don't have to redo the manual editing portion
% %apparently, there is drogue information from the microstars that I don't use
% %dirname=[sourcedir '/dwde'];
% filename=[sourcedir '/dwde/DRIFTERS_DWDE-CIGoM.nc'];
% lat=ncread(filename,'LAT');
% lon=ncread(filename,'LON');
% id=double(ncread(filename,'Drifters'));
% type=ncread(filename,'Drifter_type')';
% time=datenum(1900,1,1,0,0,0)+double(ncread(filename,'TIME'))/24/60/60;%seconds since 1900-01-01 00:00:00
% u=100*ncread(filename,'U');%convert to cm/s
% v=100*ncread(filename,'V');%convert to cm/s
% row_size=double(ncread(filename,'row_size'));
% filled=double(ncread(filename,'F3'));
% drogue=~double(ncread(filename,'F5')); %Flag5=1 means drogue is lost 
% 
% datestr(min(time)) %this agrees with the documentation
% %time_coverage_start = "21-Jun-2016 19:00:00"
% 
% ids=cell(0,1);
% for i=1:length(row_size)
%     ids{i}=id(i)+zeros(row_size(i),1);
% end
% ids=cell2col(ids);
% ids=ids(isfinite(ids));
%    
% allall(size(ids)==size(lat))  
% vsize(ids,time,lat,lon,u,v,filled,drogue)
% colbreaks(ids,time,lat,lon,u,v,filled,drogue);
% col2cell(ids,time,lat,lon,u,v,filled,drogue);
% length(cell2col(lon)) %416659, same as earlier version
% %and plotting the two datasets, they look identical
%\************************************************************************

%/*************************************************************************
%Splash 
disp('Processing SPLASH dataset... ')
%--------------------------------------------------------------------------
% Splash-raw... looks horrible
%
% Column 1:  drifter ID string (like S_XXXX, XXXX = 4-digit integer)
% Column 2:  date (yyyy-mm-dd)
% Column 3:  time (HH:MM)
% Column 4:  longitude (decimal degrees)
% Column 5:  latitude (decimal degrees)
%
% filename=[sourcedir '/Data/splash/SPLASH_drifters_raw.dat'];
% fid=fopen(filename);
% %Forward over header
% for i=1:5,tline = fgetl(fid);end
% data=textscan(fid,'%2s %f %{yyyy-MM-dd}D %{HH:mm}D %f %f');
% 
% data=data(2:end);
% [id,num,time,lon,lat]=deal(data{:});
% num=datenum(num)+datenum(time)-floor(datenum(time));
% colbreaks(id,num,lat,lon);
% col2cell(id,num,lat,lon);
%--------------------------------------------------------------------------
% Splash-5min
%
% Column 1:  drifter ID string (like S_XXXX, XXXX = 4-digit integer)
% Column 2:  date (yyyy-mm-dd)
% Column 3:  time (HH:MM)
% Column 4:  longitude (decimal degrees)
% Column 5:  latitude (decimal degrees)
% Column 6:  temperature (degrees C; NaN if not available)
% Column 7:  salinity (psu; NaN if not available)
%
% ------------------------------------------------------------------------------
% ASCII file format for the launch info file
% ------------------------------------------------------------------------------
% 
% Column 1:  drifter ID string (like S_XXXX, XXXX = 4-digit integer)
% Column 2:  drifter type
%              drogued     = standard drifters
%              undrogued   = standard drifters but without drogue
%              temperature = standard drifters with temperature (HOBO) sensor
%              salinity    = standard drifters with temperature + salinity (HOBO) sensors
% Column 3:  launch date (yyyy-mm-dd)
% Column 3:  launch time (HH:MM:SS)
% Column 4:  vessel (WS for Walton Smith; the others are small boats ARGUS, Alex,
%            Arnoldo, CB for Calibou Boca, SAND_CRAB, TATIANA, UCLA, and WHISKEY_PASS)
% Column 5:  launch group
%              LSS1a, LSS1b, LSS2           = Large Scale Survey lines
%              DDD1, DDD2                   = 12 x 12 Dense Drifter Deployments,
%                                             spaced nominally 1 km apart
%              Drip1                        = consecutive release at one location
%              SSDD1                        = 5 x 5 Small Scale Drifter Deployment,
%                                             spaced nominally 0.25 km apart
%              SSDD2, SSDD4                 = 5 x 5 Small Scale Drifter Deployments,
%                                             spaced nominally 0.5 km apart
%              SSDD3                        = drogued/undrogued alternating ring
%                                             inside SSDD2
%              Bastian, Caminada, Ronquille = inlet deployments
%              Taylor                       = 2 drifters released near Taylor oil leak
%              Misc                         = 2 additional releases
% Column 6:  launch node (See note below.)
% Column 7:  longitude (decimal degrees)
% Column 8:  latitude (decimal degrees)
% Column 9:  end date (yyyy-mm-dd) of the entire record deemed credible
% Column 10: end time (HH:MM:SS) of the entire record deemed credible
%
% The CARTHE-style drifters (drogued at 0.5 m), tracked in real-time using SPOT
% GPS units, were released in multiple groups, including several
% survey lines and multiple gridded arrays with various spacing, in the 
% Louisiana Bight of the northern Gulf of Mexico in April-May 2017.  Some of the 
% drifters had temperature and salinity sensors attached.  Many drifters 
% (including those with sensors) were released and recaptured for repeated releases.
% 
% The following processing was performed on the raw trajectory data:
%  (1) Duplicate entries, out-of-range positions, and data before launch and after 
%      best guess end time were removed.
%  (2) Positions implying velocities exceeding 2.62 m/s or accelerations exceeding 
%      0.01 m/s^2 were removed.
%  (3) Positions were interpolated to regular 5-minute intervals using a piecewise 
%      Hermitian cubic algorithm, but gaps exceeding 3 hours were retained.
%  (4) Temperature and salinity measurements as obtained from HOBO sensors were 
%      interpolated (piecewise Hermitian cubic) to the same time sampled by the 
%      trajectories and appended to the trajectory data.  Clearly non-sensical 
%      salinity and temperature data were deleted.  Note also that some sensors
%      failed, so that drifters of type ?salinity? or ?temperature? may not in fact
%      have these data appended.

filename=[sourcedir '/splash-5min/SPLASH_drifters_TS.dat'];
fid=fopen(filename);
%Forward over header
for i=1:5,tline = fgetl(fid);end
data=textscan(fid,'%2s %f %{yyyy-MM-dd}D %{HH:mm}D %f %f %f %f');

data=data(2:end);
[id,num,time,lon,lat,t,s]=deal(data{:});
num=datenum(num)+datenum(time)-floor(datenum(time));
colbreaks(id,num,lat,lon,t,s);
col2cell(id,num,lat,lon,t,s);

%reading off from launch file, splash_spot_drifter_launch_info.prn, the IDs
%of drifters that were undrogued at launch
drogued_ids=[523 91 34 26 56 176 43 40 30 46 45 335 314 33 28];
drogue=num;
for i=1:length(num)
    if anyany(drogued_ids==id{i}(1))
        drogue{i}=zeros(size(num{i}));
    else
        drogue{i}=inf+ones(size(num{i}));
    end
end

matsave gulfdrifters_splash id num lat lon t s drogue
%\*************************************************************************

%/*************************************************************************
%Hercules
disp('Processing Hercules dataset... ')
%
% These drifters used two different experimental drifter body designs:
% 
%   type-a:  plastic, tubular body roughly 50 cm high.  
%   type-b:  wood body (approx. 75 cm tall) consisting of four fins 
%            arranged in a plus sign shape.  
% 
% Thirteen type-a (tube style) drifters were launched.
% Six type-b (wood cross style) drifters were launched.
%
% In the position data files, drifters are uniquely identified by
% an ID string, formatted like:
% 
%   HERCULES_[a,b]_XXX,  where 
% 
%    [a,b] = drifter body type (a or b, see above)
%      XXX = three-digit integer ID number
%--------------------------------------------------------------------------
%Type A
filename=[sourcedir '/hercules/carthe_hercules_spot_a.dat'];
fid=fopen(filename);
%Forward over header
for i=1:5,tline = fgetl(fid);end
data=textscan(fid,'%11s %f %{yyyy-MM-dd}D %{HH:mm:ss}D %f %f');
data=data(2:end);
[id,num,time,lat,lon]=deal(data{:});
num=datenum(num)+datenum(time)-floor(datenum(time));
colbreaks(id,num,lat,lon);
col2cell(id,num,lat,lon);
matsave gulfdrifters_hercules id num lat lon
%--------------------------------------------------------------------------
%Type B ... not worth it ... small, poorly sampled, short duration
% filename=[sourcedir '/hercules/carthe_hercules_spot_b.dat'];
% fid=fopen(filename);
% %Forward over header
% for i=1:5,tline = fgetl(fid);end
% data=textscan(fid,'%11s %f %{yyyy-MM-dd}D %{HH:mm:ss}D %f %f');
% data=data(2:end);
% [idb,numb,timeb,latb,lonb]=deal(data{:});
% numb=datenum(numb)+datenum(timeb)-floor(datenum(timeb));
% kindb=nan*numb;
% figure,plot(num,lon,'.'),hold on. plot(numb,lonb,'.')
% id=[id;idb];
% num=[num;numb];
% lat=[lat;latb];
% lon=[lon;lonb];
% kind=[kind;kindb];
% colbreaks(id,num,lat,lon,kind);
% col2cell(id,num,lat,lon,kind);
% %--------------------------------------------------------------------------
%\*************************************************************************

%/*************************************************************************
%AOML
disp('Processing AOML dataset... ')
dirname=[sourcedir '/aoml'];
files=findfiles(dirname,'txt'); 
clear id yr mo dy hr mm ss lat lon t 
for i=1:length(files)
    i
    fid=fopen([dirname '/' files{i}]);
    %Forward over header
    for j=1:33,tline = fgetl(fid);end
    if ~tline(1)=='-'
        data=textscan(fid,'%f %f %f %f %f %f %f %f %f');
        [id{i,1},yr{i,1},mo{i,1},dy{i,1},hr{i,1},mm{i,1},ss{i,1},lat{i},lon{i}]=deal(data{:});
        t{i}=inf*lon{i};
    else
        tline = fgetl(fid);
        data=textscan(fid,'%f %f %f %f %f %f %f %f %f %f');
        [id{i,1},yr{i,1},mo{i,1},dy{i,1},hr{i,1},mm{i,1},ss{i,1},lat{i},lon{i},t{i}]=deal(data{:});
    end
end
    
cell2col(id,yr,mo,dy,hr,mm,ss,lat,lon,t);
t(t<10|t>50)=inf;
num=datenum(yr,mo,dy,hr,mm,ss);
col2cell(id,num,lat,lon,t);
[lat,lon,id,num,t]=trajextract(lat,lon,id,num,t,{region1,region2,region3,region4});


% What is going on with AOML
%[dt,sigdt,meddt,maxdt]=sampletimes(num);
%plot(dt*24) %irregular sampling between 1/2 hour and 3.5 hours
%plot(sort(diff(cell2col(num)))*24),hlines(1),hlines(2)

length(find(sort(diff(cell2col(num)))*24<1))/length((sort(diff(cell2col(num))~=nan)))*100
length(find(sort(diff(cell2col(num)))*24<2))/length((sort(diff(cell2col(num))~=nan)))*100
%about 87 % are less than 2 hrs, 56 % less than 1 hour ... but these are
%very strongly clustered
for i=1:length(num)
        a=min(num{i});
        b=max(num{i});
    if length(a:1/24:b)>2
        [lat{i},xmid,ymid]=twodstats(zeros(size(num{i})),num{i},lat{i},[-1 1],a:1/24:b);
        [lon{i},xmid,ymid]=twodstats(zeros(size(num{i})),num{i},lon{i},[-1 1],a:1/24:b);
        [t{i},xmid,ymid]=twodstats(zeros(size(num{i})),num{i},t{i},[-1 1],a:1/24:b);
        num{i}=ymid;
        id{i}=id{i}(1)+zeros(size(num{i}));
        lon{i}=vswap(lon{i},nan,inf);
        lat{i}=vswap(lat{i},nan,inf);
        t{i}=vswap(t{i},nan,inf);
    else
        num{i}=[];
        lat{i}=[];
        lon{i}=[];
        id{i}=[];
        t{i}=[];
    end
end
cellprune(id,num,lat,lon,t);
%how big are the gaps
xx=isinf(cell2col(lat));
[N,a,b]=blocknum(xx);
%allall(xx(a(1:2:end))==1)
linfs=b(1:2:end)-a(1:2:end)+1;
%linear interpolation of up to six hours
for i=1:length(num)
    i
    lon{i}=fillbad(lon{i},inf,6);
    lat{i}=fillbad(lat{i},inf,6);
    t{i}=fillbad(t{i},inf,6);
end

%plot(sort(linfs)/24),hlines(2),hlines(1),hlines(1/2)
%next, interpolate
%most (99.5%) of the gaps are less than 1/2 day
length(find(linfs/24<1/2))./length(linfs)

matsave gulfdrifters_aoml id num lat lon t
%\*************************************************************************

% %/*************************************************************************
% %Taylor
% disp('Taylor dataset... ')
% %possibly not worth it.  Also, no idea what this is
% 
% filename=[sourcedir '/taylor/Taylor_drifters_all_15min_Apr_2017.txt'];
% fid=fopen(filename);
% data=textscan(fid,'%f %f %f %f %2s');
% %Trajectories from 5 undrogued and 7 drogued drifters were tracked in real-time using SPOT GPS units, 
% %launched near the Taylor Energy site just offshore of the Mississippi Delta on April 18-20, 2017 
% %as a collaboration of two GoMRI projects: ?Influence of river induced fronts on hydrocarbon transport?
% %and ?CARTHE?. Positions are interpolated to uniform 15 minute intervals starting from launching time
% %through early June.
% [lon,lat,time,id,str]=deal(data{:});
% flag=true(size(lon));
% for i=1:length(lon)
%     if strcmpi(str{i},'ud')
%         flag(i)=false;
%     end
% end
% num=datenum(2017,1,1)-1+time;
% colbreaks(id,num,lat,lon,flag);
% col2cell(id,num,lat,lon,flag);
% make gulfdrifters_taylor id num lat lon flag
% %cellplot(lon)  %this looks horrible
% %\*************************************************************************

%/*************************************************************************
%LATEX
disp('LATEX dataset... ')

dirname=[sourcedir '/latex'];
files=findfiles(dirname,'drf'); 
clear num id mo dy yr hr mm ss lat lon t
for i=1:length(files)
    i
    fid=fopen([dirname '/' files{i}]);
    %Forward over header
    bdone=false;
    while ~bdone
        tline = fgetl(fid);
        if ~isempty(tline)
            if strcmpi(tline(1:5),'*END*')
                bdone=true;
               % tline = fgetl(fid);
            end
        end
    end
    data=textscan(fid,'%f %1s %f %1s %f %f %1s %f %1s %f %f %1s %f %1s %f');
    [mo,~,dy,~,yr,hr,~,mm,~,ss,lat{i,1},~,lon{i,1},~,t{i,1}]=deal(data{:});
    num{i,1}=datenum(yr,mo,dy,hr,mm,ss);
    lon{i}=-lon{i};
    name=files{i}(1:end-4);
    if name(end)=='a'
           name(end)='1';
    elseif name(end)=='b'
           name(end)='1';
    end
    
    id{i,1}=str2num(name)+0*lon{i};
end

[lat,lon,id,num]=...
    trajextract(lat,lon,id,num,{region1,region2,region3,region4});
matsave gulfdrifters_latex id num lat lon 
%\*************************************************************************

% Verif
% splash=gulfdrifters_splash;
% hercules=gulfdrifters_hercules;
% aoml=gulfdrifters_aoml;
% 
% figure,
% use splash
% cellplot(lon,lat,'2g'),hold on
% use hercules
% cellplot(lon,lat,'2r'),hold on
% use aoml
% cellplot(lon,lat,'2b'),hold on
% use gulfdrifters
% cellplot(lon,lat,'0.1k')

function[]=gulfdrifters_primary

load gulfdrifters_latex
load gulfdrifters_sculp1
load gulfdrifters_sculp2
load gulfdrifters_gdp
load gulfdrifters_hargos
load gulfdrifters_hgps
load gulfdrifters_aoml
load gulfdrifters_sgom
load gulfdrifters_ngom
load gulfdrifters_ocg
load gulfdrifters_glad
load gulfdrifters_hercules
load gulfdrifters_laser
load gulfdrifters_splash
load gulfdrifters_dwde

clear gulfdrifters_experiments

gulfdrifters_experiments.latex=gulfdrifters_latex;
gulfdrifters_experiments.sculp1=gulfdrifters_sculp1;
gulfdrifters_experiments.sculp2=gulfdrifters_sculp2;
gulfdrifters_experiments.gdp=gulfdrifters_gdp;
gulfdrifters_experiments.hargos=gulfdrifters_hargos;
gulfdrifters_experiments.aoml=gulfdrifters_aoml;
gulfdrifters_experiments.sgom=gulfdrifters_sgom;
gulfdrifters_experiments.ngom=gulfdrifters_ngom;
gulfdrifters_experiments.ocg=gulfdrifters_ocg;
gulfdrifters_experiments.glad=gulfdrifters_glad;
gulfdrifters_experiments.hercules=gulfdrifters_hercules;
gulfdrifters_experiments.hgps=gulfdrifters_hgps;
gulfdrifters_experiments.laser=gulfdrifters_laser;
gulfdrifters_experiments.dwde=gulfdrifters_dwde;
gulfdrifters_experiments.splash=gulfdrifters_splash;

%gulfdrifters=original_gulfdrifters;
%/************************************************************************
%Visual quality control
names=fieldnames(gulfdrifters_experiments);
%--------------------------------------------------------------------------
%Add NaNs for missing fields
varnames={'cv','t','drogue','gap','filled'};
for i=1:length(names)
    i
    for j = 1:length(varnames)
        clear(varnames{j})
    end
    eval(['use gulfdrifters_experiments.' names{i}])
    for j=1:length(varnames)
        if exist(varnames{j})~=1
            eval([varnames{j} '=num;']);
            for k=1:length(num)
                eval([varnames{j} '{k}=nan*num{k};']);
            end
        end
        eval(['make  gulfdrifters_experiments.' names{i} ' ' varnames{j}])
    end
end
%--------------------------------------------------------------------------
%Add depth to all of these
[topo,lato,lono]=readtopo([-99   -80    18    31]);
lat=[];lon=[];%Have to initialize these to keep Matlab from complaining
for i=1:length(names)
    i
    eval(['use gulfdrifters_experiments.' names{i}])
    depth=lat;
    parfor j=1:length(lon)
        depth{j}=interp2(lato,lono,-topo',lat{j},lon{j});
    end
    eval(['make gulfdrifters_experiments.' names{i} ' depth'])
end
%--------------------------------------------------------------------------
%Adding drifter type
for i=1:length(names)
    eval(['N=length(gulfdrifters_experiments.' names{i} '.id);'])
    if i==14  %for dwde I'm keeping track of type
        %do nothing
    else
        switch names{i}
            case 'latex'
                str='WOC';
            case {'sculp1','sculp2','aoml','ocg','glad'}
                str='COD';
            case {'gdp','hargos','hgps'}
                str='SVP';
            case {'sgom','ngom'}
                str='FHD';
            case {'hercules'}
                str='TUB';
            case {'laser','splash'}
                str='CAR';
        end
        type=vrep(str,N,1);
        eval(['gulfdrifters_experiments.' names{i} '.type=type;'])
    end
end
%--------------------------------------------------------------------------

% for i=1:length(names)
%     eval(['use gulfdrifters_experiments.' names{i}])
%   %  eval(['use original_gulfdrifters_experiments.' names{i}])
%     [dt,sigdt,meddt,maxdt]=sampletimes(num);
%     disp([int2str(i) ' ' names{i} ':' num2str(median(meddt*24)) ' ' num2str(maxmax(dt*24))  ])
%     %Note OCG has a mixture of hourly and half-hourly
% end
% 1 latex:6 6
% 2 sculp1:1.5 1.5
% 3 sculp2:1.5 1.5
% 4 gdp:6 6
% 5 hargos:1 1.6136
% 6 aoml:1 1
% 7 sgom:1 1.0083
% 8 ngom:1 3.101
% 9 ocg:1 1.9121
% 10 glad:0.25 0.25
% 11 hercules:0.083333 0.33189
% 12 hgps:1 1.1738
% 13 laser:0.25 0.33095
% 14 dwde:1 1.0007
% 15 splash:0.083333 0.31442
%
% clear mint
% for i=1:length(names)   
%     eval(['use gulfdrifters_experiments.' names{i}])
%     mint(i)=minmin(cell2col(num));
% end
%minmin(diff(mint))   %ok, in chronological order

% figure
% for i=1:length(names)
%     i
%     eval(['use gulfdrifters_experiments.' names{i}])
%     plot(sort(cell2col(depth))),hold on
% end
% figure
% for i=1:length(names)
%     i
%     eval(['use gulfdrifters_experiments.' names{i}])
%     plot(sort(abs(cell2col(cv)))),hold on
% end
% 
% for i=1:length(names)
%     i
%     eval(['use gulfdrifters_experiments.' names{i}])
%     numer=length(find(abs(cell2col(cv))<1/100));
%     denom=length(find(isfinite(cell2col(cv))));
%     [i, numer]
%     numer./denom*100
% end
% 
% 
% %there are still some gaps, so, shouldn't compute velocity yet
% for i=1:length(names)
%     i
%     eval(['use gulfdrifters_experiments.' names{i}])
%     [dt,sigdt,meddt,maxdt]=sampletimes(num);
%     [i,mean(dt),maxmax(maxdt)]
% end

id=[];%Need to inform Matlab that id is a variable or it gets confused
%with a function on my path called id
for i=1:length(names)
    i
    eval(['use gulfdrifters_experiments.' names{i}])
    %figure,cellplot(lon,'r'),hold on
    [dt,sigdt,meddt,maxdt]=sampletimes(num);
    %--------------------------------------------------------------------------
    %Set drifters on less than 0 m depth to infs ...
    for j=1:length(lon)
        lon{j}(depth{j}<0)=inf;
        %lat{j}(depth{j}<10/1000)=inf;
        %lon{j}(depth{j}<10/1000)=inf;
    end
    %Manual inspection
    %cellplot(lon)
%     ['gulfdrifters_experiments.' names{i}]
%     gulfdrifters_experiments.latex
%     size(id)
%     iscell(id)
%    return
    id1=cellfirst(id);
    switch i
         case 1 %LATEX
            lon{10}(354:415)=inf;
            lon{19}(1:5)=inf;
         case 2 %SCULP1
            lon{5}(:)=inf;
            lon{6}(720:end)=inf;
            lon{74}(:)=inf;
            lon{81}(:)=inf;
            lon{82}(:)=inf;
            lon{83}(:)=inf;
            lon{84}(:)=inf;
            lon{85}(:)=inf;
            lon{86}(:)=inf;
            lon{87}(:)=inf;
            lon{118}(:)=inf;
            lon{135}(:)=inf;
            lon{138}(1400:end)=inf;
            lon{145}(:)=inf;
            lon{158}(200:340)=inf;
            lon{164}(:)=inf;
            lon{180}(648:end)=inf;
            lon{265}(:)=inf;
            lon{293}(:)=inf;
            lon{300}(:)=inf;
            lon{306}(392:end)=inf;
            lon{308}(:)=inf;
            lon{318}(156:end)=inf;
            lon{320}(:)=inf;
            lon{328}(710:end)=inf;
            lon{336}(350:end)=inf;
            lon{341}(160:end)=inf;
            lon{348}(900:end)=inf;
            lon{349}(:)=inf;
            lon{351}(1:12)=inf;
            lon{351}(600:end)=inf;
            lon{357}(:)=inf;
            lon{358}(400:end)=inf;
            lon{364}(90:end)=inf;
            lon{365}(330:end)=inf;
            lon{366}(:)=inf;
            lon{368}(1:65)=inf;
            lon{369}(845:end)=inf;
            lon{370}(1:118)=inf;            
            lon{370}(:)=inf;            
            lon{375}(920:end)=inf;            
            lon{380}(1:300)=inf;            
            lon{383}(1:10)=inf;            
            lon{393}(:)=inf;            
            lon{395}(:)=inf;            
            lon{399}(460:end)=inf;  
            %examining some suspicious segments based on preliminary maps
            %ids 94 39 and 226 in SCULP1
            %ii=find(cellfirst(id)== 39);
            ii=find(cellfirst(id)== 94);lon{ii}(674:730)=inf;
            %ii=find(cellfirst(id)== 226);
        case 3 %SCULP2
            lon{5}(:)=inf;            
            lon{6}(1:680)=inf;            
            lon{10}(:)=inf;            
            lon{19}(190:end)=inf;            
            lon{33}(270:end)=inf;            
            lon{42}(:)=inf;            
            lon{44}(:)=inf;            
            lon{45}(:)=inf;            
            lon{59}(:)=inf;            
            lon{65}(:)=inf;            
            lon{66}(:)=inf;            
            lon{90}(1525:end)=inf;            
            lon{106}(:)=inf;            
            lon{110}(2000:end)=inf;            
            lon{115}(600:end)=inf;            
            lon{124}(:)=inf;            
            lon{126}(1500:end)=inf;            
            lon{131}(1580:end)=inf;            
            lon{142}(:)=inf;            
            lon{143}(550:end)=inf;            
            lon{148}(1060:end)=inf;            
            lon{165}(460:end)=inf;            
            lon{176}(1300:1550)=inf; 
            lon{178}(:)=inf;            
            lon{212}(780:end)=inf;            
            lon{214}(:)=inf;            
            lon{215}(:)=inf;            
            lon{219}(1350:1550)=inf;            
            lon{224}(135:end)=inf;            
            lon{226}(694:end)=inf;            
            lon{244}(240:end)=inf;            
            lon{247}(250:end)=inf;            
            lon{254}(730:end)=inf; 
            lon{256}(:)=inf; 
            %examining some suspicious segments based on preliminary maps
            %ids 406 in SCULP2
            ii=find(cellfirst(id)== 406);lon{ii}(1:850)=inf;
        case 4  %GDP
            lon{49}(99:end)=inf;
            lon{54}(765:end)=inf;
            lon{55}(1:191)=inf;
            lon{59}(173:end)=inf;
            lon{66}(:)=inf;
         case 5 %HARGOS
            %--------------------------------------------------------------------------
            %Correct some obviously visually bad sections or trajectories .. flatlining
            ii=find(cellfirst(id)== 57941);  lon{ii}(253:end)=inf;
            ii=find(cellfirst(id)== 70328);  lon{ii}(975:2724)=inf;
            %ii=find(cellfirst(id)== 132538); lon{ii}(:)=inf;
            ii=find(cellfirst(id)== 98915);  lon{ii}(9500:9700)=inf;
            ii=find(cellfirst(id)== 98916);  lon{ii}(5050:5320)=inf;
        case 6 %AOML
%             lon{1}(:)=inf;
%             lon{3}(:)=inf;
%             lon{5}(:)=inf;
%             lon{7}(:)=inf;
%             lon{18}(:)=inf;
            lon{26}(770:end)=inf;
            lon{37}(972:975)=inf;
            lon{74}(575:800)=inf;
            lon{65}(2162:end)=inf;
        case 7  %SGOM
            ii=find(cellfirst(id)== 80097); lon{ii}(1490:end)=inf;
            ii=find(cellfirst(id)== 80100); lon{ii}(2200:end)=inf;
            %examining some suspicious segments based on preliminary maps
            %ids 80253 in SGOM | ids 80125 80216 80224 in SGOM
            ii=find(cellfirst(id)== 80253);lon{ii(1)}(620:end)=inf;
            ii=find(cellfirst(id)== 80125);lon{ii}(140:end)=inf;
            %ii=find(cellfirst(id)== 80216);
            %ii=find(cellfirst(id)== 80224);
        case 8 %NGOM 
            lon{80}(1050:1633)=inf;
            lon{85}(836:1585)=inf;
            lon{110}(1391:2097)=inf;
            lon{172}(1528:1826)=inf;
            lon{177}(758:959)=inf;
            lon{193}(3495:4033)=inf;
            lon{201}(807:1442)=inf;
            lon{241}(2466:2857)=inf;
            lon{242}(1436:end)=inf;
            lon{250}(1790:2597)=inf;
            lon{272}(617:1489)=inf;
            lon{285}(880:1609)=inf;
            lon{300}(881:1681)=inf;
            lon{305}(4180:5240)=inf;
            lon{335}(470:end)=inf;
            lon{336}(2153:2401)=inf;
            lon{339}(1096:1562)=inf;
            lon{369}(97:end)=inf;
            lon{370}(2030:2257)=inf;
        case 9 %OCG
            lon{1}(86)=inf;
            lon{17}(1:2)=inf;
            lon{26}(1:2)=inf;
            lon{46}(640)=inf;
            lon{48}(360)=inf; 
            lon{50}(1226:end)=inf; 
        case 10 %GLAD 
            lon{24}(3630:end)=inf;
            lon{42}(3725:end)=inf;
            lon{79}(3645:end)=inf;
            lon{85}(3593:end)=inf;
            lon{91}(3631:end)=inf;
            lon{100}(3605:end)=inf;
            lon{123}(3229:end)=inf;
            lon{130}(3799:end)=inf;
            lon{150}(3228:end)=inf;
            lon{159}(3990:end)=inf;
            lon{177}(3249:end)=inf;
            lon{180}(3241:end)=inf;
            lon{194}(3242:end)=inf;
            lon{294}(1670:end)=inf;
            %examining some suspicious segments based on preliminary maps
            %ids 101 and 317 in GLAD
            %ii=find(cellfirst(id)== 101);
            ii=find(cellfirst(id)== 317);lon{ii}(1:32)=inf;
       case 11 %Hercules
            lon{1}(:)=inf;
            lon{2}(3650:end)=inf;
            lon{3}(4065:end)=inf;
            lon{4}(3670:3685)=inf;
            lon{4}(4354:end)=inf;
            lon{5}(3680:3900)=inf;
            lon{6}(2350:end)=inf;
            lon{7}(4818:end)=inf;
            lon{8}(7170:end)=inf;
            lon{9}(3780:end)=inf;
            lon{10}(7788:end)=inf;
            lon{11}(9470:end)=inf;
            lon{12}(3180:end)=inf;
            lon{13}(3700:end)=inf;
       case 12 %HGPS  
            lon{33}(2085:end)=inf;
            lon{37}(310:end)=inf;
            %examining some suspicious segments based on preliminary maps
            % ids 63701930
            %ii=find(cellfirst(id)==  63701930)
       case 13 %LASER %Looks fine
       case 14 %DWDE
             lon{100}(190:282)=inf;
             lon{103}(800:860)=inf;
             lon{103}(1780:1940)=inf;
             lon{116}(1700:1950)=inf;
             lon{138}(1900:1970)=inf;
             lon{138}(2110:2288)=inf;
             lon{141}(1529:1708)=inf;
             lon{144}(1897:2082)=inf;
             lon{148}(1370:1517)=inf;
             lon{148}(1560:end)=inf;
             lon{149}(239:end)=inf;
             lon{160}(1490:1609)=inf;
             lon{160}(2129:end)=inf;
             lon{165}(1305:1351)=inf;
             lon{165}(1430:1466)=inf;
             lon{165}(1512:1599)=inf;
             lon{170}(1326:1386)=inf;
             lon{170}(1435:end)=inf;
             lon{171}(1440:end)=inf;
             lon{177}(1349:1375)=inf;
             lon{177}(1449:1506)=inf;
             lon{181}(1275:end)=inf;
             lon{182}(1540:end)=inf;
             lon{183}(1888:2128)=inf;
             lon{185}(1012:end)=inf;
             lon{188}(2057:end)=inf;
             lon{190}(1154:end)=inf;
             lon{198}(1415:end)=inf;
             lon{199}(1551:end)=inf;
             lon{204}(2729:end)=inf;
             %examining some suspicious segments based on preliminary maps
             %ids 10075 10087 18320 10082 10029 17060 10026
             ii=find(cellfirst(id)== 10075);lon{ii}(600:end)=inf;
             ii=find(cellfirst(id)== 10087);lon{ii}(783:end)=inf;
             ii=find(cellfirst(id)== 18320);lon{ii}(490:end)=inf;
             %ii=find(cellfirst(id)== 10082); 
             ii=find(cellfirst(id)== 10029); 
             ii=find(cellfirst(id)== 17060);lon{ii(1)}(566:end)=inf;
             %ii=find(cellfirst(id)== 10026);%ok
        case 15 %SPLASH   %looks ok ... some obvious stationary regions we can flag later
             %examining some suspicious segments based on preliminary maps
             %ids 307
             %ii=find(cellfirst(id)==307);
    end 
%    eval(['make gulfdrifters_experiments.' names{i} ' id num lat lon cv t drogue gap laterr lonerr cverr depth err' ])
    eval(['make gulfdrifters_experiments.' names{i} ' id type num lat lon cv t drogue gap depth filled' ])
end
%\*************************************************************************


%/*************************************************************************
%convert type to a number
for i=1:length(names)
    eval(['id=gulfdrifters_experiments.' names{i} '.id;'])
    eval(['type=gulfdrifters_experiments.' names{i} '.type;'])
    typenum=zeros(size(id));
    for j=1:size(type,1)
        switch type(j,:)
            case 'WOC'
                typenum(j)=1;
            case 'COD'
                typenum(j)=2;
            case 'SVP'
                typenum(j)=3;
            case 'FHD'
                typenum(j)=4;
            case 'TUB'
                typenum(j)=5;
            case 'CAR'
                typenum(j)=6;
            case 'MIC'
                typenum(j)=7;
            case 'DOR'
                typenum(j)=8;
        end
    end
    eval(['gulfdrifters_experiments.' names{i} '.typenum=typenum;'])
end
%\*************************************************************************


for i=1:length(names)
    names{i}
    eval(['gulfdrifters_experiments.' names{i}])
end


clear gulfdrifters_experiments_primary
%/*************************************************************************
%Primary processing
for i=1:length(names)           
     
    i
    num=[];lat=[];lon=[];%Have to initialize these to keep Matlab from complaining
    eval(['use gulfdrifters_experiments.' names{i}])
    %eval(['gulfdrifters_experiments.' names{i}])
    %--------------------------------------------------------------------------
    %Remove all data points marked with infs
    cellpack(lon,lat,num,id,cv,gap,drogue,filled,t,depth);
    %--------------------------------------------------------------------------
    %Split at gaps
    %[dt,sigdt,meddt,maxdt]=sampletimes(num);
       
    %vsize(num,lat,lon,id,laterr,lonerr,cv,cverr,gap,drogue,err,t);
    %cellsplit(num,lat,lon,id,laterr,lonerr,cv,cverr,gap,drogue,err,t,1); %Split with longest gap at 1 days
    %vsize(num,lat,lon,id,laterr,lonerr,cv,cverr,gap,drogue,err,t);
    
    %[dt,sigdt,meddt,maxdt]=sampletimes(num);
    [dt,sigdt,meddt,maxdt]=sampletimes(cell2col(num));
    %--------------------------------------------------------------------------
    %remove segments less than 1 day in length
    vsize(id,typenum,num,lat,lon,cv,gap,drogue,filled,t,depth);
    %vindex(id,num,lat,lon,laterr,lonerr,cv,cverr,gap,drogue,err,t,(cellmax(num)-cellmin(num))>=7,1);
    vindex(id,typenum,num,lat,lon,cv,gap,drogue,filled,t,depth,(cellmax(num)-cellmin(num))>=1,1);
    vsize(id,typenum,num,lat,lon,cv,gap,drogue,filled,t,depth);
    %--------------------------------------------------------------------------
    %interpolate to a uniform hourly time base
    numo=num;
    wasfilled=filled;
    [dt,sigdt,meddt,maxdt]=sampletimes(num);
    
    cellgrid(num,lat,lon,id,cv,gap,drogue,t,depth,1/24);
    nearest=num;
    parfor j=1:length(num)
        %j
        numoj=numo{j}(wasfilled{j}~=1);
        nearest{j}=min(abs(vrep(num{j},size(numoj,1),2)...
            -vrep(numoj',size(num{j},1),1)),[],2);
    end
    filled=num;%overwriting filled with a new hourly filled value
    for j=1:length(num)
        %true if filled
        %filled{j}=nearest{j}>(2/3)*meddt(j);
        filled{j}=nearest{j}>(6+5/60)/24; %flag as filled for gaps >  6hrs 5min
    end
    filled{j}=filled{j}|gap{j}>6;  %for when we have gap, which is just Shane's
    %--------------------------------------------------------------------------
    %recompute velocity
    cv=latlon2uv(num,lat,lon);
    %--------------------------------------------------------------------------
    %remove those containing no valid data
    bool=sum(isfinite(col2mat(cell2col(lat))),1)>0;
    vindex(id,typenum,num,lat,lon,cv,gap,drogue,t,filled,depth,bool,1);
    %--------------------------------------------------------------------------
    %remove more than 10% filled
    %fillfrac=cellsum(filled)./cellength(filled);
    %vindex(id,num,lat,lon,laterr,lonerr,cv,cverr,gap,drogue,err,t,filled,fillfrac<0.1,1);
    %--------------------------------------------------------------------------
    %Sort by decreasing length
    %[len,sorter]=sort(cellength(num),'descend');
    %vindex(id,num,lat,lon,laterr,lonerr,cv,cverr,gap,drogue,err,t,depth,filled,sorter,1);
    
    eval(['make gulfdrifters_experiments_primary.' names{i} ' id typenum num lat lon cv t drogue gap depth filled' ])
end
matsave gulfdrifters_experiments_primary
%\*************************************************************************

%/*************************************************************************
%Merging into one dataset
names=fieldnames(gulfdrifters_experiments_primary);

for i=1:length(names)
    eval(['use gulfdrifters_experiments_primary.' names{i}])
    source=id;
    for j=1:length(id)
        source{j}=i+0*source{j};
    end
    eval(['make gulfdrifters_experiments_primary.' names{i} ' source;'])
end

varnames={'id','typenum','source','num','lat','lon','cv','t','depth','drogue','gap','filled'};

for j = 1:length(varnames)
    eval([varnames{j} '=[];'])
end

for i=1:length(names)
    i
    %eval(['use gulfdrifters.' names{i}])
    for j=1:length(varnames)
        if strcmpi(varnames{j},'typenum')
          eval([varnames{j} '=[' varnames{j} ';' ...
              'gulfdrifters_experiments_primary.' names{i} '.' varnames{j} '];'])
        else            
          eval([varnames{j} '=[' varnames{j} ';' ...
              'cell2col(gulfdrifters_experiments_primary.' names{i} '.' varnames{j} ')];'])
        end
    end
end
gulfdrifters.names=names;

for j=1:length(varnames)
    if ~strcmpi(varnames{j},'typenum')
        eval([varnames{j} '=col2cell(' varnames{j} ');'])
    end
    eval(['make gulfdrifters_primary ' varnames{j}])
end
matsave gulfdrifters_primary
%\*************************************************************************

function[]=gulfdrifters_secondary

gulfdrifters_primary.lat=[];
gulfdrifters_secondary.lat=[];
gulfdrifters_primary.filled=[];
gulfdrifters_secondary.filled=[];

load gulfdrifters_primary
varnames=fieldnames(gulfdrifters_primary);
%/*************************************************************************
%Secondary processing
%--------------------------------------------------------------------------
%Set apparently stuck or motionless drifters to infs ...
use gulfdrifters_primary
cell2col(lat,lon,source,cv,depth,filled);
% for i=1:13
%    index=find(source==i);
%    plot([1:length(index)]./length(index),abs(sort(cv(index)))),hold on
% end
% for i=1:13
%     [n,x]=hist(log10(abs(sort(cv(source==i)))),-2:.025:3);
%     plot(x,n./sum(n)),hold on
% end
bool=abs(cv)<1/10;
%figure,plot(lon,lat,'.'),hold on,plot(lon(bool),lat(bool),'mo')
%figure,plot(lon(bool),lat(bool),'.')
length(find(bool))./length(find(isfinite(lon)))*100  %0.11  percent
lon(bool)=inf;
filled(bool)=true;
%--------------------------------------------------------------------------
%removing very high speed drifters from everywhere
%figure
bool=abs(cv)>250;
%plot(lon(bool),lat(bool),'ro'),hold on,eval(tweakmap)
length(find(bool))./length(find(isfinite(lon)))*100  %0.007 percent
lon(bool)=inf;
filled(bool)=true;
%--------------------------------------------------------------------------
%specifically addressing some problems I see in Sculp 1&2
%>150 cms inside of 1000m isobath, excluding region around the mississippi
%outflow where high speeds are observed in other experiments
region=[-90.034151482582089 -88.668331789970395  28.341654151112468 29.375137894780067];
%bool=abs(cv)>150&depth<1&~(source==2|source==3);
bool=abs(cv)>150&depth<1&(source==2|source==3)&~inregion(region,lat,lon);
%plot(lon(bool),lat(bool),'ro'),hold on,eval(tweakmap)
length(find(bool))./length(find(isfinite(lon)&(source==2|source==3)))*100  %0.1 percent of SCULP
lon(bool)=inf;
filled(bool)=true;
%--------------------------------------------------------------------------
%specifically addressing some problems I see in Sculp 1&2
%>100 cms inside of 7m isobath, also excluding region around the mississippi
bool=abs(cv)>100&depth<7/1000&(source==2|source==3)&~inregion(region,lat,lon);
%plot(lon(bool),lat(bool),'m+'),
length(find(bool))./length(find(isfinite(lon)&(source==2|source==3)))*100  %0.006 percent of SCULP
lon(bool)=inf;
filled(bool)=true;
%--------------------------------------------------------------------------
col2cell(lat,lon,source,cv);
%--------------------------------------------------------------------------
%Remove also low and high acceleration portions
use gulfdrifters_primary
[~,ca]=latlon2uv(num,lat,lon);
cabar=ca;
bool=ca;
N=floor(7*24);%That's a weeklong filter
for j=1:length(ca)
 %   cabar{j}=abs(vfilt(ca{j},N));%weeklong filter
    cabar{j}=vfilt(abs(ca{j}),N);%weeklong filter
end
cell2col(num,lat,lon,source,cv,ca,cabar,filled);

ca(isnan(lon))=nan;
cabar(isnan(lon))=nan;

% figure
% lon1=lon;lat1=lat;num1=num;
% %bool=abs(ca)<1e-2;
% bool=abs(cabar)<1e-2|isinf(cabar);
% lon1(bool)=inf;lat1(bool)=inf;
% cellplot(col2cell(lon1),col2cell(lat1))
% eval(tweakmap_thin)
% bool=abs(cabar)>10^(-2.15)&isfinite(cabar);
% % length(find(abs(bool)))./length(find(isfinite(lom)))*100  %0.063  percent
% %lon(bool)=inf;
% plot(lon(bool),lat(bool),'.')
 
% figure
% lon1=lon;lat1=lat;
% %bool=abs(ca)>1e-4;
% bool=abs(cabar)>1e-4|isinf(cabar);
% lon1(bool)=inf;lat1(bool)=inf;
% cellplot(col2cell(lon1),col2cell(lat1))
% %This one definitely looks good
% %figure,cellplot(col2cell(lon1))
%eval(tweakmap_thin)
bool=abs(cabar)<10^(-4);
%bool=abs(ca)<10^(-7);
%plot(lon(bool),lat(bool),'.')
length(find(abs(bool)))./length(find(isfinite(lon)))*100  %0.63  percent
lon(bool)=inf;
filled(bool)=true;

col2cell(num,lat,lon,source,cv,ca,cabar,filled);

% for i=1:13
%    index=find(source==i);
%    plot([1:length(index)]./length(index),abs(sort(ca(index)))),hold on
% end
% for i=1:13
%    index=find(source==i);
%    plot([1:length(index)]./length(index),abs(sort(cabar(index)))),hold on
% end
%--------------------------------------------------------------------------
% I don't like this way anymore
% cabar=ca;
% bool=ca;
% N=floor(7*24);%That's a weeklong filter
% for j=1:length(ca)
%     cabar{j}=abs(vfilt(ca{j},N));%weeklong filter
%     boolj=cabar{j}<1e-6;%find accelerations less than cutoff
%     L=blocklen(boolj);
%     boolj=boolj&(L>1/7*N);%find low-acceleration blocks more than one day long
%     [L,ia,ib]=blocklen(boolj);
%     bool{j}=false(size(bool{j}));
%     for k=1:length(ia)%lengthen these blocks by 3.5 days in both direction
%         if boolj(ia(k))
%             ia(k)=max(ia(k)-floor(N/2),1);
%             ib(k)=min(ib(k)+floor(N/2),length(boolj));
%             bool{j}(ia(k):ib(k))=true;
%         end
%     end
% end
% 
% 
% for j=1:length(ca)
%     lat{j}(bool{j})=inf;
%     lon{j}(bool{j})=inf;
% end
%--------------------------------------------------------------------------
% Dont' actually need to do all this stuff
% %Flag as missing if we already filled it before, since I don't want to
% %interpolate from previously interpolated points
% %Remove all data points marked with infs
% cell2col(lon,filled);lon(find(filled==1))=inf;col2cell(lon,filled);
% cellpack(lon,lat,num,id,source,laterr,lonerr,cv,cverr,gap,drogue,err,t,depth,filled);
% %--------------------------------------------------------------------------
% %remove segments less than 1 day in length
% vsize(id,source,num,lat,lon,laterr,lonerr,cv,cverr,gap,drogue,err,t,depth,filled);
% %vindex(id,num,lat,lon,laterr,lonerr,cv,cverr,gap,drogue,err,t,(cellmax(num)-cellmin(num))>=7,1);
% vindex(id,source,num,lat,lon,laterr,lonerr,cv,cverr,gap,drogue,err,t,depth,filled,(cellmax(num)-cellmin(num))>=1,1);
% vsize(id,source,num,lat,lon,laterr,lonerr,cv,cverr,gap,drogue,err,t,depth,filled);
% %--------------------------------------------------------------------------
% %interpolate to a uniform hourly time base
% numo=num;
% [dt,sigdt,meddt,maxdt]=sampletimes(num);
% cellgrid(num,lat,lon,id,source,laterr,lonerr,cv,cverr,gap,drogue,err,t,depth,filled,1/24);
% nearest=num;
% parfor j=1:length(num)
%     %j
%     nearest{j}=min(abs(vrep(num{j},size(numo{j},1),2)...
%         -vrep(numo{j}',size(num{j},1),1)),[],2);
% end
% filled=num;
% for j=1:length(num)
%     %true if filled
%     filled{j}=nearest{j}>(6+5/60)/24; %flag as filled for gaps >  6hrs 5min
% end
% %--------------------------------------------------------------------------
% %recompute velocity
% cv=latlon2uv(num,lat,lon);
% %--------------------------------------------------------------------------
% %remove those containing no valid data
% bool=sum(isfinite(col2mat(cell2col(lat))),1)>0;
% vindex(id,source,num,lat,lon,laterr,lonerr,cv,cverr,gap,drogue,err,t,filled,depth,bool,1);
% %--------------------------------------------------------------------------
% %Sort by decreasing length
% [len,sorter]=sort(cellength(num),'descend');
% vindex(id,source,num,lat,lon,laterr,lonerr,cv,cverr,gap,drogue,err,t,depth,filled,sorter,1);
%--------------------------------------------------------------------------
clear gulfdrifters_secondary
for j=1:length(varnames)
    eval(['make gulfdrifters_secondary ' varnames{j}])
end
length(cell2col(gulfdrifters_primary.lat)) %4500258
length(cell2col(gulfdrifters_secondary.lat)) %4500258
length(find(isnan(cell2col(gulfdrifters_primary.lon)))) %=3960
length(find(isnan(cell2col(gulfdrifters_secondary.lon)))) %=3960
length(find(isinf(cell2col(gulfdrifters_primary.lon)))) %=0
length(find(isinf(cell2col(gulfdrifters_secondary.lon)))) %=28474
length(find(cell2col(gulfdrifters_primary.filled)==1))  %=34358
length(find(cell2col(gulfdrifters_secondary.filled)==1)) %=42241
42241./4460442*100 %= 0.947%
%all that has changed is that more points have been re-interpolated over
%and hence flagged as filled 
matsave gulfdrifters_secondary
%--------------------------------------------------------------------------
%remove more than 10% filled
%fillfrac=cellsum(filled)./cellength(filled);
%vindex(id,source,num,lat,lon,laterr,lonerr,cv,cverr,gap,drogue,err,t,filled,fillfrac<0.1,1);
%\*************************************************************************


function[]=gulfdrifters_tertiary

load gulfdrifters_experiments_primary
names=fieldnames(gulfdrifters_experiments_primary);

%/*************************************************************************
%Final step, mitigating SCULP issues
load gulfdrifters_secondary
use gulfdrifters_secondary
cell2col(cv,source,id,lat,lon,filled);
spd=abs(cv);
%--------------------------------------------------------------------------
%try the median-based robustification of Cleveland
bool=(filled==0);
[mspd,xmid,ymid]=twodstats(lon(bool),lat(bool),abs(cv(bool)),[-99:1/4:-78],[18:1/4:31]);

spdi=interp2(xmid,ymid,mspd,lon,lat,'linear');
spd0=interp2(xmid,ymid,mspd,lon,lat,'nearest');
spdi(~isfinite(spdi))=spd0(~isfinite(spdi));

res=abs(spd-spdi);
[s,xmid,ymid]=twodmed(lon,lat,res,[-99:1/4:-78],[18:1/4:31]);
%[s,xmid,ymid]=twodmed(lon,lat,res,[-99:1/2:-78],[18:1/2:31]);

si=interp2(xmid,ymid,s,lon,lat,'linear');
s0=interp2(xmid,ymid,s,lon,lat,'nearest');
si(~isfinite(si))=s0(~isfinite(si));

eps=abs(frac(res,si));
%delta=squared(1-squared(frac(res,6*si)));
%delta(frac(res,6*si)>1)=0;  %set extreme outliers to zero
%make gulfdrifters delta
% %--------------------------------------------------------------------------
% bool1=eps>10&(source==1|source==2);
% [N,a,b]=blocknum(bool1);a=a(2:2:end);b=b(2:2:end);
% figure,plot(sort(b-a+1))
% bool2=eps>10&~(source==1|source==2);
% [N,a,b]=blocknum(bool2);a=a(2:2:end);b=b(2:2:end);
% hold on, plot(sort(b-a+1))
% %in general the lengths are much shorter for non-SCULP ... the largest 
% %being a thirty-hour gap
%--------------------------------------------------------------------------
%re-do the whole above interpolation and chunking etc  
epscutoff=7.5;
bool=eps>epscutoff&(source==2|source==3)&(isfinite(filled)&filled==0);length(find(bool))  %2427
length(find((source==2|source==3)&(isfinite(filled)&filled==0)))%956986
2415/956986*100 %0.25%
filled(bool)=true;
lon(bool)=inf;
lat(isinf(lon))=inf;
length(find(lon)) %4500258
length(find(isnan(lon))) %3960
col2cell(id,source,lat,lon,filled,cv);
for i=1:length(lat)
    index=find(isfinite(lat{i}),1,'first'):find(isfinite(lat{i}),1,'last');
    [num{i},id{i},source{i},lat{i},lon{i},id{i},cv{i},gap{i},drogue{i},t{i},filled{i},depth{i}]=...
        vindex(num{i},id{i},source{i},lat{i},lon{i},id{i},cv{i},gap{i},drogue{i},t{i},filled{i},depth{i},index,1);
end
cellprune(num,lat,lon,id,typenum,source,cv,gap,drogue,t,filled,depth);
vsize(num,lat,lon,id,typenum,source,cv,gap,drogue,t,filled,depth) %now change
for i=1:length(lat)
    lat{i}=fillbad(lat{i},inf);
    lon{i}=fillbad(lon{i},inf);
end
%--------------------------------------------------------------------------
%remove segments less than 1 day in length
vsize(id,typenum,source,num,lat,lon,cv,gap,drogue,t,filled,depth)
vindex(id,typenum,source,num,lat,lon,cv,gap,drogue,t,filled,depth,(cellmax(num)-cellmin(num))>=1,1);
vsize(id,typenum,source,num,lat,lon,cv,gap,drogue,t,filled,depth)
%--------------------------------------------------------------------------
%recompute velocity
cv=latlon2uv(num,lat,lon);
%length(find(~isfinite(cell2col(cv))))  %3308, good
%--------------------------------------------------------------------------
%remove more than 10% filled
%fillfrac=cellsum(filled)./cellength(filled);  %note!!! 
%vindex(id,num,lat,lon,laterr,lonerr,cv,cverr,gap,drogue,err,t,filled,fillfrac<0.1,1);
%--------------------------------------------------------------------------
%Sort by decreasing length
%[len,sorter]=sort(cellength(num),'descend');
%vindex(id,source,num,lat,lon,laterr,lonerr,cv,cverr,gap,drogue,err,t,filled,depth,sorter,1);
%--------------------------------------------------------------------------
id=cellfirst(id);
source=cellfirst(source);

for i=1:length(names)
    if strcmpi(names{i},'hercules')
        names{i}='Hercules';
    else 
        names{i}=upper(names{i});
    end
end

%sort by start date
[~,sorter]=sort(cellmin(num));
vindex(id,typenum,source,num,lat,lon,cv,t,drogue,gap,depth,filled,sorter,1);
%sort by source
[~,sorter]=sort(source);
vindex(id,typenum,source,num,lat,lon,cv,t,drogue,gap,depth,filled,sorter,1);

%recompute depth
[topo,lato,lono]=readtopo([-99   -80    18    31]);
depth=lat;
parfor j=1:length(lon)
    depth{j}=interp2(lato,lono,-topo',lat{j},lon{j});
end

% %which of these code-type drifters actually has a drogue flag?
% for i=1:length(id)
%     if typenum(i)==2&&anyany(isfinite(drogue{i}))
%         source(i)
%     end
% end
%just dwde

%set drogue status to 'present' for all of the CODE-type drifters + Tube
%that's typenum 2 or 5
for i=1:length(lon)
    if typenum(i)==2||typenum(i)==5
        drogue{i}=1+zeros(size(drogue{i}));
    end
end

%length(find(cell2col(drogue(typenum==2))))
%length(find(cell2col(drogue(typenum==2))==inf))
%length(find(cell2col(drogue(typenum==2))==0))
%ok

clear gulfdrifters
about='For more information, type ''about_gulfdrifters''.';
matsave gulfdrifters about names id source typenum num lat lon cv t drogue gap depth filled
%\*************************************************************************

length(lon)%3768
length(cell2col(lon))%  4505546
length(find(isinf(cell2col(lon))))


% /*************************************************************************
% finding a small number of drifters that give suspiciously large velocities
% use gulfflow_onetwelfth
% use gulfflow_onequarter
% spd=vmean(abs(u+1i*v),3);
% figure,jpcolor(lon,lat,spd) %zoom on on suspicious regions 
% region=axis;
% use gulfdrifters
% sources=ids;
% for i=1:length(sources)
%     sources{i}=source(i)+0*sources{i};
%     ids{i}=drifter_id(i)+0*ids{i};
% end
% vsize(ids,time,lat,lon,sources)
% cell2col(ids,time,lat,lon,sources);
% bool=inregion(region,lat,lon);
% sources(bool)
% ids(bool)
% ids 80253 in SGOM | ids 80125 80216 80224 in SGOM
% ids 10075 10087 18320 in DWDE | ids 10082 in DWDE | ids 10029 in DWDE | 17060* in DWDE | 10075* |10026*
% ids 406 in SCULP2
% ids 101 and 317 in GLAD
% ids 94 39 and 226 in SCULP1
% ids 307* in splash
% ids 63701930* in hgps
% \*************************************************************************

function[]=gulfdrifters_netcdf(netcdfdir)

% 1 latex:6 6
% 2 sculp1:1.5 1.5
% 3 sculp2:1.5 1.5
% 4 gdp:6 6
% 5 hargos:1 1.6136
% 6 aoml:1 1
% 7 sgom:1 1.0083
% 8 ngom:1 3.101
% 9 ocg:1 1.9121
% 10 glad:0.25 0.25
% 11 hercules:0.083333 0.33189
% 12 hgps:1 1.1738
% 13 laser:0.25 0.33095
% 14 dwde:1 1.0007
% 15 splash:0.083333 0.31442


for ii=1:3
%/*************************************************************************
load gulfdrifters
use gulfdrifters

L=cellength(num);
traj_id=[1:length(id)]';
ids=num;
for i=1:length(id)
    ids{i}=i+0*num{i};
end

if ii==2
    bool=(source~=7)&(source~=8)&(source~=14);
    vindex(id,L,traj_id,source,typenum,ids,num,lat,lon,cv,t,drogue,gap,depth,filled,bool,1);
elseif ii==3
    bool=(source~=7)&(source~=8);
    vindex(id,L,traj_id,source,typenum,ids,num,lat,lon,cv,t,drogue,gap,depth,filled,bool,1);
end

cell2col(ids,num,lat,lon,cv,t,drogue,gap,depth,filled);

%remove nans marking tails
bool=isnan(num);%length(find(bool))
vindex(ids,num,lat,lon,cv,t,drogue,gap,depth,filled,~bool,1);

% length(find(~isfinite(num)))
% length(find(~isfinite(lat)))
% length(find(~isfinite(lon)))
% length(find(~isfinite(cv)))
% length(find(~isfinite(t)))
% length(find(~isfinite(drogue)))
% length(find(~isfinite(gap)))
% length(find(~isfinite(depth)))
% length(find(~isfinite(filled)))
vswap(t,drogue,gap,inf,nan);

names_mat=zeros(length(names),max(cellength(names)));
for i=1:size(names,1)
    names_mat(i,1:length(names{i}))=names{i};
end
names_mat=setstr(names_mat);

types_mat=['WOCE     ';'CODE     ';'SVP      ';'FHD      ';'Tube     ';'CARTHE   ';'Microstar';'DORIS    '];

%Following these convention
%https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3
%http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html#trajectory-data
%List of names is at http://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
%see http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html#standard-name

writedir=['/Users/lilly/Desktop/Dropbox/NetCDF/'];

if ii==1
    name_nc=['GulfDriftersAll.nc'];
elseif ii==2
    name_nc=['GulfDriftersOpen.nc'];
elseif ii==3
    name_nc=['GulfDriftersDWDE.nc'];
end
% netCDF
% Define the global attributes
ncid = netcdf.create([writedir name_nc],'NETCDF4');
varid = netcdf.getConstant('GLOBAL');
if ii==1
    netcdf.putAtt(ncid,varid,'title','GulfDriftersAll');
    netcdf.putAtt(ncid,varid,'summary','GulfDriftersAll is a consolidated dataset containing all available surface drifter data for the Gulf of Mexico, subjected to a uniform processing and quality control methodology, and interpolated onto hourly resolution.  It is not for public distribution.')
elseif ii==2
    netcdf.putAtt(ncid,varid,'title','GulfDriftersOpen');
    netcdf.putAtt(ncid,varid,'summary','GulfDriftersOpen is a consolidated dataset containing all publicly-available surface drifter data for the Gulf of Mexico, subjected to a uniform processing and quality control methodology, and interpolated onto hourly resolution.')
elseif ii==3
    netcdf.putAtt(ncid,varid,'title','GulfDriftersDWDE');
    netcdf.putAtt(ncid,varid,'summary','GulfDriftersDWDE is a consolidated dataset containing all publicly-available surface drifter data for the Gulf of Mexico, plus the DWDE dataset that is avaiable for noncommercial use only, subjected to a uniform processing and quality control methodology, and interpolated onto hourly resolution.')
end
netcdf.putAtt(ncid,varid,'product_version','1.1.0');
netcdf.putAtt(ncid,varid,'keywords','Gulf of Mexico, surface currents, surface drifters, Lagrangian data, Bay of Campeche, Loop Current');
netcdf.putAtt(ncid,varid,'Conventions','CF-1.6, ACDD-1.3');
netcdf.putAtt(ncid,varid,'references','Lilly and Perez-Prunius (2021a). A gridded surface current product for the Gulf of Mexico from consolidated drifter measurements.  Earth Syst. Sci. Data.');

netcdf.putAtt(ncid,varid,'id','gulfdrifters');
netcdf.putAtt(ncid,varid,'naming_authority','http://www.jmlilly.net');
netcdf.putAtt(ncid,varid,'history',[datestr(now,1) ' v. 1.1.1 [next update] | ' ...
                                    '04-Jan-2021 v. 1.1.0 updated to include Global Drifter Program data through mid-2020, added drifter type information, incorporated additional drogue flags | ' ...  
                                    '05-May-2020 v. 1.0.0 initial version']);
netcdf.putAtt(ncid,varid,'source','surface drifter');
netcdf.putAtt(ncid,varid,'platform','surface drifter');
netcdf.putAtt(ncid,varid,'instrument','WOCE drifter, CODE drifter, SVP drifter, Far Horizon Drifter (FHD), CARTHE drifter, Microstar drifter, DORIS drifter');

netcdf.putAtt(ncid,varid,'cdm_data_type','trajectory');
netcdf.putAtt(ncid,varid,'featureType','trajectory');
%netcdf.putAtt(ncid,varid,'data_type','trajectory data');

netcdf.putAtt(ncid,varid,'processing_level','Variable levels of processing by original investigators, followed by uniform processing and quality control as described in Lilly and Perez-Brunius (2021a).');
if ii==1
    netcdf.putAtt(ncid,varid,'comment',['Thanks to Paula Garcia Carrillo for help with the NetCDF conversion.  ' ...
    'To convert this dataset to cell array form in Matlab, with the jLab toolbox installed from https://github.com/jonathanlilly/jLab, use ncload(''GulfDriftersAll.nc'')']);
elseif ii==2 
    netcdf.putAtt(ncid,varid,'comment',['Thanks to Paula Garcia Carrillo for help with the NetCDF conversion.  ' ...
    'To convert this dataset to cell array form in Matlab, with jLab installed, use ncload(''GulfDriftersOpen.nc'')']);
elseif ii==3
    netcdf.putAtt(ncid,varid,'comment',['Thanks to Paula Garcia Carrillo for help with the NetCDF conversion.  ' ...
    'To convert this dataset to cell array form in Matlab, with jLab installed, use ncload(''GulfDriftersDWDE.nc'')']);
end

netcdf.putAtt(ncid,varid,'acknowledgment','A complete list of all acknowledgements for the various datasets incorporated herein can be found in Lilly and Perez-Prunius (2021a).');
if ii~=3
    netcdf.putAtt(ncid,varid,'license','Creative Commons Attribution 4.0 International');
elseif ii==3
    netcdf.putAtt(ncid,varid,'license','This dataset is licensed for non-commercial academic use only, and it is prohibited to share it with third parties, as well as to profit or sell products derived from it.');
end
netcdf.putAtt(ncid,varid,'standard_name_vocabulary','CF Standard Name Table');
netcdf.putAtt(ncid,varid,'date_created','05-May-2020');
netcdf.putAtt(ncid,varid,'date_modified',datestr(now,1));
netcdf.putAtt(ncid,varid,'creator_name','Jonathan M. Lilly');
netcdf.putAtt(ncid,varid,'creator_email','jmlilly@psi.edu');
netcdf.putAtt(ncid,varid,'creator_url','http://www.jmlilly.net');
netcdf.putAtt(ncid,varid,'creator_institution','Planetary Science Institute');
netcdf.putAtt(ncid,varid,'institution','Planetary Science Institute');
netcdf.putAtt(ncid,varid,'contributor_name','Paula Perez-Brunius');
netcdf.putAtt(ncid,varid,'contributor_email','brunius@cicese.mx');
netcdf.putAtt(ncid,varid,'project','Lagrangian Data Analysis for the Gulf of Mexico');
netcdf.putAtt(ncid,varid,'publisher_name','Jonathan M. Lilly')
netcdf.putAtt(ncid,varid,'publisher_email','eponym@jmlilly.net');
netcdf.putAtt(ncid,varid,'publisher_url','http://www.jmlilly.net')
netcdf.putAtt(ncid,varid,'geoespatial_lat_min',min(lat));
netcdf.putAtt(ncid,varid,'geoespatial_lat_max',max(lat));
netcdf.putAtt(ncid,varid,'geoespatial_lat_units','degree_north');
netcdf.putAtt(ncid,varid,'geoespatial_lon_min',min(lon));
netcdf.putAtt(ncid,varid,'geoespatial_lon_max',max(lon));
netcdf.putAtt(ncid,varid,'geoespatial_lon_units','degree_east');
netcdf.putAtt(ncid,varid,'time_coverage_start',datestr(min(num)));
netcdf.putAtt(ncid,varid,'time_coverage_end',datestr(max(num)));
netcdf.putAtt(ncid,varid,'time_coverage_resolution','1h');

%netcdf.putAtt(ncid,varid,'site_code','Gulf of Mexico');
%netcdf.putAtt(ncid,varid,'platform_code','drifter');
%netcdf.putAtt(ncid,varid,'area','Gulf of Mexico');

% Define the dimensions
dimiobs   = netcdf.defDim(ncid,'obs',length(num));
dimifloat = netcdf.defDim(ncid,'traj',length(id));

dimiexp  = netcdf.defDim(ncid,'exp',size(names_mat,1));
dimistr  = netcdf.defDim(ncid,'exp_str_length',size(names_mat,2));
dimitypes  = netcdf.defDim(ncid,'type',size(types_mat,1));
dimistr2 = netcdf.defDim(ncid,'type_str_length',size(types_mat,2));

% Define IDs for the dimension variables (lon,lat, depth, time...)
exp_str_ID = netcdf.defVar(ncid,'exp_names','char',[dimiexp dimistr]);
netcdf.putAtt(ncid,exp_str_ID,'long_name','experiment names');

type_str_ID = netcdf.defVar(ncid,'drifter_types','char',[dimitypes dimistr2]);
netcdf.putAtt(ncid,type_str_ID,'long_name','drifter types');

traj_ID = netcdf.defVar(ncid,'id','int',dimifloat);
netcdf.putAtt(ncid,traj_ID,'long_name','trajectory segment ID');
netcdf.putAtt(ncid,traj_ID,'cf_role','trajectory_id');

float_ID = netcdf.defVar(ncid,'drifter_id','int',dimifloat);
netcdf.putAtt(ncid,float_ID,'long_name','drifter ID');

exp_ID = netcdf.defVar(ncid,'source','int',dimifloat);
netcdf.putAtt(ncid,exp_ID,'long_name','number designating originating experiment name within exp_names');

type_ID = netcdf.defVar(ncid,'type','int',dimifloat);
netcdf.putAtt(ncid,type_ID,'long_name','number designating drifter design type within drifter_types');

row_size_ID = netcdf.defVar(ncid,'row_size','int',dimifloat);
netcdf.putAtt(ncid,row_size_ID,'long_name','number of observations for each trajectory segment');
netcdf.putAtt(ncid,row_size_ID,'sample_dimension','obs');

% Define the main variables

ids_ID = netcdf.defVar(ncid,'ids','int',dimiobs);
netcdf.putAtt(ncid,ids_ID,'long_name','trajectory segment IDs for all data points');
netcdf.defVarDeflate(ncid,ids_ID,true,true,9);

time_ID = netcdf.defVar(ncid,'time','double',dimiobs);
netcdf.putAtt(ncid,time_ID,'standard_name','time');
netcdf.putAtt(ncid,time_ID,'long_name','time');
%netcdf.putAtt(ncid,time_ID,'units','Days since 00-Jan-0000 00:00:00');
%tried this: update time units to fix error in Python time conversion, but didn't help
netcdf.putAtt(ncid,time_ID,'units','days since 0000-00-00 00:00:00');
netcdf.putAtt(ncid,time_ID,'axis','T')
netcdf.putAtt(ncid,time_ID,'calendar','standard');
netcdf.defVarDeflate(ncid,time_ID,true,true,9);

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

U_ID = netcdf.defVar(ncid,'u','double',dimiobs);
netcdf.putAtt(ncid,U_ID,'standard_name','eastward_sea_water_velocity');
netcdf.putAtt(ncid,U_ID,'long_name','eastward component of sea surface current velocity');
netcdf.putAtt(ncid,U_ID,'units','m/s');
%netcdf.defVarFill(ncid,U_ID,false,NaN);
%netcdf.putAtt(ncid,U_ID,'processing_level','Data Manually reviewed');
netcdf.defVarDeflate(ncid,U_ID,true,true,9);

V_ID = netcdf.defVar(ncid,'v','double',dimiobs);
netcdf.putAtt(ncid,V_ID,'standard_name','northward_sea_water_velocity');
netcdf.putAtt(ncid,V_ID,'long_name','northward component of sea surface current velocity');
netcdf.putAtt(ncid,V_ID,'units','m/s');
%netcdf.defVarFill(ncid,V_ID,false,NaN);
%netcdf.putAtt(ncid,V_ID,'processing_level','Data Manually reviewed');
netcdf.defVarDeflate(ncid,V_ID,true,true,9);

T_ID = netcdf.defVar(ncid,'temp','double',dimiobs);
netcdf.putAtt(ncid,T_ID,'standard_name','sea_surface_temperature');%jml changed from 'sea_water_temperature'
netcdf.putAtt(ncid,T_ID,'long_name','sea surface temperature');
netcdf.putAtt(ncid,T_ID,'units','degree_Celsius');
netcdf.defVarFill(ncid,T_ID,false,nan);%since there are undefined values
netcdf.defVarDeflate(ncid,T_ID,true,true,9);

drogue_ID = netcdf.defVar(ncid,'drogue','double',dimiobs);
netcdf.putAtt(ncid,drogue_ID,'long_name','drogue on / off flag');
netcdf.putAtt(ncid,drogue_ID,'valid_range','0, 1');
netcdf.putAtt(ncid,drogue_ID,'flag_values','0, 1, NaN');
netcdf.putAtt(ncid,drogue_ID,'flag_meanings','drogue_off drogue_on unknown');
netcdf.defVarFill(ncid,drogue_ID,false,nan);%since there are undefined values
netcdf.defVarDeflate(ncid,drogue_ID,true,true,9);

gap_ID = netcdf.defVar(ncid,'gap','double',dimiobs);
netcdf.putAtt(ncid,gap_ID,'long_name','upstream interpolation gap');
netcdf.putAtt(ncid,gap_ID,'units','hours');
netcdf.defVarFill(ncid,gap_ID,false,nan);
netcdf.defVarDeflate(ncid,gap_ID,true,true,9);

depth_ID = netcdf.defVar(ncid,'depth','double',dimiobs);
netcdf.putAtt(ncid,depth_ID,'long_name','ocean bottom depth from interpolation within one-minute Smith and Sandwell global database version 19.1');
netcdf.putAtt(ncid,depth_ID,'units','kilometers');
netcdf.defVarDeflate(ncid,depth_ID,true,true,9);

filled_ID = netcdf.defVar(ncid,'filled','double',dimiobs);
netcdf.putAtt(ncid,filled_ID,'long_name','filled on / off flag');
netcdf.putAtt(ncid,filled_ID,'valid_range','0, 1');%jml removed b's from this
netcdf.putAtt(ncid,filled_ID,'flag_values','0, 1');%jml removed b's from this
netcdf.putAtt(ncid,filled_ID,'flag_meanings','at_least_one_valid_original_datapoint_within_plus_or_minus_six_hours_and_five_minutes otherwise');
netcdf.defVarDeflate(ncid,filled_ID,true,true,9);

% We are done defining the NetCdf
netcdf.endDef(ncid);
% Then store the dimension variables in
netcdf.putVar(ncid,exp_str_ID,names_mat);
netcdf.putVar(ncid,type_str_ID,types_mat);
netcdf.putVar(ncid,traj_ID,traj_id);
netcdf.putVar(ncid,float_ID,id);
netcdf.putVar(ncid,exp_ID,source);
netcdf.putVar(ncid,type_ID,typenum);
netcdf.putVar(ncid,row_size_ID,L);

% Then store my main variable
netcdf.putVar(ncid,ids_ID,ids);
netcdf.putVar(ncid,time_ID,num);
netcdf.putVar(ncid,Y_ID,lat);
netcdf.putVar(ncid,X_ID,lon);
netcdf.putVar(ncid,U_ID,real(cv)/100);%convert cm/s to m/s
netcdf.putVar(ncid,V_ID,imag(cv)/100);%convert cm/s to m/s
netcdf.putVar(ncid,T_ID,t);
netcdf.putVar(ncid,drogue_ID,drogue);
netcdf.putVar(ncid,gap_ID,gap);
netcdf.putVar(ncid,depth_ID,depth);
netcdf.putVar(ncid,filled_ID,filled);
% We're done, close the netCDF
netcdf.close(ncid)
%\*************************************************************************
end

function[]=gulfdrifters_augmented

%/*************************************************************************
%Interpolating model velocities, aviso, and winds onto drifter positions
readdir='/Volumes/Alfheim/';
netcdfdir='/Users/lilly/Desktop/Dropbox/NetCDF';

%copy GulfDriftersAll to writedir and rename GulfDriftersAugmented.nc
eval(['!cp ' netcdfdir '/GulfDriftersAll.nc ' netcdfdir '/GulfDriftersAllAugmented.nc' ])
filename=[netcdfdir '/GulfDriftersAllAugmented.nc'];

num=ncread(filename,'time');
lat=ncread(filename,'lat');
lon=ncread(filename,'lon');

names={'cmems_u','cmems_v','cmems_h','hycom_u','hycom_v','hycom_h','nemo_u','nemo_v','nemo_h','roms_u','roms_v','roms_h'};
longnames={'CMEMS altimetry','CMEMS altimetry','CMEMS altimetry' ...
    'HYCOM model output',  'HYCOM model output',  'HYCOM model output' ...
    'NEMO model output',  'NEMO model output',  'NEMO model output' ...
    'ROMS model output',  'ROMS model output',  'ROMS model output'};

ncid = netcdf.open(filename,'write');
varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid,varid,'title','GulfDriftersAllAugmented');

clear varid
for i=1:length(names)
    varid{i} = netcdf.defVar(ncid,names{i},'double',0); %The obs dimensions is dimension zero
    switch names{i}(end)
        case 'u'
            netcdf.putAtt(ncid,varid{i},'long_name',['eastward component of sea surface current velocity in ' longnames{i}]);
            netcdf.putAtt(ncid,varid{i},'units','m/s');
        case 'v'
            netcdf.putAtt(ncid,varid{i},'long_name',['northward component of sea surface current velocity in ' longnames{i}]);
            netcdf.putAtt(ncid,varid{i},'units','m/s');
        case 'h'
            netcdf.putAtt(ncid,varid{i},'long_name',['sea surface height in ' longnames{i}]);
            netcdf.putAtt(ncid,varid{i},'units','m');
    end
    netcdf.defVarFill(ncid,varid{i},false,NaN);         %since there are undefined values
    netcdf.defVarDeflate(ncid,varid{i},true,true,9);
end

% We are done defining the NetCdf
netcdf.endDef(ncid);


tic;[cmems_u,cmems_v,cmems_h]=ncinterp([readdir 'gom_aviso_madt.nc'],num,lat,lon,'u','v','sla');toc

%Did I not finish this?
%Iamhere

datestr(maxmax(num))
numb=datenum(2013,1,1);
numa=datenum(2000,1,1);
bool=num>=numb;
num(bool)=num(bool)+(numa-numb);
datestr(maxmax(num))

tic;[hycom_u,hycom_v,hycom_h]=ncinterp([readdir 'hycom_surface.nc'],num,lat,lon,'u','v','ssh');toc
tic;[nemo_u,nemo_v,nemo_h]=ncinterp([readdir 'nemo_surface.nc'],num,lat,lon,'u','v','ssh');toc
tic;[roms_u,roms_v,roms_h]=ncinterp([readdir 'roms_surface.nc'],num,lat,lon,'u','v','ssh');toc


for i=1:length(names),eval(['maxmax(' names{i} ')']),end

%need to divide by 100 for all except the cmems variables to convert cm to m
for i = 1:length(names)
    if i>=4
        eval([ names{i} '=' names{i} '/100;'])
    end
end
for i = 1:length(names)
    evalstr=[ 'netcdf.putVar(ncid,' int2str(varid{i}) ',' names{i} ');']
    eval(evalstr)
end  
netcdf.close(ncid);
%--------------------------------------------------------------------------
%winds ... for future reference, using my old format 
%num=ncread(filename,'time');
%
% tic;[merra_uwnd,merra_vwnd]=ncinterp('gom_merra_winds.nc',num,lat,lon,'uwnd','vwnd');toc
% ncwrite(filename,'merra_uwnd',merra_uwnd,1);
% ncwrite(filename,'merra_vwnd',merra_vwnd,1);
% 
% tic;[merra_uflx,merra_vflx]=ncinterp('gom_merra_stress.nc',num,lat,lon,'uflx','vflx');toc
% ncwrite(filename,'merra_uflx',merra_uflx,1);
% ncwrite(filename,'merra_vflx',merra_vflx,1);
% 
% tic;[ccmp_uwnd,ccmp_vwnd]=ncinterp('gom_ccmp_winds.nc',num,lat,lon,'uwnd','vwnd');toc
% ncwrite(filename,'ccmp_uwnd',ccmp_uwnd,1);
% ncwrite(filename,'ccmp_vwnd',ccmp_vwnd,1);
%\*************************************************************************
