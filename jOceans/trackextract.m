function[varargout]=trackextract(varargin)
%TRACKEXTRACT  Extracts alongtrack altimetry segments within given region.
%
%   TRACKEXTRACT is designed to extract alongtrack altimeter data from 
%   the file TPJAOS.MAT within a specified latigude and longitude region. 
%   See ABOUT_TPJAOS for more details on that dataset.
% 
%   [LAT,LON]=TRACKEXTRACT(LAT,LON,REGION) where the input LAT and LON are
%   matrices of latitude and longitude, outputs the sub-matrices LAT and 
%   LON that are contained entirely within REGION.
% 
%   REGION is an array with the format [WEST EAST SOUTH NORTH]. Longitudes 
%   may either be specified on the interval [-180, 180] or on [0, 360].
%
%   Locations of the input LAT and LON within REGION are found on a column-
%   by-column basis.  Those columns having data within the region are
%   included in the output LAT and LON, while the others are excluded.  
%
%   The columns of the output are padded with NaNs at their trailing end
%   to make them all the same length.  Their length is such that the column
%   having the most number of data points in REGION will contain no NaNs.
%
%   TRACKEXTRACT with no output arguments overwrites the original named
%   output variables. 
%   __________________________________________________________________
%
%   Additional input arguments
%
%   [LAT,LON,SSH]=TRACKEXTRACT(LAT,LON,SSH,REGION) where SSH is a 3D array
%   with its first two dimensions matching the input LAT and LON, also
%   extracts the corresponding portions of this array.  
%
%   In this case, the output SSH will have the same number of rows and 
%   columns as the output LAT and LON, and the same number of "pages" or 
%   elements in the third dimension as the input SSH.    
%
%   [LAT,LON,X1,X2,...,XN]=TRACKEXTRACT(LAT,LON,X1,X2,...,XN,REGION) with
%   multiple input arguments also works.
%   __________________________________________________________________
%
%   Composite regions
%
%   REGION is passed directly to the function INREGION.  Thus, REGION can
%   also be a cell array built up from several region boxes.  See INREGION
%   for more details on this option.
%   __________________________________________________________________
%
%   Specifying length cutoff
%
%   [LAT,LON,...]=TRACKEXTRACT(LAT,LON,,...,REGION,LMIN) only returns those 
%   track segments containing LMIN or more points. The default is LMIN=0.
%   __________________________________________________________________
%
%   See also INREGION, REGIONPLOT, TRAJEXTRACT.
%
%   'trackextract --f' generates a sample figure for the Labrador Sea.  
%
%   Usage: [lat,lon]=trackextract(lat,lon,region);
%          [lat,lon,ssh]=trackextract(lat,lon,ssh,region);
%          [lat,lon,ssh,mss,atd]=trackextract(lat,lon,ssh,mss,atd,region,200);
%          trackextract(lat,lon,ssh,region,200);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015--2016 J.M. Lilly --- type 'help jlab_license' for details

if strcmpi(varargin{1}, '--f')
    if  exist('tpjaos.mat')~=2
        disp('Sorry, this demo figure requires that TPJAOS.MAT is on the search path.')
        disp('You can get it from http://www.jmlilly.net/jmldata.html.')
    else
        disp('Generating sample figure, this may take a moment...')
        type makefigs_trackextract
        makefigs_trackextract;
    end
    return
end

if length(varargin{end})==1
    Lmin=varargin{end};
    varargin=varargin(1:end-1);
else
    Lmin=0;
end
region=varargin{end};
varargin=varargin(1:end-1);

lato=varargin{1};
lono=varargin{2};

bool=inregion(region,lato,lono);
if Lmin~=0
    for i=1:size(bool,2);
        if sum(bool(:,i))<Lmin
            bool(:,i)=false;
        end
    end
end


[lat,lon]=vzeros(max(sum(bool)),length(sum(bool)~=0),nan);
for i=1:size(lat,2)
    lat(1:length(find(bool(:,i))),i)=lato(bool(:,i),i);
    lon(1:length(find(bool(:,i))),i)=lono(bool(:,i),i);
end
varargout{1}=lat;
varargout{2}=lon;

for j=3:length(varargin)
    temp=vzeros(max(sum(bool)),length(sum(bool)~=0),size(varargin{j},3),nan);
    for i=1:size(lat,2)
        temp(1:length(find(bool(:,i))),i,:)=varargin{j}(bool(:,i),i,:);
    end
    varargout{j}=temp;
end
varargout{end+1}=bool;

index=find(sum(~isnan(lat))>0);
for i=1:length(varargin)
    varargout{i}=varargout{i}(:,index,:,:,:);
end

eval(to_overwrite(length(varargin)));




