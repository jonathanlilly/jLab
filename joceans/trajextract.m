function[varargout]=trajextract(varargin)
%TRAJEXTRACT  Extracts Lagrangian trajectory segments within given region.
%
%   [LAT,LON]=TRAJEXTRACT(LAT,LON,REGION) where the input LAT and LON are 
%   cell arrays of columns vectors of latitude and longitude, outputs new
%   cell arrays with LAT and LON of the segments contained within REGION.
%
%   This format for LAT and LON is used by two Lagrangian datasets
%   distributed with JLAB, FLOATS.MAT and DRIFTERS.MAT.  Type 'help floats'
%   or 'help drifters' for more details.
%
%   REGION is an array with the format [WEST EAST SOUTH NORTH]. Longitudes 
%   may either be specified on the interval [-180, 180] or on [0, 360].
%
%   The region may overlap the prime meridian (LON=0) or the dateline
%   (LON=180).  Region boundaries are interpreted to exclude the poles. 
%
%   If an input LAT and LON record contains more than one segment passing 
%   through REGION, these are split into separate entries in the output.
%
%   TRAJEXTRACT with no output arguments overwrites the original named
%   output variables. 
%   __________________________________________________________________
%
%   Multiple input arguments
%
%   [LAT,LON,X1,X2,...,XN]=TRAJEXTRACT(LAT,LON,X1,X2,...,XN,REGION) also
%   returns the portions of the cell array variables X1,X2,...,XN within
%   the region.  The input XN are the same size as the input LAT and LON.
%
%   Any of the input XN may also be numerical arrays of the same length as 
%   LAT and LON, rather than cell arrays.  The corresponding output 
%   variable will then also be a numerical array.  An example of such a
%   field is the identification number used in FLOATS.MAT and DRIFTERS.MAT.  
%   __________________________________________________________________
%
%   Specifying length cutoff
%
%   [LAT,LON,...]=TRAJEXTRACT(LAT,LON,,...,REGION,LMIN) only returns those 
%   segments containing LMIN or more points. The default is LMIN=2.
%   __________________________________________________________________
%
%   REGION is passed directly to the function INREGION.  Thus, REGION can
%   also be a cell array built up from several region boxes.  See INREGION
%   for more details on this option.
%
%
%   See also INREGION, REGIONPLOT, TRACKEXTRACT.
%
%   'trajextract --t' runs a test.
%   'trajextract --f' generates a sample figure.
%
%   Usage: [lat,lon]=trajextract(lat,lon,region);
%          [lat,lon,num,id]=trajextract(lat,lon,num,id,region);
%          [lat,lon,num,id]=trajextract(lat,lon,num,id,region,200);
%          trajextract(lat,lon,num,id,region,200);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2016 J.M. Lilly --- type 'help jlab_license' for details

if strcmpi(varargin{1}, '--f')
    type makefigs_trajextract
    makefigs_trajextract;
    return
end
if strcmpi(varargin{1}, '--t')
    trajextract_test,return
end

if length(varargin{end})==1
    Lmin=varargin{end};
    varargin=varargin(1:end-1);
else
    Lmin=2;
end
region=varargin{end};
varargin=varargin(1:end-1);

lat=varargin{1};
lon=varargin{2};

%/**************
%Put numerical array input into cell arrays
bid=false(size(varargin));
for i=3:length(bid) %Lat and lon are not allowed to be arrays
    if ~iscell(varargin{i})
        bid(i)=true;
        varargin{i}=celladd(varargin{i},cellmult(0,varargin{1})); 
    end
end
%\**************

bool=inregion(region,lat,lon);

%Keep first in to last out
for i=1:length(bool)
    a=find(bool{i},1,'first');
    b=find(bool{i},1,'last');
    index{i,1}=(a:b)';
    
    if length(bool{i})>1     
       bool{i}=bool{i}+0;  %Convert logical to numeric
       %Mark locations where the next point is out of the region
       bool{i}(find(bool{i}(1:end-1)&~bool{i}(2:end))+1)=nan;
    end
end

[lat,lon,bool]=cellindex(lat,lon,bool,index);

varargout{1}=lat;
varargout{2}=lon;

for i=3:length(varargin)
    varargout{i}=cellindex(varargin{i},index);
end

varargout{end+1}=bool;

%Strip short trajectories, to make next sorting quicker
L=cellength(varargout{1});
bool=bool(L>Lmin);
for i=1:length(varargin)
    varargout{i}=varargout{i}(L>Lmin);
end

%L=cellength(varargout{1});
%figure,plot(L)

%Reject individual trajectory segments outside of region
for i=1:length(varargin)
    [temp,varargout{i}]=cell2col(bool,varargout{i});
    temp(isinf(temp))=nan;  %Because cell2col replaces pre-existing nans
    [temp,varargout{i}]=col2mat(temp,varargout{i});
    
    %Remove segments that lie outside
    [temp,varargout{i}]=vindex(temp,varargout{i},temp(1,:)==1,2);
    [temp,varargout{i}]=mat2col(temp,varargout{i});
    [temp,varargout{i}]=col2cell(temp,varargout{i});
end

%Strip short trajectories again
L=cellength(varargout{1});
for i=1:length(varargin)
    varargout{i}=varargout{i}(L>Lmin);
end

%Return numerical array input back to numerical arrays
for i=1:length(bid)
    if bid(i)
        varargout{i}=cellfirst(varargout{i});
    end
end

eval(to_overwrite(length(varargin)));

function[]=trajextract_test
load ebasnfloats
use ebasnfloats

region=[-30 -21 24 35];
[lat1,lon1]=trajextract(lat,lon,region,200);
reporttest('TRAJEXTRACT length cutoff',min(cellength(lat1))>=200)

use ebasnfloats
lon=celladd(-155,lon);

region=[deg180(-30-155) deg180(-21-155) 24 35];
[lat1,lon1]=trajextract(lat,lon,region);
lon1=celladd(155,lon1);

use ebasnfloats
region=[-30 -21 24 35];
[lat2,lon2]=trajextract(lat,lon,region);

reporttest('TRAJEXTRACT crossing dateline',aresame(cell2col(lat1),cell2col(lat2),1e-8)&&aresame(cell2col(lon1),cell2col(lon2),1e-8))

use ebasnfloats
region=[-30 -21 24 35];
[lat1,lon1,lon2,lat2]=trajextract(lat,lon,lon,lat,region);

reporttest('TRAJEXTRACT additional input arguments',aresame(cell2col(lat1),cell2col(lat2),1e-8)&&aresame(cell2col(lon1),cell2col(lon2),1e-8))
 

 
