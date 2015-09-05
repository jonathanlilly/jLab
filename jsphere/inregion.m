function[bool]=inregion(varargin)
% INREGION  Tests whether lat/lon points lie within a specified box.
%
%   BOOL = INREGION(REGION,LAT,LON), where LAT and LON are arrays having
%   the same size, returns an array BOOL which is true (=1) for all
%   LAT/LON pairs which lie inside REGION, and false (=0) otherwise.
%
%   REGION is an array with the format [WEST EAST SOUTH NORTH]. Longitudes 
%   may either be specified on the interval [-180, 180] or on [0, 360].
%
%   The region may overlap the prime meridian (LON=0) or the dateline
%   (LON=180).  Region boundaries are interpreted to exclude the poles. 
%   ______________________________________________________________________
%
%   Cell array lat/lon input 
%
%   INREGION also works if LAT and LON are cell arrays. 
%
%   In this case, BOOL is also a cell array, with the Kth element of BOOL, 
%   BOOL{K}, being a boolean array of the same size as LAT{K} and LON{K}.
%   ______________________________________________________________________
%
%   Multiple boxes 
%
%   INREGION can detect if points are within any of several boxes.  To do
%   so, just make REGION a cell array, REGION={R1,R2,...,RN} where the RN
%   are numeric arrays, each with the format [WEST EAST SOUTH NORTH].
%  
%   In this way one can build up a more general region from the union of
%   several rectangular regions.  
%   ______________________________________________________________________
%
%   See also REGIONPLOT, FLOATREGION, INELLIPSE.
%
%   'inregion --t' runs a test.
%
%   Usage:  bool=inregion(region,lat,lon);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2014 J.M. Lilly --- type 'help jlab_license' for details        

if strcmpi(varargin{1},'--t')
     inregion_test,return
end

region=varargin{1};
lat=varargin{2};
lon=varargin{3};


if iscell(lat)
    for i=1:length(lat)
        bool{i,1}=inregion_one_loop(region,lat{i},lon{i});
    end
else
   bool=inregion_one_loop(region,lat,lon);
end

function[bool]=inregion_one_loop(region,lat,lon)

if ~iscell(region)
    bool=inregion_one(region,lat,lon);
elseif iscell(region)
    for i=1:length(region)
        %Just choosing an unused dimension
        bool(:,:,:,:,i)=inregion_one(region{i},lat,lon);
    end
    bool=any(bool,5);
end
    
    
    
function[bool]=inregion_one(region,lat,lon)

west=region(1);
east=region(2);
south=region(3);
north=region(4);

west=deg180(west);
east=deg180(east);
lon=deg180(lon);

boollat=(lat>=south)  & (lat<=north);
if west<east
    %Dateline is not inside region
    boollon=(lon>=west)  & (lon<=east);
elseif east<west 
    boollon=(lon>=west) | (lon<=east);
end

bool=boollat&boollon;

function[]=inregion_test
lat=[58 76];
lon=[-52.5 -52.5];
region=[-63 -41 52 65];
ans1=inregion(region,lat,lon);
reporttest('INREGION Labrador Sea', aresame(ans1,[1 0]))

lat=[58 76];
lon=360+[-52.5 -52.5];
region=[-63 -41 52 65];
ans1=inregion(region,lat,lon);
reporttest('INREGION Labrador Sea, differing longitude conventions', aresame(ans1,[1 0]))

lat=[58 58];
lon=[170 150];
region=[160 -165 52 65];
ans1=inregion(region,lat,lon);
reporttest('INREGION enclosing dateline', aresame(ans1,[1 0]))

load ebasnfloats
use ebasnfloats
region=[-30 -21 24 35];
bool1=inregion(region,lat,lon);
region={[-30 -25 24 35],[-25 -21 24 35]};
bool2=inregion(region,lat,lon);
reporttest('INREGION multiple regions', aresame(cell2col(bool1),cell2col(bool2)))

