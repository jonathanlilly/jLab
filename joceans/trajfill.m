function[lat,lon,filled]=trajfill(varargin)
%TRAJFILL  Fills float or drifter trajectories with linear interpolation.
%
%   [LAT,LON]=TRAJFILL(LAT,LON) fills any gaps in the input vectors LAT
%   and LON, as indicated by INF values, using linear interpolation.
%
%   [LAT,LON,FILLED]=TRAJFILL(LAT,LON) also returns a boolean array 
%   FILLED of the same size as LAT and LON.  FILLED is true if that data 
%   point has been filled, and false otherwise. 
% 
%   TRAJFILL also works if the input LAT and LON are cell arrays of 
%   numeric arrays, as in the dataset FLOATS.MAT.  In this case the 
%   output fields are all cell arrays of the same size as the input.
%
%   Specifically TRAJFILL works by converting latitude and longitude
%   values to 3-vectors representing a position on the surface of the
%   earth, interpolating, and then converting back. 
%
%   Usage: [lat,lon]=trajfill(lat,lon);
%          [lat,lon,filled]=trajfill(lat,lon);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2015 J.M. Lilly --- type 'help jlab_license' for details
 

latin=varargin{1};
lonin=varargin{2};
if ~iscell(latin)
    [lat,lon,filled]=trajfill_one(latin,lonin);
    disp(['TRAJFILL filled ' int2str(length(find(filled))) ' points out of ' int2str(length(filled)) '.'])
else
    lat=latin;
    lon=lonin;
    filled=latin;
    for i=1:length(lat)
        [lat{i},lon{i},filled{i}]=trajfill_one(latin{i},lonin{i});
    end
    disp(['TRAJFILL filled ' int2str(length(find(cell2col(filled)))) ' points out of ' int2str(length(cell2col(filled))) '.'])
end


function[lat,lon,filled]=trajfill_one(latin,lonin)
[x,y,z]=latlon2xyz(latin,lonin);

filled=false(size(x));
filled(isinf(x)|isinf(y)|isinf(z))=true;
%vswap(x,y,z,nan,inf);

x=fillbad(x,inf);y=fillbad(y,inf);z=fillbad(z,inf);
alpha=sqrt(x.^2+y.^2+z.^2)./radearth;
x=x./alpha;y=y./alpha;z=z./alpha;
[lat,lon]=xyz2latlon(x,y,z);



