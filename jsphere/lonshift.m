function[varargout]=lonshift(varargin)
%LONSHIFT  Shifts longitude origin for plotting purposes.
%
%   LON=LONSHIFT(LONO,LON) shifts the array LON such that its first
%   point its just greater than LONO. This is used for controlling
%   where the longitude break occurs on a map of the earth.
%
%   LON may be given either on the interval [-180,+180] or on [0,360],
%   and LONO may take on any value, e.g. -340, 20, 380, etc.
%   
%   [LON,X1,X2,...XN]=LONSHIFT(LONO,LON,X1,X2,...XN) also shifts the
%   N matrices X1...XN.  These matrices must be oriented such that 
%   their number of columns is the same as the length of LON. 
%
%   Upon output, LON will increase from a minimum value of LONO.   
%
%   'lonshift --t' runs a test.
%
%   Usage:  lon=lonshift(lono,lon);
%           [lon,x1,x2,x3]=lonshift(lono,lon,x1,x2,x3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2013 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    lonshift_test,return
end

clono=varargin{1};
lon=varargin{2};
na=nargin-2;
varargin=varargin(3:end);

lon=deg360(lon);
clon=deg360(clono);

if ~isvector(lon)
    error('LON must be a one-dimensional array.')
else
    if size(lon,2)>size(lon,1)
        brow=1;
    else
        brow=0;
        lon=lon';
    end
end

clon=clon-1e-15;

if clon<minmin(lon)||clon>maxmax(lon)
    [temp,index]=min(lon);
else
    index=find(vshift(lon,1,2)>clon&lon<=clon,1);
end

lon=vshift(lon,index,2);

if ~brow
    lon=lon';
end

lon=lon-(clon-clono);

varargout{1}=frac(360,2*pi)*unwrap(lon*frac(2*pi,360));

for i=1:na;
   varargout{i+1}=vshift(varargin{i},index,2);
end


function[]=lonshift_test

lon=(-179.5:1:179.5)';
lat=(-59.5:1:59.5)';
[lonmat,latmat]=meshgrid(lon,lat);

[lon1,lonmat1,latmat1]=lonshift(30,deg360(lon),lonmat,latmat);
[lon2,lonmat2,latmat2]=lonshift(180,lon1,lonmat1,latmat1);

tol=1e-10;
reporttest('LONSHIFT reversibility',aresame(deg360(lon),deg360(lon2),tol)&&aresame(latmat,latmat2)&&aresame(deg360(lonmat),deg360(lonmat2)))
