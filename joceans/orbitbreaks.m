function[varargout]=orbitbreaks(varargin)
%ORBITBREAKS  Separate orbit into passes based on turning points.
%
%   ORBITBREAKS is a level function called by ABOUT_ALONGTRACK.
%
%   [LAT,LON]=ORBITBREAKS(LAT,LON), where LAT and LON are column vectors
%   representing a satellite orbital groundtrack, is used to separate 
%   different passes based on the turning points in latitude.
%
%   NANs are inserted into both LAT and LON wherever LAT has a turning
%   point.  Note that no data is removed.
%
%   All latitudes and lontidues are in degrees.  The output longitude has 
%   range [-180, 180].
%
%   Calling COL2MAT(LAT,LON) will then reshape LAT and LON into matrices, 
%   with each orbital pass in its own column.
%   
%   [LAT,LON,X1,X2,... XN]=ORBITBREAKS(LAT,LON,X1,X2,... XN) will insert
%   NANs into the other arguments in the same locations as for LAT and LON.
%
%   The input arguments should contain no NANs.  For missing data, use INFs
%   instead.
%
%   ORBITBREAKS(LAT,LON, ...) with no output arguments overwrites the 
%   original variables.
%
%   See also COLBREAKS, COL2MAT.
%
%   'orbitbreak --t' runs a test.
%
%   Usage: [lat,lon]=orbitbreaks(lat,lon);
%          [lat,lon,x1,x2,x3]=orbitbreaks(lat,lon,x1,x2,x3);
%          orbitbreaks(lat,lon);
%          orbitbreaks(lat,lon,x1,x2,x3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2015 J.M. Lilly --- type 'help jlab_license' for details
 
%This functionality is removed: 
%Also, INFs are inserted into LON whenever LON jumps by more 
%than 180 degrees, to be able to PLOT(LON,LAT) without spurious lines.

if strcmpi(varargin{1}, '--t')
    orbitbreaks_test,turningpoint_test,return
end
 
if nargin<2
    error('At least LAT and LON must be given.')
end

for i=1:nargin
    if ~iscolumn(varargin{i})
        error('All input arguments must be column vectors.')
    end
end

x=varargin{1};
index=find(isfinite(x));
if ~isempty(index)
    if maxmax(abs(x(index)))>90
         error('Looks like you have LAT and LON arguments swapped.')
    end 
end
mat=zeros(size(varargin{1},1),nargin);

for i=1:nargin
   mat(:,i+1)=varargin{i};
end
if anyany(isnan(mat))
    error('Input arguments should contain no NaNs --- use INFs for missing data.')
end

mat=orbitbreaks_turningpoints(mat);
%mat=orbitbreaks_lonjumps(mat);

for i=1:nargin
   varargout{i}=mat(:,i+1);
end

if nargout==0
   eval(to_overwrite(nargin));
end


function[mat]=orbitbreaks_turningpoints(mat)
lat=vswap(mat(:,2),inf,nan);
bool=turningpoint(lat);
bool([1 end],1)=0;
bool=vshift(bool,-1,1);

mat(:,1)=cumsum(bool);
mat=colbreaks(mat);


function[mat]=orbitbreaks_lonjumps(mat)

mat=vswap(mat,nan,-inf);
mat(:,3)=deg180(mat(:,3));

lon=vswap(mat(:,3),-inf,nan);
bool=abs(lon-vshift(lon,-1,1))>90;
mat(:,1)=cumsum(bool);
mat=colbreaks(mat);
vswap(mat,nan,inf);
vswap(mat,-inf,nan);

%Remove last row
mat=mat(1:end-1,:);
function[]=orbitbreaks_test


tol=1e-6;
x1=[ 4 5 6 7     6 5     6    ]';
y1=[ 4 5 6 7 nan 6 5 nan 6 nan]';

[y1a,y1b]=orbitbreaks(x1,x1);

reporttest('ORBITBREAKS simple',aresame([y1 y1],[y1a y1b],tol))

x2=[ 1 1 6 -10     2 2     1    ]';
y2=[ 1 1 6 -10 nan 2 2 nan 1 nan]';
x3=x1;
y3=y1;
[y1o,y2o,y3o]=orbitbreaks(x1,x2,x3);

reporttest('ORBITBREAKS multiple',aresame([y1o y2o y3o],[y1 y2 y3],tol))

orbitbreaks(x1,x2,x3);
 
reporttest('ORBITBREAKS multiple overwriting',aresame([x1 x2 x3],[y1 y2 y3],tol))

% x1=[ 4     5    6          7      8    9       10   11  12]';
% y1=[ 4     5    6    inf   7      8    9  inf  10   11  12 nan]';
% x2=[ 100  110  120        -170  -160 -150     100  110 120  ]';
% y2=[ 100  110  120   inf  -170  -160 -150 inf 100  110 120 nan]';
% x3=x1;
% y3=y1;
% [y1o,y2o,y3o]=orbitbreaks(x1,x2,x3);
% 
% reporttest('ORBITBREAKS multiple with lon breaks',aresame([y1o y2o y3o],[y1 y2 y3],tol))


function[bool]=turningpoint(x)
%TURNINGPOINT  True for turning points, i.e. local extrema, along rows.
%
%   BOOL=TURNINGPOINT(X) where X is a column vector returns BOOL, a
%   vector of the same size as X containing zeros and ones.  BOOL
%   is equal to one if the corresponding element of X is a local
%   maximum or minimum, and zero otherwise.
%
%   X may also be an array of any dimensionality, in which case 
%   local maxima and minima are found with respect to differentiation 
%   along the first dimension, i.e. along rows.
%
%   Note that if any two adjacent points are identical, TURNINGPOINT
%   adds a very small amount of numerical noise to one of the points
%   in order to make a choise about which is larger.
%
%   See also CROSSINGS.
%
%   'turningpoint --t' runs a test.
%
%   Usage: bool=turningpoint(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2015 J.M. Lilly --- type 'help jlab_license' for details
 
%if strcmpi(x, '--t')
%    turningpoint_test,return
%end

%Add very small amout of noise to keep identical points from occurring
index=find(x==vshift(x,-1,1));
if ~isempty(index)
    x(index)=x(index)+randn(size(x(index)))*1e-10;
end

boolmin=~(x>vshift(x,1,1)|x>vshift(x,-1,1));
boolmax=x>vshift(x,1,1)&x>vshift(x,-1,1);
bool=boolmin|boolmax;
bool(1,:)=0;
bool(end,:)=0;

function[]=turningpoint_test

x=[ 4 5 6 7 6 5 6]';
b=[ 0 0 0 1 0 1 0]';

reporttest('TURNINGPOINT',aresame(b,turningpoint(x)))

x= [ 4 5 6 7 7 6 5 6]';
b1=[ 0 0 0 1 0 0 1 0]';
b2=[ 0 0 0 0 1 0 1 0]';

y1=turningpoint(x);
reporttest('TURNINGPOINT with repeated entry',aresame(y1,b1) || aresame(y1,b2))

