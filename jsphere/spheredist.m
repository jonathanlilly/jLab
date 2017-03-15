function[d]=spheredist(varargin)
%SPHEREDIST  Computes great circle distances on a sphere.
%
%   D=SPHEREDIST(LAT1,LON1,LAT2,LON2) computes arc length in kilometers 
%   along a great circle on the earth from the point with latitude and
%   longitude coordinates LAT1 and LON1 to the point LAT2 and LON2. 
%  
%   All input arguments should be arrays of the same size, or else either
%   LAT1,LON1 or LAT2,LON2 should be scalars.
%   __________________________________________________________________
%
%   Cumulative distance
%
%   With two input arguments, SPHEREDIST instead returns the cumulative
%   distance from the first point in each column.
%
%   D=SPHEREDIST(LAT,LON) where LAT and LON are 2D arrays of the same size,
%   computes the cumulative great circle distance along each row.
%
%   The great circle distance between adjacent LAT/LON points are cumputed,
%   and these are then cumsummed to generate the cumulative distance D.
%
%   In this case, SPHEREDIST can accommodate NaNs at the top of each row,
%   and will begin with the first non-NaN value.  Any interior NaNs however
%   will cause all subsequent output values in that row to also be NaNs.
%   __________________________________________________________________
%
%   The Earth is approximated as a sphere of radius RADEARTH.
%
%   D=SPHEREDIST(...,R) instead uses a sphere of radius R, in kilometers.
%
%   'spheredist --t' runs a test.
%
%   Usage:  d=spheredist(lat,lon,lato,lono); 
%           d=spheredist(lat,lon); 
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2015 J.M. Lilly --- type 'help jlab_license' for details    

if strcmpi(varargin{1}, '--t')
  spheredist_test,return
end

if ischar(varargin{end})
    flag=varargin{end};
    varargin=varargin(1:end-1);
else
    flag='haversine';
end

if length(varargin)==3
   R=varargin{3};
   varargin=varargin(1:end-1);
elseif length(varargin)==5
   R=varargin{5};
else
   R=radearth;
end

if length(varargin)==2
    if ~isvector(varargin{1})||~isvector(varargin{2})
        if ~(aresame(size(varargin{1}),size(varargin{2}))&&(ndims(varargin{1})==2))
            error('Sorry, with two or three input arguments, LAT and LON must be 1D arrays or 2D matrices.')
        end
    end
    lat1=varargin{1}(1:end-1,:);
    lon1=varargin{2}(1:end-1,:);
    lat2=varargin{1}(2:end,:);
    lon2=varargin{2}(2:end,:);
    %d=spheredist(lat1,lon1,lat2,lon2,R);
else
    lat1=varargin{1}; 
    lon1=varargin{2};
    lat2=varargin{3};
    lon2=varargin{4};
end

if ~aresame(size(lat1),size(lon1))
   error('LAT1 and LON1 must be the same size.')
end
if ~aresame(size(lat2),size(lon2))
   error('LAT2 and LON2 must be the same size.')
end
%vsize(lat1,lon1,lat2,lon2) 
if (~aresame(size(lat1),size(lat2)) && (~isscalar(lat1) && ~isscalar(lat2)))
   error('LAT1 and LON1 must be the same size as LAT1 and LON1, or one pair must be scalars.')
end
if anyany(isfinite(lat1)&abs(lat1)>90)||anyany(isfinite(lat2)&abs(lat2)>90)
    error('Sorry, latitudes must be between -90 and +90')
end
%if anyany(abs(deg180(lon1))>180)||anyany(abs(deg180(lon2))>180)
%    error('Sorry, longitudes must be between -180 and +180, or 0 and 360.')
%end

d=spheredist_inner(lat1,lon1,lat2,lon2,R,flag);

if length(varargin)==2
    if ~anyany(isnan(d(1,:)))
        d=[zeros(size(d(1,:)));cumsum(d)];
    else
        %Correct for leading NANs
        d0=d;
        d=vzeros(size(d0,1)+1,size(d0,2));
        for i=1:size(d,2)
            m=find(~isnan(d0(:,i)),1,'first');
            d(m+1:end,i)=cumsum(d0(m:end,i));
            d(1:m-1,i)=nan;
        end
    end
end

function[d]=spheredist_inner(lat1,lon1,lat2,lon2,R,flag)
%Distance for each column

if (~aresame(size(lat1),size(lat2)))
    if isscalar(lat1)
        lat1=lat1+zeros(size(lat2));
        lon1=lon1+zeros(size(lat2));
    elseif isscalar(lat2)
        lat2=lat2+zeros(size(lat1));
        lon2=lon2+zeros(size(lat1));
    end
end
        
d=[];
if ~isempty(lat1)&&~isempty(lat2)
    %This saves some unnecessary computations if we have lots of NANs
    d=nan.*zeros(size(lat1));
    %index=find(~isnan(lat1.*lat2.*lon1.*lon2));
    index=find(~isnan(lat1));
    %length(index)
    if ~isempty(index)
        d(index)=spheredist_one(lat1(index),lon1(index),lat2(index),lon2(index),R,flag);
    end
end

bool=isnan(lat1)|isnan(lat2)|isnan(lon1)|isnan(lon2);
d(bool)=nan;


function[d]=spheredist_one(lat1,lon1,lat2,lon2,R,flag)

if strcmpi(flag(1:3),'sph')
    %disp('Spherical law of cosines')
    %  From http://mathworld.wolfram.com/GreatCircle.html 
    %     Numerically bad small distances due to acos of numerical noise!!
    d=R.*acos(cosd(lat1).*cosd(lat2).*cosd(lon1-lon2)+sind(lat1).*sind(lat2));
else
    %disp('Haversine formula')

    %This gives bad values when lon1 and lon2 are both near zero...
    %Haversin formula from http://www.movable-type.co.uk/scripts/latlong.html
    a1=squared(sind(frac(lat2-lat1,2)));
    a2=cosd(lat1).*cosd(lat2).*squared(sind(frac(lon2-lon1,2)));
    a=squared(sind(frac(lat2-lat1,2)))+cosd(lat1).*cosd(lat2).*squared(sind(frac(lon2-lon1,2)));
    %maxmax(a)-1
    a=min(a,1);%Sometimes a can exceed unity at the level of numerical noise
    d=R.*2.*atan2(sqrt(a),sqrt(1-a));
    %a3=squared(sind(frac(lon2-lon1,2)));
    %a3=frac(1,2)*(cosd(lat1-lat2)+cosd(lat1+lat2)).*squared(sind(frac(lon2-lon1,2)));
    %figure,plot([a(1:end-1) a2(1:end-1) a3(1:end-1) ])
    %figure,plot([lon2-lon1])
end

%cosAcosB=cos(A-B)+cos(A+B)
%e A+B = eA eB = cA+isA cB+isB = cAcB-sAsB = cA+B
%e A-B = eA e-B = cA+isA cB-isB = cAcB-sAsB =cA-B
%cA+B/2 +cA-B/2 = cAcB


%Sometimes there is a very small complex portion
d=real(d);

function[]=spheredist_test

try
    s=which('sw_dist');
    
    N=1000;
    tol=1;  %1 km
    
    lat1=180*rand(N,1)-90;    
    lon1=360*rand(N,1);   
    lat2=lat1+randn(N,1)-.5;       
    lon2=lon1+randn(N,1)-.5;
    
    index=find(lat2>-90&lat2<90);
    vindex(lat1,lon1,lat2,lon2,index,1);
            
    %Note that SW_DIST compares badly for large displacements 
    %lat2=180*rand(N,1)-90;    
    %lon2=360*rand(N,1);

    d1=spheredist(lat1,lon1,lat2,lon2);
    d2=0*lat1;
    for i=1:length(lat1)
        d2(i)=sw_dist([lat1(i) lat2(i)],[lon1(i) lon2(i)],'km');
    end
    b=aresame(d1,d2,tol);
    reporttest('SPHEREDIST versus SW_DIST for small displacements',b)
    
catch
    disp('SPHEREDIST test not run because SW_DIST not found.')
end

reporttest('SPHEREDIST equals zero for zero displacement', spheredist( 23.876, -63.397, 23.876, -63.397)==0);
%reporttest('SPHEREDIST equals zero for zero displacement', spheredist( 23.876, -63.397, 23.876, -63.397,'sphere')==0);
reporttest('SPHEREDIST equals pi * radearth for opposing points to within one meters',aresame(spheredist( 23.876, 180-63.397, -23.876, -63.397),radearth*pi,1e-3));


N=10000;
lat1=180*rand(N,1)-90;    
lon1=360*rand(N,1);   
lat2=lat1+randn(N,1)-.5;       
lon2=lon1+randn(N,1)-.5;
index=find(lat2>-90&lat2<90);
vindex(lat1,lon1,lat2,lon2,index,1);
d1=spheredist(lat1,lon1,lat2,lon2);
d2=spheredist(lat1,lon1,lat2,lon2,'sphere');
bool=(d1>1);
reporttest('SPHEREDIST haversine formula and spherical law of cosines match to within 1 mm for > 1 km distance',aresame(d1(bool),d2(bool),1e-6));


