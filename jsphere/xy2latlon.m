function[lat,lon]=xy2latlon(varargin)
% XY2LATLON  Converts tangent plane coordinates into latitude and longitude.
%
%   [LAT,LON]=XY2LATLON(X,Y,LATO,LONO) converts (X,Y) position with units 
%   of kilometers, specifying a location in a plane tangent to the earth 
%   at the point (LATO,LONO), into latitude and longitude locations.  
%
%   LONO may each either be specified on the interval [-180, 180] or on the
%   interval [0, 360].  Output longitudes LON are defined to be within the
%   interval [-180, 180].
%
%   LAT and LON may be computed with either the full spherical geometry,
%   the default, or using a small angle approximation.  To specify the
%   small angle approximation use XY2LATLON(...,'small').
%
%   LAT and LON are defined to be NAN for points with SQRT(X^2+Y^2)
%   exceeding the radius of the earth.
%
%   XY2LATLON(CX,LATO,LONO) with three input arguments, where CX is the 
%   complex-valued displacement CX=X+SQRT(-1)*Y, also works
%
%   The radius of the earth is given by the function RADEARTH.
%   ___________________________________________________________________
%
%   Cell array input / output
%
%   XY2LATLON returns cell array output given cell array input.  
%
%   That is, if the input arguments are all cell arrays of length K, 
%   containing K different numerical arrays, then the output will also be 
%   cell arrays of length K.  
%
%   This also works if X and Y, or alternatively CX, are cell arrays but 
%   LATO and LONO are scalars.
%   ___________________________________________________________________
%
%   XY2LATLON is inverted by LATLON2XY.
%
%   See also LATLON2XY, LATLON2UV.
%
%   Usage:  [lat,lon]=xy2latlon(x,y,lato,lono);
%           [lat,lon]=xy2latlon(cx,lato,lono);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2005--2018 J.M. Lilly --- type 'help jlab_license' for details        
  
na=nargin;
if ischar(varargin{end})
    str=varargin{end};
    na=na-1;
else 
    str='tangent';
end

if na==3
   if ~iscell(varargin{1})
       x=real(varargin{1});
       y=imag(varargin{1});
   else
       x=cellreal(varargin{1});
       y=cellimag(varargin{1});
   end
   lato=varargin{2};
   lono=varargin{3};
elseif na==4
   x=real(varargin{1});
   y=real(varargin{2});
   lato=varargin{3};
   lono=varargin{4};
end

if ~iscell(x)
   [lat,lon]=xy2latlon_one(x,y,lato,lono,str);
else
   if ~iscell(lato)
       lato1=lato;
       clear lato
       for i=1:length(x)
           lato{i}=lato1;
       end
   end
   if ~iscell(lono)
       lono1=lono;
       clear lono
       for i=1:length(x)
           lono{i}=lono1;
       end
   end
   lat=lato;
   lon=lono;
   for i=1:length(x)
       [lat{i},lon{i}]=xy2latlon_one(x{i},y{i},lato{i},lono{i},str);
   end
end

function[lat,lon]=xy2latlon_one(x,y,lato,lono,str)

[lato,lono]=jdeg2rad(lato,lono);

R=radearth;
index=find(sqrt(x.^2+y.^2)>=R);
if ~isempty(index)
    x(index)=nan;
    y(index)=nan;
end

if strcmpi(str(1:3),'sma')
    [lat,lon]=xy2latlon_cartesian(x,y,lato,lono,R);
elseif strcmpi(str(1:3),'tan')
    [lat,lon]=xy2latlon_sphere(x,y,lato,lono,R);
end


function[lat,lon]=xy2latlon_cartesian(x,y,lato,lono,R)

r1=R*cos(lato);

lat=y./R+lato;
lon=x./r1+lono;

%figure,plot(lat)

[lat,lon]=jrad2deg(lat,lon);

function[lat,lon]=xy2latlon_sphere(x,y,lato,lono,R)

r=sqrt(x.^2+y.^2);
%r  = distance in tangent plane

r1 = R - sqrt(R.^2 - r.^2);
%r1 = tangential distance from tangent plane to surface of earth
%     choosing smaller root, corresponding to near side of earth

%Now choose an xyz coordinate system with x=east, z= north

%tangent point = point on sphere at which plane is tangent
%contact point = point on sphere directly underneath (x,y) point in plane

R1 = sqrt((R-r1).^2+y.^2);
%R1 = distance from center of earth to contact point 
%projected onto the xz plane, i.e., looking down the y-axis

gamma=asin(frac(y,R1));
%gamma = angle spanned between contact point and tangent point
%projected onto the xz plane, i.e., looking down the y-axis

phi=lato+gamma;
%gamma = angle spanned between contact point and x-axis
%projected onto the xz plane, i.e., looking down the y-axis

xo=R1.*cos(phi);
zo=R1.*sin(phi);

yo=sqrt(R.^2-xo.^2-zo.^2);

index=find(x<0);
if ~isempty(x)
    yo(index)=-yo(index);
end
[lat,lon]=xyz2latlon(xo,yo,zo);
lon=jrad2deg(angle(rot(jdeg2rad(lon)+lono)));


