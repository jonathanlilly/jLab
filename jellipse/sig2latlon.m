function[latres,lonres,lathat,lonhat]=sig2latlon(x,y,lato,lono)
%SIG2LATLON  Converts an oscillatory signal to lat/lon displacements.
%
%   [LATRES,LONRES]=SIG2LATLON(X,Y,LAT,LON) where X and Y are zonal and
%   meridional displacements, in kilometers, about some time-varying 
%   position with latitude LAT and longitude LON, returns the residuals
%   between the signal and the position curve.
%
%   [LATRES,LONRES,DLAT,DLON]=SIG2LATLON(X,Y,LAT,LON) also returns the
%   latitude and longitude displacements associated with X and Y.  
%
%   The residuals are output first because this is usually what is desired.
%   These are essentially LAT-DLAT and LON-DLON, apart from an unwrapping 
%   correction for longitude.
%
%   All input arguments are arrays of the same size.  
%   ____________________________________________________________________
%   
%   Cell array input/output
%
%   If SIG2LATLON is given cell array input, it returns cell array output.
%
%   Thus X, Y, LAT, and LON may each be cell arrays of the same size, 
%   where each element in the cell array is a numerical array. The output
%   arguments will then also be cell arrays of this size.
%   ____________________________________________________________________
%
%   'sig2latlon --t' runs a test.
%
%   Usage: [latres,lonres]=sig2latlon(x,y,lat,lon);
%          [latres,lonres,dlat,dlon]=sig2latlon(x,y,lat,lon);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2016 J.M. Lilly --- type 'help jlab_license' for details
 
if ~isempty(x)
    if ~iscell(x)
        [lathat,lonhat,latres,lonres]=sig2latlon_one(x,y,lato,lono);
    else
        for i=1:length(x)
            [lathat{i,1},lonhat{i,1},latres{i,1},lonres{i,1}]=sig2latlon_one(x{i},y{i},lato{i},lono{i});
        end
    end
else
    latres=x;
    lonres=x;
    lathat=x;
    lonhat=x;
end

function[lathat,lonhat,latres,lonres]=sig2latlon_one(x,y,lato,lono)
x=real(x);
y=real(y);

lathat=frac(360,2*pi)*frac(y,radearth);
lonhat=deg180(frac(360,2*pi)*frac(x,radearth*cosd(lato)));
latres=lato-lathat;
lonres=deg180(frac(360,2*pi)*angle(rot(frac(2*pi,360)*(lono-lonhat))));

    
    