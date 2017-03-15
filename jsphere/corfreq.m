function[fc]=corfreq(lat)
%CORFREQ  Coriolis frequency in radians per hour.
%
%   FC=CORFREQ(LAT) returns the Coliolis frequency at latitude LAT in
%   *radians* per hour.  FC is positive in the Northern Hemisphere and
%   negative in the Southern Hemisphere. 
%
%   The input argument LAT may also be a cell array of numeric arrays.  In 
%   this case FC will also be a cell array of the same size.
%
%   'corfreq --t' runs a test.
%
%   Usage: fc=corfreq(lat);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(lat, '--t')
    corfreq_test,return
end
 
if iscell(lat)
    fc=lat;
    for i=1:length(lat)
        fc{i}=corfreq_one(lat{i});
    end
else
   fc=corfreq_one(lat);
end

function[fc]=corfreq_one(lat)
omega=7.2921159e-5;
fc=2*sind(lat).*omega.*(3600);

function[]=corfreq_test
 
reporttest('CORFREQ at 30 degrees is about 2*pi/24',aresame(1./corfreq(30),24/2/pi,1e-1))
