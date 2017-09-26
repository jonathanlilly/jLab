function[y]=degunwrap(varargin)
%DEGUNWRAP  Unwraps arrays given in degrees.
%
%   DEGUNWRAP(D) unwraps array D, given in degrees, along its rows. 
%
%   DEGUNWRAP(D,TOL) uses tolerance TOL, also in degrees.  The default
%   value of TOL is 180, corresponding to the value PI used by UNWRAP.
%
%   The unwrapping may be reversed by DEG180 or DEG360.
%
%   The input argument D may also be a cell array of numerical arrays. 
%   In this case, each element of D is unwrapped along its rows. 
%
%   This is useful for unwrapping longitudes, for example, the surface
%   drifter trajectories in DRIFTERS.MAT described in ABOUT_DRIFTERS.
%   __________________________________________________________________
%   Parallelization
%
%   DEGUNWRAP(D,'parallel'), when D is a cell array, parallelizes the 
%   computation using a PARFOR loop.  This requires that Matlab's Parallel 
%   Computing Toolbox be installed, and is useful for very large datasets.
%   __________________________________________________________________
%
%   See also JDEG2RAD, JRAD2DEG, DEG180, DEG360.
%
%   Usage: lon=degunwrap(lon);
%          lon=degunwrap(lon,tol);
%          lon=degunwrap(lon,tol,'parallel');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
 
cores='serial';

if ischar(varargin{end})
    if strcmpi(varargin{end}(1:3),'par')||strcmpi(varargin{end}(1:3),'ser')
        cores=varargin{end};
    end
    varargin=varargin(1:end-1);
end

if strcmpi(cores(1:3),'par')
    if exist('parpool')~=2
        disp('Sorry, parallel algorithm requires the Parallel Computing Toolbox.')
        disp('Defaulting to the standard algorithm.')
        cores='serial';
    end
end

x=varargin{1};
if length(varargin)==1
    tol=180;
else
    tol=varargin{2};
end
tol=frac(2*pi,360)*tol;

if ~iscell(x)
    y=frac(360,2*pi)*unwrap(frac(2*pi,360)*x,tol);
else
    if strcmpi(cores(1:3),'par')
        parfor i=1:length(x)
            y{i,1}=frac(360,2*pi)*unwrap(frac(2*pi,360)*x{i},tol);
        end
    else
        for i=1:length(x)
            y{i,1}=frac(360,2*pi)*unwrap(frac(2*pi,360)*x{i},tol);
        end
    end 
end
        
