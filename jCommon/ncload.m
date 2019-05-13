function[varargout]=ncload(varargin)
%NCLOAD  Load all variables from a NetCDF file and convert columns to cells.
%
%   NCLOAD(FILENAME) loads all variables from the NetCDF file FILENAME and
%   places them into a structure of the same name in the calling workspace. 
%   The extension '.nc' in FILENAME is optional.
%
%   If FILENAME contains a full pathname, only the final portion after the
%   last '/' or '\' is used for the structure name.
%
%   For example, if FILENAME='/Home/data' contains variables 'num', 'lat', 
%   and 'lon', then the result of calling NCLOAD('data') will be a 
%   structure 'data' with fields 'data.num', 'data.lat', and 'data.lon'.  
%
%   Then USE can be used to map these variables into the main workspace.  
%
%   NCLOAD(FILENAME,VAR1,VAR2,...,VARN) will only load the variables with
%   the names VAR1, VAR2,...VARN.
%
%   In addition, if any variables in FILENAME are of type 'column', these
%   will be interpreted as being in the NaN-separated column format used by
%   COL2CELL, which will be called to put then in cell array format. 
%
%   The following formats also work:
%
%          ncload filename
%          ncload filename var1 var2 ... varN
%   
%   In this format, the input strings have to be the actual names of the
%   file and variables, as opposed to variables containing those names.
%
%   This is basically a way to conveniently work with small NetCDF files as
%   if they were mat-files.  Typically, NetCDF files load much faster. 
%
%   Usage: ncload(filename);
%          ncload filename 
%          ncload(filename,'num','lat','lon');
%          ncload filename num lat lon 
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2019 J.M. Lilly --- type 'help jlab_license' for details


filename=varargin{1};
structname=filename;
ii=findstr(structname,'.');
if isempty(ii)
    filename=[filename '.nc'];
else
    structname=structname(1:ii-1);
end    
index=findstr(structname,'/');
if ~isempty(index)
    structname=structname(index(end)+1:end);
end
index=findstr(structname,'\');
if ~isempty(index)
    structname=structname(index(end)+1:end);
end

vars=varargin(2:end);
infostruct=ncinfo(filename);

%infostruct.Dimensions

for i=1:length(infostruct.Variables)
    allvars{i}=infostruct.Variables(i).Name;
    alldimensions{i}=infostruct.Variables(i).Dimensions.Name;
end

%alldimensions

dimensions=[];
if isempty(vars)
    vars=allvars;
    dimensions=alldimensions;
else
    for i=1:length(vars)
        for j=1:length(allvars)
            if strcmp(vars{i},allvars{j})
                 dimensions{i}=alldimensions{j};
            end
        end
    end
end

%vars
%dimensions

if ~isempty(dimensions)
    for i=1:length(vars)
        if ~isempty(dimensions{i})
            %infostruct.Variables(i).Dimensions.Name
            if strcmp(dimensions{i},'column')
                str=[structname '.' vars{i} '=col2cell(ncread(''' filename ''',''' vars{i} '''));'  ];
            else
                str=[structname '.' vars{i} '=ncread(''' filename ''',''' vars{i} ''');'  ];
            end
            evalin('caller',str)
        end
    end
end


