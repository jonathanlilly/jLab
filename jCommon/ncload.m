function[varargout]=ncload(varargin)
%NCLOAD  Load all variables from a NetCDF file and convert trajectories to cells.
%
%   NCLOAD(FILENAME) loads all variables from the NetCDF file FILENAME and
%   places them into a structure of the same name in the calling workspace. 
%   The extension '.nc' in FILENAME is optional.
%
%   If FILENAME contains a full pathname, only the final portion after the
%   last '/' or '\' is used for the structure name.
%
%   Note that any hyphens '-' in the filename are replaced with underscores
%   in the structure name, as the former are not allowed in variable names.
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
%   NCLOAD will convert any variables that are numeric but not doubles to
%   doubles.  This is because experience has shown using other data types 
%   in code that is expecting doubles can lead to errors in Matlab.
%   __________________________________________________________________
%
%   Convert trajectory data
%
%   In addition, if FILENAME has the global attribute  
%
%         featureType    = 'trajectory'
%
%   then any variables in FILENAME are have 'Dimensions: obs' will be
%   interpreted as concatentated trajectory data, following the NetCDF 
%   Climate and Forecast (CF) Metadata Conventions Appendix H.4.3
%
%         http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html
%
%   NCLOAD adopts the convention that a variable named 'ids', having 
%   'Dimensions: obs', provides the id for all other variables.  
%
%   All variables of dimension obs will be converted to cell arrays by
%   first calling COLBREAKS using ids as the first argument, then COL2CELL.
%
%   Alternatively, all variables with 'Dimensions: obs'  will be 
%   interpreted as being in the NaN-separated column format used by
%   COL2CELL, which will be called to put then in cell array format. This
%   format is supported to provide for reverse compatibility; it is 
%   recommended to use the CF trajectory convention instead. 
%   __________________________________________________________________
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
%   (C) 2019--2020 J.M. Lilly --- type 'help jlab_license' for details


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
structname(real(structname)==real('-'))='_';  %Replace minus signs with underscores

            
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

try
    attvalue=ncreadatt(filename,'/','featureType');
catch
    attvalue='nontrajectory';
end

convert=false;
if strcmp(attvalue,'trajectory')
    convert=true;
    isobs=false(length(dimensions),1);
    isids=false(length(dimensions),1);
    for i=1:length(isobs)
        isobs(i)=strcmp(dimensions{i},'obs');
        isids(i)=strcmp(vars{i},'ids');
    end
    if anyany(isobs)
        if ~anyany(isids)
            convert=false;
            disp('NCLOAD finding trajectory data but no "ids" variable, so not converting.')
        end
    end
    
end

if ~isempty(dimensions)
    for i=1:length(vars)
        if ~isempty(dimensions{i})
            %infostruct.Variables(i).Dimensions.Name
            if strcmp(dimensions{i},'column')
                %converting my older 'column' format
                str=[structname '.' vars{i} '=doubifnum_internal(col2cell_internal(ncread(''' filename ''',''' vars{i} ''')));'  ];
            else
                str=[structname '.' vars{i} '=doubifnum_internal(ncread(''' filename ''',''' vars{i} '''));'  ];
            end
            evalin('caller',str)
        end
    end
end

if convert
    strlist=[structname '.ids'];%make sure ids comes firsts
    for i=1:length(vars)
        if strcmp(dimensions{i},'obs')&&~strcmp(vars{i},'ids')
            strlist=[strlist ',' structname '.' vars{i} ];
        end
    end
    evalin('caller',['[' strlist ']=colbreaks_internal(' strlist ');'])
    evalin('caller',['[' strlist ']=col2cell_internal(' strlist ');'])
end

function[x]=doubifnum_internal(x)
%DOUBIFNUM  Converts a variable to double precision if it is numeric.
%
%   X=DOUBIFNUM
%
%   Usage: x=doubifnum(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2020 J.M. Lilly --- type 'help jlab_license' for details
 
if isnumeric(x)
    x=double(x);
end
    
function[varargout]=colbreaks_internal(varargin)
%COLBREAKS  Insert NANs into discontinuties in a vector.
%
%   O1=COLBREAKS(I1), where I1 is a column vector, inserts NANs into I1 at
%   its discontinuities (and puts one at the end), for example:
%		
%	  COLBREAKS([2;2;2;3;3])=[2;2;2;NAN;3;3;NAN];
%
%   [O1,O2,O3,...]=COLBREAKS(I1,I2,...), where I1,I2,... are all column 
%   vectors of the same length, uses I1 as a reference for the other 
%   vectors and inserts NANs into all the vectors at locations where I1 is 
%   discontinous.
%	
%   For instance, using station number for I1 will sort a column of CTD 
%   data into a NAN-padded matrix.
%
%   Note that for complex-valued input arrays, a complex-valued NAN, 
%   NAN+SQRT(-1)*NAN, is used to indicate the break locations. 
%
%   MAT2=COLBREAKS(MAT1), where MAT1 and MAT2 are matrices of the same 
%   size, also works.  In this case the first column of MAT1 is used as the
%   reference vector.
%
%   COLBREAKS, COL2MAT, and MAT2COL together form a system for moving data
%   with segments of nonuniform length rapidly back and forth between a 
%   column format and a padded-matrix format. CTD or float data, for 
%   instance, can be stored in the (usually much smaller) column format and
%   converted into the matrix format upon loading.
%  
%   COLBREAKS(C1,C2,...); with no output arguments overwrites the original
%   input variables.  
%
%   See also COL2MAT, MAT2COL, COL2CELL, CELL2COL, ORBITBREAKS.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2000--2013 J.M. Lilly --- type 'help jlab_license' for details

bmat=0;
nargs=nargin;
ref=varargin{1};
if nargin==1 && size(ref,2)>1
    temp=ref;
    ref=temp(:,1);
    c=cell(size(temp,2),1);
    for i=1:size(temp,2)-1
        c{i}=temp(:,i+1);
    end
    nargs=size(temp,2);
    clear temp
    bmat=1;
else
    c=varargin(2:end);
end



ii=0;
index=find(diff(ref)~=0)+1;
bool=false(size(ref));
bool(index)=ones(size(index));
index=(1:length(ref))'+cumsum(bool);

while ii<nargs - 1
    ii=ii+1;
    col=c{ii};
    
    if isreal(col)
        colout=nan*ones(max(index),1);
        colout(index)=col;
        colout=[colout;nan];
    else
        colout=(nan+sqrt(-1)*nan)*ones(max(index),1);
        colout(index)=col;
        colout=[colout;(nan+sqrt(-1)*nan)];
    end
    
    varargout{ii+1}=colout;
end

if isreal(ref)
    refout=nan*ones(max(index),1);
    refout(index)=ref;
    refout=[refout;nan];
else
    refout=(nan+sqrt(-1)*nan)*ones(max(index),1);
    refout(index)=ref;
    refout=[refout;(nan+sqrt(-1)*nan)];
end

varargout{1}=refout;
if bmat
    for i=1:nargs-1
        varargout{1}(:,i+1)=varargout{i+1};
    end
    varargout=varargout(1);
end

function[varargout]=col2cell_internal(varargin)
%COL2CELL  Converts 'column-appended' data into cell arrays of numeric arrays.
%
%   X=COL2CELL(COL) converts the array COL, a 2D array having blocks of 
%   data separated by rows of all NaNs, into a cell array, with each cell
%   containing one of the data blocks, excluding the trailing NaNs.
%
%   Since NANs mark the end of each block, missing or bad data in COL
%   should instead be indicated by INFs.
%
%   [X1,X2,...,XN]=COL2CELL(C1,C2,...,CN) also works, where the input 
%   fields C1,...,CN are all the same size.  In this case the locations of 
%   NANs in C1 are used to break all the other CN into cells.  
%
%   COL2CELL(C1,C2,...,CN); with no output arguments overwrites the 
%   original variables.
%
%   COL2CELL is inverted by CELL2COL.
%   __________________________________________________________________
%
%   See also CELL2COL, COL2MAT, MAT2COL, COLBREAKS.
%
%   'col2cell --t' runs a test.
%
%   Usage: x=col2cell(col);
%          [x1,x2,...,xN]=col2cell(c1,c2,...,cN);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2018 J.M. Lilly --- type 'help jlab_license' for details

i=0;
if ~isempty(varargin{1})
    while i<nargin
        i=i+1;
        col=varargin{i};
        if i==1
            %Note minor changes to make this work on matrices rather than
            %just columns. 
            ia=[1;find(isnan(col(1:end-1,1)))+1];
            ib=find(isnan(col(:,1)))-1;
            %[L,ia,ib]=blocklen(col);
            ray=cell(length(ia),1);
        end
        for j=1:length(ia)
            if ~isempty(col)
                ray{j}=col(ia(j):ib(j),:);
            else
                ray{j}=[];
            end
        end
        varargout{i}=ray;
    end
else
    for i=1:nargout
        varargout{i}=varargin{1};
    end
end



