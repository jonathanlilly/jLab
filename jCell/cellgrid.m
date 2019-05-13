function[varargout]=cellgrid(varargin)
%CELLGRID  Interpolate a cell array of numeric arrays onto a regular grid. 
%
%   [TO,Y]=CELLGRID(T,X,DT) where T is a cell array of time arrays, and
%   X is of the same size as T, linearly interpolates the elements of X 
%   within T at times TO, which are regularly spaced with interval DT.  
%   
%   DT may either be a scalar, or an array of the same length as T and X.
%
%   CELLGRID does not modify bad data points, such as those marked with
%   NaNs.  See CELLFILL to interpolate over bad data points.  
%
%   As an example, with T and X given by
%
%       T{1} = [1 2 3  5  6]';   T{2}=[3 7  9 10]';   
%       X{1} = [2 4 6 10 12]';   X{2}=[5 9 11 12]';
%   
%   [TO,Y]=CELLGRID(T,X,1) will return
%
%       TO{1} = [1 2 3 4  5  6]';   TO{2} = [3 4 5 6 7 8  9  10]';   
%       Y{1}  = [2 4 6 8 10 12]';   Y{2}  = [5 6 7 8 9 10 11 12]';
%
%   By default, CELLGRID uses INTERP1 with the 'pchip' method of
%   interpolation.  CELLGRID(...,STR) instead uses the method specified by
%   STR, e.g. STR='linear'.  See INTERP1 for details. 
%
%   [TO,Y1,Y2,...YN]=CELLGRID(T,X1,X2,...XN,DT) with multiple input
%   arguments also works provided the XN are all the same size. 
%
%   CELLGRID(X1,X2,...XN); with no output arguments overwrites the 
%   original input variables. 
%   __________________________________________________________________
%
%   Specifying interpolated times
%
%   [TO,Y]=CELLGRID(T,X,DT,A,B) will set TO to TO=A:DT:B, where A and B 
%   are either scalars or an arrays of the same size as X. 
%
%   The default behavior is equivalent to choosing A and B as the first and
%   last elements of TO, that is, to setting A=CELLMIN(T) and B=CELLMAX(T).
%   __________________________________________________________________
%
%   Parallelization
%
%   CELLGRID(...,'parallel') parallelizes the computation using a PARFOR 
%   loop over the various input variables.  This requires that Matlab's 
%   Parallel Computing Toolbox be installed. 
%   __________________________________________________________________
%
%   'cellgrid --t' runs a test.
%
%   Usage: [to,y]=cellgrid(t,x,dt);
%          [to,y1,y2,y3]=cellgrid(t,x1,x2,x3,dt);
%          [to,y1,y2,y3]=cellgrid(t,x1,x2,x3,dt,a,b);
%          cellgrid(t,x1,x2,x3,dt);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015--2019 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    cellgrid_test,return
end

str='pchip';
cores='serial';

for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'ser')||strcmpi(varargin{end}(1:3),'par')
            cores=varargin{end};
        else
            str=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end

if strcmpi(cores(1:3),'par')
    if exist('parpool')~=2
        disp('Sorry, parallel algorithm requires the Parallel Computing Toolbox.')
        disp('Defaulting to the serial algorithm.')
        str='serial';
    end
end

t=varargin{1};

if ~iscell(varargin{end-2})&&~ischar(varargin{end-2})
    a=varargin{end-1};
    b=varargin{end};
    if length(a)==1
        a=a+zeros(size(t));
    end
    if length(b)==1
        b=b+zeros(size(t));
    end
    varargin=varargin(1:end-2);
else 
    a=cellmin(t);
    b=cellmax(t);
end
    
dt=varargin{end};

varargin=varargin(1:end-1);
 
%Size check
for i=2:length(varargin)
    if size(varargin{i})~=size(varargin{1});
        error('All input arguments must be the same size.')
    end
end

if length(dt)==1
    dt=dt+0*a;
end

for i=1:length(t)
   to{i,1}=[a(i):dt(i):b(i)]';
end

if strcmpi(cores(1:3),'par')
    disp('CELLGRID employing parallel algorithm.')
    parfor j=2:length(varargin)
        for i=1:length(t)
            if length(t{i})>1
                bool=~isnan(t{i})&~isnan(varargin{j}{i});
                if length(find(bool))~=0
                    args{j}{i,1}=interp1(t{i}(bool),varargin{j}{i}(bool),to{i},str);
                else
                    args{j}{i,1}=nan*to{i};
                end
            else
                args{j}{i,1}=varargin{j}{i};
            end
        end
    end
else
    for j=2:length(varargin)
        for i=1:length(t)
            if length(t{i})>1
                bool=~isnan(t{i})&~isnan(varargin{j}{i});
                if length(find(bool))>=2
                    %i,j,length(find(bool))
                    args{j}{i,1}=interp1(t{i}(bool),varargin{j}{i}(bool),to{i},str);
                else
                    args{j}{i,1}=nan*to{i};
                end
            else
                args{j}{i,1}=varargin{j}{i};
            end
        end
    end
end

varargout=args;
if ~isempty(varargout{2})
    for i=1:length(to)
        if isempty(varargout{2}{i})
            to{i}=[];
        end
    end
else
    to=[];
end
varargout{1}=to;

eval(to_overwrite(length(varargin)))

function[]=cellgrid_test
 
t{1}=[1 2 3  5  6]';
t{2}=[3 7  9 10]';
x{1}=[2 4 6 10 12]';
x{2}=[5 9 11 12]';

 
to{1}=[1 2 3 4  5  6]';
to{2}=[3 4 5 6 7 8  9  10]';
y{1}=[2 4 6 8 10 12]';
y{2}=[5 6 7 8 9 10 11 12]';

cellgrid(t,x,1);

reporttest('CELLGRID',aresame(x,y)&&aresame(to,t))
