function[varargout]=cellgrid(varargin)
%CELLGRID  Interpolate a cell array of numeric arrays onto a regular grid. 
%
%   [TO,Y]=CELLGRID(T,X,DT) where T is a cell array of time arrays, and
%   X is of the same size as T, linearly interpolates the elements of X 
%   within T at times TO, which are regularly spaced with interval DT.  
%
%   CELLGRID does not modify bad data points, such as those marked with
%   NaNs.  See CELLFILL to interpolate over bad data points.  
%
%   As an example, with T and X given by
%
%       T{1} = [1 2 3  5  6]';   T{2}=[3 7  9 10]';   
%       X{1} = [2 4 6 10 12]';   X{2}=[5 9 11 12]';
%   
%   [TO,Y]=CELLGRID(1,T,X) will return
%
%       TO{1} = [1 2 3 4  5  6]';   TO{2} = [3 4 5 6 7 8  9  10]';   
%       Y{1}  = [2 4 6 8 10 12]';   Y{2}  = [5 6 7 8 9 10 11 12]';
%
%   By default, CELLGRID uses INTERP1 with linear interpolation. 
%   CELLGRID(...,STR) instead uses the method specified by STR, e.g. 
%   'pchip'.  See INTERP1 for details. 
%
%   [TO,Y1,Y2,...YN]=CELLGRID(T,X1,X2,...XN,DT) with multiple input
%   arguments also works provided the XN are all the same size. 
%
%   CELLGRID(X1,X2,...XN); with no output arguments overwrites the 
%   original input variables. 
%
%   'cellgrid --t' runs a test.
%
%   Usage: [to,y]=cellgrid(t,x,dt);
%          [to,y1,y2,y3]=cellgrid(t,x1,x2,x3,dt);
%          cellgrid(t,x1,x2,x3,dt);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    cellgrid_test,return
end

str='linear';
if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
end

dt=varargin{end};
varargin=varargin(1:end-1);

t=varargin{1};
 
%Size check
for i=2:length(varargin)
    if size(varargin{i})~=size(varargin{1});
        error('All input arguments must be the same size.')
    end
end

for i=1:length(t)
   to{i,1}=[min(t{i}):dt:max(t{i})]';
end

for j=2:length(varargin)
   for i=1:length(t)
        if length(t{i})>1
            varargout{j}{i,1}=interp1(t{i},varargin{j}{i},to{i},str);
        else
            varargout{j}{i,1}=varargin{j}{i};
        end
   end
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
