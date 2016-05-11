function[varargout]=cellpack(varargin)
%CELLPACK  Removes all NaN values from cell arrays of numeric arrays.
%
%   Y=CELLPACK(X) where X is a cell array of numeric arrays, removes any
%   NaN values occuring within the cells and returns the result in Y.  
%
%   An example of an input cell array X and the output Y is as follows:
%
%       X{1} = [1 2 nan 4 nan]';  X{2}=[nan 2 3]';  X{3}=[1 2];   X{4}=[];
%       Y{1} = [1 2 4]';          Y{2}=[2 3]';      Y{3}=[1 2];   Y{4}=[];
%
%   CELLPACK does not remove empty cells.  This is done by CELLPRUNE.
%
%   [Y1,Y2,...YN]=CELLPACK(X1,X2,...XN) with multiple input arguments, all
%   of the same size, removes locations of NaNs that occur in *any* of the 
%   input variables.  The output variables Y1,Y2,...YN will all be the same
%   size, and none of them will contain any NaNs.
%
%   CELLPACK(X1,X2,...XN); with no output arguments overwrites the 
%   original input variables.   
%
%   'cellpack --t' runs a test.
%
%   Usage: y=cellpack(x);
%          [y1,y2,y3]=cellpack(x1,x2,x3);
%          cellpack(x1,x2,x3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015--2016 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    cellpack_test,return
end

%Size check
for i=2:length(varargin)
    if size(varargin{i})~=size(varargin{1});
        error('All input arguments must be the same size.')
    end
end
 
%Find all nonfinite data
for i=1:length(varargin)
    for j=1:length(varargin{1});
        if isempty(varargin{i}{j})
            bool{j}(:,i)=false;
        else
            bool{j}(:,i)=~isnan(varargin{i}{j});
        end
    end
end

n=vzeros(length(varargin{1}),1);
for j=1:length(varargin{1})
    n(j)=length(find(prod(bool{j},2)==0));
    index{j}=find(prod(bool{j},2)~=0);  %Finds all rows with zero bad columns
end

varargout=varargin;

nempty=0;
%Strip all nans top and bottom
for j=1:length(varargin{1})
    if isempty(index{j})
        nempty = nempty+1;
    end
    for i=1:length(varargin)
        if isempty(index{j})
            varargout{i}{j}=[];
        else
            varargout{i}{j}=varargin{i}{j}(index{j});
        end
    end
end
%n1,n2
disp(['Removing ' int2str(sum(n)) ' bad points.'])

if nempty>0
    disp(['Finding ' int2str(nempty) ' cells with no valid data points.'])
end

eval(to_overwrite(nargin))  

 
function[]=cellpack_test

x{1}=[1 4 3 2]';
x{2}=[1 3 2]';
x{3}=[1 7 10 nan 4 5]';
x{4}=[1 2 3 nan nan]';
x{5}=[nan nan nan]';
x{6}=[]';

y=cellmult(2,x);
y{2}(1)=nan;

x1{1}=[1 4 3 2]';
x1{2}=[3 2]';
x1{3}=[1 7 10 4 5]';
x1{4}=[1 2 3]';
x1{5}=[]';
x1{6}=[]';

y1=cellmult(2,x1);

cellpack(x,y);
reporttest('CELLPACK',aresame(x1,x)&&aresame(y1,y))


