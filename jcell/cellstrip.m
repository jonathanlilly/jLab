function[varargout]=cellstrip(varargin)
%CELLSTRIP  Strips NaN values from the beginnings or ends of cell arrays.
%
%   Y=CELLSTRIP(X) where X is a cell array of numeric arrays, strips any
%   NaN values at the beginnings or the ends of the cells and returns the 
%   result in Y.  Thus Y will begin and end with non-NaN values.
%
%   An example of an input cell array X and the output Y is as follows:
%
%       X{1} = [1 2 nan 4 nan]';  X{2}=[nan 2 3]';  X{3}=[1 2];   X{4}=[];
%       Y{1} = [1 2 nan 4]';      Y{2}=[2 3]';      Y{3}=[1 2];   Y{4}=[];
%
%   CELLSTRIP does not remove empty cells.  This is done by CELLPRUNE.
%
%   [Y1,Y2,...YN]=CELLSTRIP(X1,X2,...XN) with multiple input arguments, all
%   of the same size, strips locations of leading or trailing NaNs that 
%   occur in *any* of the input variables.  Thus the output variables will
%   all be the same size, and none of them will begin or end with a NaN.
%
%   CELLSTRIP(X1,X2,...XN); with no output arguments overwrites the 
%   original input variables.   
%
%   'cellstrip --t' runs a test.
%
%   Usage: y=cellstrip(x);
%          [y1,y2,y3]=cellstrip(x1,x2,x3);
%          cellstrip(x1,x2,x3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    cellstrip_test,return
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

[n1,n2]=vzeros(length(varargin{1}),1);
for j=1:length(varargin{1})
    index{j}=find(prod(bool{j},2)~=0);  %Finds all rows with zero bad columns
    if isempty(index{j})
        n1(j)=0;
        n2(j)=0;
    else
        n1(j)=index{j}(1)-1;
        n2(j)=length(varargin{1}{j})-index{j}(end);
    end
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
            varargout{i}{j}=varargin{i}{j}(index{j}(1):index{j}(end));
        end
    end
end
%n1,n2
disp(['Removing ' int2str(sum(n1)) ' bad points from beginnings of cells.'])
disp(['Removing ' int2str(sum(n2)) ' bad points from ends of cells.'])
disp(['Finding ' int2str(nempty) ' cells with no valid data points.'])

eval(to_overwrite(nargin))

function[]=cellstrip_test
 

x{1}=[1 4 3 2]';
x{2}=[1 3 2]';
x{3}=[1 7 10 nan 4 5]';
x{4}=[1 2 3 nan nan]';
x{5}=[nan nan nan]';
x{6}=[]';

y=x;
y{2}(1)=nan;

x1=x;
y1=y;

x1{2}=x{2}(2:3);
y1{2}=y{2}(2:3);
x1{4}=x{4}(1:3);
y1{4}=y{4}(1:3);
x1{5}=[];
y1{5}=[];
x1{6}=[];
y1{6}=[];

cellstrip(x,y);

reporttest('CELLSTRIP',aresame(x1,x)&&aresame(y1,y))

