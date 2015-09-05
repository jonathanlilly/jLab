function[varargout]=cellprune(varargin)
%CELLPRUNE  Removes all empty cells, or cells less than a specified length.
%
%   Y=CELLPRUNE(X), where X is a cell array, returns Y which is the same 
%   as X but excluding any empty cells.  
%
%   An example of an input cell array X and the output Y is as follows:
%
%       X{1} = [1 2 3 4 nan]';    X{2}=[]';    X{3}=[1 2];
%       Y{1} = [1 2 3 4]';                     Y{2}=[1 2];
%
%   Y=CELLPRUNE(X,M) will remove all cells with a length of less than M 
%   points. In the above example, Y=CELLPRUNE(X,3) gives Y{1}=[1 2 3 4]'.
%
%   [Y1,Y2,...YN]=CELLPRUNE(X1,X2,...XN) or CELLPRUNE(X1,X2,...XN,M) with
%   multiple input arguments, all of the same size, also works.  
%
%   CELLPRUNE(X1,X2,...XN); or CELLPRUNE(X1,X2,...XN,M);  with no output 
%   arguments overwrites the original input variables.  
%
%   'cellprune --t' runs a test.
%
%   Usage: y=cellprune(x);
%          [y1,y2,y3]=cellprune(x1,x2,x3);
%          [y1,y2,y3]=cellprune(x1,x2,x3,M);
%          cellprune(x1,x2,x3,M);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    cellprune_test,return
end

if ~iscell(varargin{end})
    L=varargin{end};
    varargin=varargin(1:end-1);
else
    L=1;
end
 
%Size check
for i=2:length(varargin)
    if size(varargin{i})~=size(varargin{1});
        error('All input arguments must be the same size.')
    end
end

len=cellength(varargin{1});
index=find(len>L); 
nempty=length(find(len==0));

%Find all nonfinite data
for i=1:length(varargin)
    varargout{i}=varargin{i}(index);
end

if nempty>0
    disp(['Excluding ' int2str(nempty) ' empty cells.'])
end

eval(to_overwrite(length(varargin))) 


function[]=cellprune_test
 
x{1}=[1 4 3 2]';
x{2}=[]';
x{3}=[1 7 10 nan 4 5]';
x{4}=[1 2 3 nan nan]';
x{5}=[nan nan nan]';
x{6}=[]';

y=cellmult(2,x);

x1=x([1 3 4 5]);
y1=y([1 3 4 5]);

cellprune(x,y);

reporttest('CELLPRUNE',aresame(x1,x)&&aresame(y1,y))

x1=x(2:3);
y1=y(2:3);

cellprune(x,y,4);

reporttest('CELLPRUNE',aresame(x1,x)&&aresame(y1,y))



