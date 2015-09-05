function[varargout]=cellfirst(varargin)
%CELLFIRST  Returns the first (or last) element of each entry in a cell array.
%
%   Y=CELLFIRST(X), where X is a cell array of N numeric arrays, returns a
%   length N array Y containing the first element of each of the N cells.
%
%   An example of an input cell array X and the output Y is as follows:
%
%       X{1} = [1 2 3]';    X{2} = [3 4]';     X{3} = [5 6 7 ]';  
%       Y = [1 3 5]';
%
%   [Y1,Y2,...YN]=CELLFIRST(X1,X2,...XN) with multiple input arguments, all
%   of the same size, also works.  
%
%   CELLFIRST(X1,X2,...XN); with no output arguments overwrites the 
%   original input variables.   
%
%   CELLFIRST(...,'last') optionally returns the last element in each cell.
%
%   'cellfirst --t' runs a test.
%
%   Usage: y=cellfirst(x);
%          [y1,y2]=cellfirst(x1,x2);
%          [y1,y2]=cellfirst(x1,x2,'last');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    cellfirst_test,return
end
 
str='first';
if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
end

%Size check
for i=2:length(varargin)
    if size(varargin{i})~=size(varargin{1});
        error('All input arguments must be the same size.')
    end
end

%Find all nonfinite data
for i=1:length(varargin)
    varargout{i}=zeros(length(varargin{i}),1);
    for j=1:length(varargin{i})
        if isempty(varargin{i}{j})
            varargout{i}(j,1)=nan;
        else
            if strcmp(str(1:3),'fir')
                varargout{i}(j,1)=varargin{i}{j}(1);
            elseif strcmp(str(1:3),'las')
                varargout{i}(j,1)=varargin{i}{j}(end);
            end
        end
    end
end

eval(to_overwrite(nargin))  

function[]=cellfirst_test
 

x{1}=[1 4 3 2]';
x{2}=[2 3 4 5 6]';
x{3}=[5 7 10 nan 4 5]';
x{4}=[3 2 3 nan nan]';
x{5}=[nan nan nan]';
x{6}=[]';

y=cellmult(2,x);
x1=[1 2 5 3 nan nan]';
y1=2*[1 2 5 3 nan nan]';

cellfirst(x,y);

reporttest('CELLFIRST',aresame(x1,x)&&aresame(y1,y))

clear x
x{1}=[1 4 3 2]';
x{2}=[2 3 4 5 6]';
x{3}=[5 7 10 nan 4 5]';
x{4}=[3 2 3 nan nan]';
x{5}=[nan nan nan]';
x{6}=[]';

y=cellmult(2,x);
x1=[2 6 5 nan nan nan]';
y1=2*[2 6 5 nan nan nan]';

cellfirst(x,y,'last');


reporttest(['CELLFIRST with ''last'' flag'],aresame(x1,x)&&aresame(y1,y))
