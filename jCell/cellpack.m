function[varargout]=cellpack(varargin)
%CELLPACK  Removes all INF values from cell arrays of numeric arrays.
%
%   Y=CELLPACK(X) where X is a cell array of numeric arrays, removes any
%   INF values occuring within the cells and returns the result in Y.  
%
%   An example of an input cell array X and the output Y is as follows:
%
%       X{1} = [1 2 inf 4 inf]';  X{2}=[inf 2 3]';  X{3}=[1 2];   X{4}=[];
%       Y{1} = [1 2 4]';          Y{2}=[2 3]';      Y{3}=[1 2];   Y{4}=[];
%
%   CELLPACK does not remove empty cells.  This is done by CELLPRUNE.
%
%   [Y1,Y2,...YN]=CELLPACK(X1,X2,...XN) with multiple input arguments, all
%   of the same size, removes locations of INFs that occur in *first* input
%   variable X1 only, while ignoring INFs in the others. The output
%   variables Y1,Y2,...YN will all be the same size.
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
%   (C) 2015--2020 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    cellpack_test,return
end

%Size check
for i=2:length(varargin)
    if size(varargin{i})~=size(varargin{1});
        error('All input arguments must be the same size.')
    end
end
 
% %Find all nonfinite data
% for i=1:length(varargin)
%     for j=1:length(varargin{1})
%         if isempty(varargin{i}{j})
%             bool{j}(:,i)=false;
%         else
%             bool{j}(:,i)=~isinf(varargin{i}{j});
%         end
%     end
% end

%Find all nonfinite data in only the first argument
for j=1:length(varargin{1})
    if isempty(varargin{1}{j})
        bool{j}=false;
    else
        bool{j}=~isinf(varargin{1}{j});
    end
end

for i=1:length(varargin)
    for j=1:length(varargin{1})
        if ~anyany(bool{j})
            varargout{i,1}{j,1}=[]; 
        else
            varargout{i,1}{j,1}=varargin{i}{j}(bool{j},:);
        end
    end
end
% n=vzeros(length(varargin{1}),1);
% for j=1:length(varargin{1})
%     n(j)=length(find(prod(bool{j},2)==0));
%     index{j}=find(prod(bool{j},2)~=0);  %Finds all rows with zero bad columns
% end
% 
% varargout=varargin;
% 
% nempty=0;
% %Strip all infs top and bottom
% for j=1:length(varargin{1})
%     if isempty(index{j})
%         nempty = nempty+1;
%     end
%     for i=1:length(varargin)
%         if isempty(index{j})
%             varargout{i}{j}=[];
%         else
%             varargout{i}{j}=varargin{i}{j}(index{j});
%         end
%     end
% end
%n1,n2
%disp(['Removing ' int2str(sum(n)) ' bad points.'])

%if nempty>0
%    disp(['Finding ' int2str(nempty) ' cells with no valid data points.'])
%end

eval(to_overwrite(nargin))  

 
function[]=cellpack_test


x{1}=[1 inf 3 2]';
x{2}=[1 3 2]';
x{3}=[1 7 10 inf 4 5]';
x{4}=[1 2 3 inf inf]';
x{5}=[inf inf inf]';
x{6}=[]';

y{1}=[1 4 3 inf]';
y{2}=[1 inf 2]';
y{3}=[1 7 10 5 4 5]';
y{4}=[1 2 3 2 4]';
y{5}=[inf inf 2]';
y{6}=[]';

x1{1}=[1 3 2]';
x1{2}=[1 3 2]';
x1{3}=[1 7 10 4 5]';
x1{4}=[1 2 3]';
x1{5}=[]';
x1{6}=[]';

y1{1}=[1 3 inf]';
y1{2}=[1 inf 2]';
y1{3}=[1 7 10 4 5]';
y1{4}=[1 2 3]';
y1{5}=[]';
y1{6}=[]';

cellpack(x,y);
reporttest('CELLPACK',aresame(x1,x)&&aresame(y1,y))

% for i=1:6
%     [x1{i} x{i}]
% end


