function[varargout]=cellimit(varargin)
%CELLIMIT  Limits the ranges of times in a cell array of numerical arrays.
%
%   NUMO=CELLIMIT(NUM,A,B), where NUM is a cell array of times, and A and B
%   are minimum and maximum times, returns NUMO in which each element is
%   limited to times not less than A and not greater than B.
%
%   As an example, with NUM given by 
%
%       NUM{1} = [1 2 3 4 5 6]';   NUM{2}=[3 4 5 6 7]';   NUM{3}=[5 6]';
%
%   calling NUMO=CELLIMIT(NUMI,4,6) will return
%
%       NUMO{1} = [4 5 6]';        NUMO{2}=[4 5 6]';      NUMO{3}=[5 6]';
%
%   such that that ranges of each element of NUMO are now between A and B.
%
%   A and B may also be arrays with the same number of elements as NUM.  In 
%   this case, the Nth entries A(N) and B(N) will set the range for NUM{N}.
%   
%   [NUMO,Y1,Y2,...,YN]=CELLIMIT(NUM,X1,X2,...,XN,A,B), with multiple input 
%   arguments, all cell arrays of numerical arrays having identical sizes
%   to NUM, will limit these to the same points as NUM. Thus Y1,Y2,...,YN
%   will be all of identical sizes to NUMO.
%
%   CELLIMIT(X1,X2,...XN); with no output arguments overwrites the 
%   original input variables.   
%
%   Usage: num=cellimit(num,a,b);
%          [num,y1,y2,y3]=cellimit(num,x1,x2,x3,a,b);
%          cellimit(num,x1,x2,x3,a,b);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    cellimit_test,return
end
if strcmp(varargin{1}, '--t')
    cellpack_test,return
end

b=varargin{end};
a=varargin{end-1};
varargin=varargin(1:end-2);

%Size check
for i=2:length(varargin)
    if size(varargin{i})~=size(varargin{1});
        error('All input arguments must be the same size.')
    end
end

num=varargin{1};
if length(a)==1
    a=a+zero(size(num));
end
if length(b)==1
    b=b+zero(size(num));
end

n0=0;
for i=1:length(num)
    n0=n0+length(num{i});
    index{i}=find(num{i}>=a(i),1,'first'):find(num{i}<=b(i),1,'last');
end

%Find all out-of-range
varargout=varargin;
for i=1:length(varargin)
    for j=1:length(varargin{1});
        varargout{i}{j}= varargin{i}{j}(index{j});
    end
end

disp(['Removing ' int2str(n0-sum(cellength(index))) ' out-of-range points.'])

nempty=length(find(cellength(index)==0));
if nempty>0
    disp(['Finding ' int2str(nempty) ' cells with no valid data points.'])
end

eval(to_overwrite(length(varargin)))

function[]=cellimit_test
 
%reporttest('CELLIMIT',aresame())
