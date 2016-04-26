function[varargout]=cellfill(varargin)
%CELLFILL  Fills missing data marked by NaNs in a cell array.
%
%   Y=CELLFILL(X) where X is a cell array of numeric arrays, interpolates 
%   to fill NaN values in the of the cells, and returns the result in Y.  
%
%   By default CELLFILL uses the 'pchip' method.  CELLFILL(X,STR) instead
%   uses the method specified by STR.  See INTERP1 for possible methods. 
%
%   An example of an input cell array X and the output Y is as follows:
%
%       X{1} = [1 2 nan 4 5]';   X{2}=[2 4 6 nan nan 12]';   
%       Y{1} = [1 2  3  4 5]';   Y{2}=[2 4 6  8  10  12]';   
%
%   If NaNs in X are at the beginnings or the ends of cells, some NaNs may 
%   still be present in Y, depending on the method used.  If all NaNs in X 
%   are interior NaNs, there will be no NaNs in Y regardless of the method.
%   See also CELLSTRIP for dealing with leading or trailing NaNs.
%
%   [Y1,Y2,...YN]=CELLFILL(X1,X2,...XN) with multiple input arguments also
%   works provided the XN are all the same size. 
%
%   CELLFILL(X1,X2,...XN); with no output arguments overwrites the 
%   original input variables.   
%
%   'cellfill --t' runs a test.
%
%   Usage: y=cellfill(x);
%          [y1,y2,y3]=cellfill(x1,x2,x3);
%          cellfill(x1,x2,x3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    cellfill_test,return
end

str='pchip';
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
    for j=1:length(varargin{1});
        temp=varargin{i}{j};
        if  length(find(~isnan(temp)))==0
            bool{j}(:,i)=false;
            varargout{i}{j,1}=[];
        else
            bool{j}(:,i)=~isnan(temp);
            index=find(bool{j}(:,i));
            if length(temp)>1
               varargout{i}{j,1}=interp1(index,temp(index),[1:length(temp)]',str);
            else
               varargout{i}{j,1}=temp(index);
            end
        end
    end
end

for j=1:length(bool)
    bool{j}=prod(bool{j},2);
end
    
eval(to_overwrite(nargin))  
varargout{length(varargin)+1}=bool;


function[]=cellfill_test
 

x{1}=[1 4 3 2]';
x{2}=[1 3 2]';
x{3}=[1 2 nan 4 5]';
x{4}=conj([2+1i 4+1i 6+1i nan nan 12+1i])';
x{5}=[nan nan nan]';
x{6}=[]';

y=x;
x1=x;
x1{3}=[1 2 3 4 5]';
x1{4}=[2 4 6 8 10 12]'+1i;
x1{5}=[];
y1=x1;

cellfill(x,y,'linear');

reporttest('CELLFILL',aresame(x1,x)&&aresame(y1,y))

