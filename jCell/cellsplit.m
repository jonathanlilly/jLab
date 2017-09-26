function[varargout]=cellsplit(varargin)
%CELLSPLIT  Splits cell arrays of numeric arrays at data gaps.
%
%   TO=CELLSPLIT(T,TOL), where T is a cell array of time arrays, splits T
%   into additional cells wherever gaps of longer than duration TOL occur.
%   Thus TO contains no gaps of longer than duration TOL. 
%
%   As an example, if T is a length two cell array 
%   
%       T{1}=[1 2 3 8 9 12 17]';  T{2} =[1 6 7]'; 
%
%   then TO=CELLSPLIT(TI,3) leads to 
%
%       TO{1}=[1 2 3]'; TO{2}=[8 9 12]'; TO{3}=17; TO{4}=1; TO{5}=[6 7]'; 
% 
%   such that there are no gaps in TO greater than 3. 
%
%   In the process of this splitting, any empty cells in T are removed.  
%
%   [TO,Y1,Y2,...,YN]=CELLSPLIT(T,X1,X2,...,XN,TOL) where the XN are cell
%   arrays of the same size as T, splits the XN at the *same locations* at
%   which T is split, and returns the results in the YN.  
%
%   The purpose of this is to separate all variables XN into separate cells 
%   whenever data gaps are longer than a specified cutoff.
%   
%   CELLSPLIT(T,X1,X2,...,XN,TOL); with no output arguments overwrites
%   the original input variables T and X1,X2,...XN.
%
%   'cellsplit --t' runs a test.
%
%   Usage: to=cellsplit(t,tol);
%          [to,y1,y2,y3]=cellsplit(t,x1,x2,x3,tol);
%          cellsplit(t,x1,x2,x3,tol);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    cellsplit_test,return
end

tol=varargin{end};
varargin=varargin(1:end-1);

%Size check
for i=2:length(varargin)
    if size(varargin{i})~=size(varargin{1});
        error('All input arguments must be the same size.')
    end
end

num=varargin{1};
n=0;
for i=1:length(num)
    index=find(diff(num{i})>tol);
    ia=[1;index+1];
    ib=[index;length(num{i})];
    for k=1:length(ia)
        n=n+1;
        for j=1:length(varargin)
            varargout{j}{n,1}=varargin{j}{i}(ia(k):ib(k));
        end
    end
end

eval(to_overwrite(length(varargin))) 


function[]=cellsplit_test

clear x
x{1}=[1 2 3 4 8 9 12 17]';
x{2}=[21 22 23 24]';
x{3}=[31 36 37 38]';

y=cellmult(2,x);

x1{1}=[1 2 3 4]';
x1{2}=[8 9 12]';
x1{3}=17;
x1{4}=[21 22 23 24]';
x1{5}=31;
x1{6}=[36 37 38]';

y1=cellmult(2,x1);

cellsplit(x,y,3);

reporttest('CELLSPLIT',aresame(x1,x)&&aresame(y1,y))


