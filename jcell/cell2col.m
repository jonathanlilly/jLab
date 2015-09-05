function[varargout]=cell2col(varargin)
%CELL2COL  Converts cell arrays of column vectors into 'column-appended' data.
%
%   COL=CELL2COL(X) concatenates X, a cell array of N column vectors, into a 
%   single column vector COL with N blocks of data, each ending with a NAN.
%
%   If the data is complex, each column will end with NAN+SQRT(-1)*NAN;
%
%   To mark missing data, any pre-exiting NANs in X will be replaced with 
%   INFs, and any complex NANs will be replaced with complex INFs.
%
%   CELL2COL(X) where X is a cell array of numeric arrays of any dimension
%   also works.  In this case each array will be convered into a column 
%   vector via X{1}=X{1}(:) and so forth prior to concatenation.
%   __________________________________________________________________
%
%   Multiple input /output arguments
%
%   [C1,C2,...,CN]=CELL2COL(X1,X2,...,XN) also works.  
%   __________________________________________________________________
% 
%   Variable overwriting
%
%   CELL2COL(X1,X2,...,XN) with no output arguments overwrites the input
%   variables.
%   __________________________________________________________________
%
%   Invertibility
%  
%   CELL2COL is inverted by COL2CELL, provided 
%
%       (i)  The cell arrays XN are all column vectors of the same size and
%       (ii) The first input argument X1 contains no NANs.
%   __________________________________________________________________
%
%   See also COL2CELL, COL2MAT, MAT2COL, COLBREAKS, COL2CELL.
%
%   Usage: col=cell2col(x);
%          [c1,c2,...,cN]=cell2col(x1,x2,...,xN);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2015 J.M. Lilly --- type 'help jlab_license' for details

if strcmpi(varargin{1}, '--t')
    cell2col_test,return
end

for i=1:nargin
    ray=varargin{i}(:);  
    %if 1
    for j=1:length(ray)
        ray{j}=[ray{j};-inf];
    end
    ray=cell2mat(ray);
    
    ray=vswap(ray,nan,inf);
    ray=vswap(ray,nan+sqrt(-1)*nan,inf+sqrt(-1)*inf);
    ray=vswap(ray,-inf,nan);

    if i==1
        ray1=ray;
    elseif i~=1
        ray(isnan(ray1))=nan;
        %ray(isinf(ray1))=nan;
    end
    if ~isreal(ray)
        ray=vswap(ray,nan,nan+sqrt(-1)*nan);
    end
        
    varargout{i}=ray;
end   
eval(to_overwrite(nargin));

function[]=cell2col_test
%function[]=cell2col_test
x{1}=[1 1]';
x{2}=3;
x{3}=[2 2 2]';
y=x;
z=[1 1 3 2 2 2]';
cell2col(x,y);
reporttest('COL2CELL column vector', aresame(nonnan(x),z) && aresame(nonnan(y),z))

clear x
x{1}=[1 1]';
x{2}=3;
x{3}=[];
x{4}=[2 2 2]';
y=x;
z=[1 1 3 2 2 2]';
cell2col(x,y);
reporttest('COL2CELL column vector with empty', aresame(nonnan(x),z) && aresame(nonnan(y),z))

clear x
x{1}=[];
x{2}=[];
x{3}=[];
x{4}=[2 2 2]';
y=x;
z=[2 2 2]';
cell2col(x,y);
reporttest('COL2CELL column vector leading empties', aresame(nonnan(x),z) && aresame(nonnan(y),z))

% clear x
% x{1}=[1 1; 2 2];
% x{2}=[3 3];
% x{3}=[2 2; 3 3 ; 4 4];
% y=x;
% z=[x{1};x{2};x{3}];
% cell2col(x);
% 
% reporttest('COL2CELL two-column array', aresame(x,z))

%Further tests for cell2col are contained in col2cell