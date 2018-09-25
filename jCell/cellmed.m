function[varargout]=cellmed(varargin)
%CELLMED  Median value of each element a cell array.
%
%   M=CELLMED(X) where X is a cell array of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   returns the N x 1 array of median values M with
%  
%      M(1)=MEDIAN(X1(:)),  M(2)=MEDIAN(X2(:)),...,  M(N)=MEDIAN(XN(:)).
%
%   In taking the median, non-finite values are ignored, as in VMEDIAN.
%   M is a column vector of the same length as X.  
%
%   [M1,M2,...,MP]=CELLMED(X1,X2,...,XP) also works for P different 
%   input arguments. 
%
%   CELLMED(X1,X2,...XP);  with no output arguments overwrites the 
%   original input variables.
%   __________________________________________________________________   
%
%   Parallelization
%
%   CELLMED(...,'parallel') parallelizes the computation using a PARFOR 
%   loop.  This requires that Matlab's Parallel Computing Toolbox be 
%   installed, and is useful for very large datasets.
%   __________________________________________________________________   
%
%   See also CELLMEAN, CELLSTD, JCELL, VTOOLS.
%
%   Usage: m=cellmed(x);
%          [m1,m2,m3]=cellmed(x1,x2,x3);
%          cellmed(x1,x2,x3);
%          cellmed(x1,x2,x3,'parallel');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2018 J.M. Lilly --- type 'help jlab_license' for details

if ~iscell(varargin{1})
    error('X must be a cell array.')
end

cores='serial';

if ischar(varargin{end})
    if strcmpi(varargin{end}(1:3),'par')||strcmpi(varargin{end}(1:3),'ser')
        cores=varargin{end};
    end
    varargin=varargin(1:end-1);
end

if strcmpi(cores(1:3),'par')
    if exist('parpool')~=2
        disp('Sorry, parallel algorithm requires the Parallel Computing Toolbox.')
        disp('Defaulting to the standard algorithm.')
        cores='serial';
    end
end

for i=1:length(varargin)
    varargout{i}=cellmed_one(varargin{i},cores);
end

eval(to_overwrite(length(varargin)))

function[y]=cellmed_one(x,cores)
y=nan*zeros(length(x),1);
if strcmpi(cores(1:3),'par')
    parfor i=1:length(x)
        y(i)=median(x{i}(:),1);
        if ~isfinite(y(i))
            y(i)=vmedian(x{i}(:),1);
        end
    end
else
    for i=1:length(x)
        y(i)=median(x{i}(:),1);
        if ~isfinite(y(i))
            y(i)=vmedian(x{i}(:),1);
        end
    end
end
