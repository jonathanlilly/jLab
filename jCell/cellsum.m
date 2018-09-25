function[varargout]=cellsum(varargin)
%CELLSUM  Mean value of each element a cell array, possibly weighted.
%
%   M=CELLSUM(X) where X is a cell array of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   returns the N x 1 array of mean values M with
%  
%      M(1)=SUM(X1(:)),  M(2)=SUM(X2(:)),...,  M(N)=SUM(XN(:)).
%
%   In taking the mean, non-finite values are ignored, as in VSUM.
%   M is a column vector of the same length as X.  
%
%   [M1,M2,...,MP]=CELLSUM(X1,X2,...,XP) also works for P different 
%   input arguments. 
%
%   CELLSUM(X1,X2,...XP);  with no output arguments overwrites the 
%   original input variables.
%   __________________________________________________________________   
%
%   Weighted sums
%
%   CELLSUM(X,'weight',W) or CELLSUM(X1,X2,...,XP,'weight',W) where W is 
%   a cell array of the same size as the other input variables, computes
%   the weighted sum, using the weighting factor ABS(W).^2.  
%   __________________________________________________________________   
%
%   Parallelization
%
%   CELLSUM(...,'parallel') parallelizes the computation using a PARFOR 
%   loop.  This requires that Matlab's Parallel Computing Toolbox be 
%   installed, and is useful for very large datasets.
%   __________________________________________________________________   
%
%   See also CELLSTD, CELLMED, JCELL, VTOOLS.
%
%   Usage: m=cellsum(x);
%          [m1,m2,m3]=cellsum(x1,x2,x3);
%          cellsum(x1,x2,x3,'weight',w);
%          cellsum(x1,x2,x3);
%          cellsum(x1,x2,x3,'parallel');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2018 J.M. Lilly --- type 'help jlab_license' for details

if ~iscell(varargin{1})
    error('X must be a cell array.')
end

weight=[];
cores='serial';

if ischar(varargin{end})
    if strcmpi(varargin{end}(1:3),'par')||strcmpi(varargin{end}(1:3),'ser')
        cores=varargin{end};
    end
    varargin=varargin(1:end-1);
end

if length(varargin)>1
    if ischar(varargin{end-1})
        if strcmpi(varargin{end-1}(1:3),'wei')
            weight=varargin{end};
        end
        varargin=varargin(1:end-2);
    end
end

if strcmpi(cores(1:3),'par')
    if exist('parpool')~=2
        disp('Sorry, parallel algorithm requires the Parallel Computing Toolbox.')
        disp('Defaulting to the standard algorithm.')
        cores='serial';
    end
end

for i=1:length(varargin)
    if isempty(weight)
        varargout{i}=cellsum_one_unweighted(varargin{i},cores);
    else
        varargout{i}=cellsum_one_weighted(varargin{i},weight,cores);
    end
end

eval(to_overwrite(length(varargin)))


function[y]=cellsum_one_unweighted(x,cores)

y=nan*zeros(length(x),1);
if strcmpi(cores(1:3),'par')
    parfor i=1:length(x)
        y(i)=sum(x{i}(:),1);
        if ~isfinite(y(i))
            y(i)=vsum(x{i}(:),1);
        end
    end
else
    for i=1:length(x)
        y(i)=sum(x{i}(:),1);
        if ~isfinite(y(i))
            y(i)=vsum(x{i}(:),1);
        end
    end
end


function[y]=cellsum_one_weighted(x,weight,cores)

y=nan*zeros(length(x),1);
if strcmpi(cores(1:3),'par')
    parfor i=1:length(x)
        w=squared(weight{i}(:));
        y(i)=sum(w.*x{i}(:),1);
        if ~isfinite(y(i))
            y(i)=vsum(w.*x{i}(:),1);
        end
    end
else
    for i=1:length(x)
        w=squared(weight{i}(:));
        y(i)=sum(w.*x{i}(:),1);
        if ~isfinite(y(i))
            y(i)=vsum(w.*x{i}(:),1);
        end
    end
end

