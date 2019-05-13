function[varargout]=cellstd(varargin)
%CELLSTD  Standard deviation of each element a cell array, possibly weighted.
%
%   S=CELLSTD(X) where X is a cell array of N arrays, is equivalent to
%  
%      S(1,1)=STD(X{1}(:)),  S(2,1)=STD(X{2}(:)), ...,  S(N,1)=STD(X{N}(:))
%
%   thus returning an N x 1 array containing the standard deviation of all
%   values in each element in the cell array.
%
%   In taking the standard deviation, non-finite values are ignored, as in
%   VSTD. The output S is a column vector of the same length as X.  
%
%   [S1,S2,...,SP]=CELLSTD(X1,X2,...,XP) also works for P different input
%   arguments. 
%
%   CELLSTD(X1,X2,...XP);  with no output arguments overwrites the 
%   original input variables.
%   __________________________________________________________________   
%
%   Weighted standard deviation
%
%   CELLSTD(X,'weight',W) or CELLSTD(X1,X2,...,XP,'weight',W) where W is 
%   a cell array of the same size as the other input variables, computes
%   the weighted standard deviation, using the weighting factor ABS(W).^2.  
%   __________________________________________________________________   
%
%   Parallelization
%
%   CELLSTD(...,'parallel') parallelizes the computation using a PARFOR 
%   loop.  This requires that Matlab's Parallel Computing Toolbox be 
%   installed, and is useful for very large datasets.
%   __________________________________________________________________   
%
%   See also CELLMEAN, CELLMED, JCELL.
%
%   Usage: s=cellstd(x);
%          [s1,s2,s3]=cellstd(x1,x2,x3);
%          [s1,s2,s3]=cellstd(x1,x2,x3,'weight',w);
%          cellstd(x1,x2,x3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015--2019 J.M. Lilly --- type 'help jlab_license' for details

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
        varargout{i}=cellstd_one_unweighted(varargin{i},cores);
    else
        varargout{i}=cellstd_one_weighted(varargin{i},weight,cores);
    end
end


eval(to_overwrite(length(varargin)))

function[y]=cellstd_one_unweighted(x,cores)

y=nan*zeros(length(x),1);
if strcmpi(cores(1:3),'par')
    parfor i=1:length(x)
        y(i)=std(x{i}(:),1);
        if ~isfinite(y(i))
            y(i)=vstd(x{i}(:),1);
        end
    end
else
    for i=1:length(x)
        y(i)=std(x{i}(:),1);
        if ~isfinite(y(i))
            y(i)=vstd(x{i}(:),1);
        end
    end
end


function[y]=cellstd_one_weighted(x,weight,cores)

y=nan*zeros(length(x),1);
if strcmpi(cores(1:3),'par')
    parfor i=1:length(x)
        w=squared(weight{i}(:));
        m=mean(w.*x{i}(:),1)./mean(w,1);
        y(i)=sqrt(mean(w.*squared(x{i}(:)-m),1)./mean(w,1));
        if ~isfinite(y(i))
            m=vmean(w.*x{i}(:),1)./vmean(w,1);
            y(i)=sqrt(vmean(w.*squared(x{i}(:)-m),1)./vmean(w,1));
        end
    end
else
    for i=1:length(x)
        w=squared(weight{i}(:));
        m=mean(w.*x{i}(:),1)./mean(w,1);
        y(i)=sqrt(mean(w.*squared(x{i}(:)-m),1)./mean(w,1));
        if ~isfinite(y(i))
            m=vmean(w.*x{i}(:),1)./vmean(w,1);
            y(i)=sqrt(vmean(w.*squared(x{i}(:)-m),1)./vmean(w,1));
        end
    end
end
