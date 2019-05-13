function[varargout]=cellmean(varargin)
%CELLMEAN  Mean value of each element a cell array, possibly weighted.
%
%   M=CELLMEAN(X) where X is a cell array of N arrays, is equivalent to
%  
%     M(1,1)=MEAN(X{1}(:)), M(2,1)=MEAN(X{2}(:)), ..., M(N,1)=MEAN(X{N}(:))
%
%   thus returning an N x 1 array containing the mean of all values in each
%   element in the cell array.
%
%   In taking the mean, non-finite values are ignored, as in VMEAN.
%   M is a column vector of the same length as X.  
%
%   [M1,M2,...,MP]=CELLMEAN(X1,X2,...,XP) also works for P different 
%   input arguments. 
%
%   CELLMEAN(X1,X2,...XP);  with no output arguments overwrites the 
%   original input variables.
%   __________________________________________________________________   
%
%   Weighted means
%
%   CELLMEAN(X,'weight',W) or CELLMEAN(X1,X2,...,XP,'weight',W) where W is 
%   a cell array of the same size as the other input variables, computes
%   the weighted mean, using the weighting factor ABS(W).^2.  
%   __________________________________________________________________   
%
%   Parallelization
%
%   CELLMEAN(...,'parallel') parallelizes the computation using a PARFOR 
%   loop.  This requires that Matlab's Parallel Computing Toolbox be 
%   installed, and is useful for very large datasets.
%   __________________________________________________________________   
%
%   See also CELLSTD, CELLMED, JCELL.
%
%   Usage: m=cellmean(x);
%          [m1,m2,m3]=cellmean(x1,x2,x3);
%          cellmean(x1,x2,x3,'weight',w);
%          cellmean(x1,x2,x3);
%          cellmean(x1,x2,x3,'parallel');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2019 J.M. Lilly --- type 'help jlab_license' for details

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
        varargout{i}=cellmean_one_unweighted(varargin{i},cores);
    else
        varargout{i}=cellmean_one_weighted(varargin{i},weight,cores);
    end
end

eval(to_overwrite(length(varargin)))


function[y]=cellmean_one_unweighted(x,cores)

y=nan*zeros(length(x),1);
if strcmpi(cores(1:3),'par')
    parfor i=1:length(x)
        y(i)=mean(x{i}(:),1);
        if ~isfinite(y(i))
            y(i)=vmean(x{i}(:),1);
        end
    end
else
    for i=1:length(x)
        y(i)=mean(x{i}(:),1);
        if ~isfinite(y(i))
            y(i)=vmean(x{i}(:),1);
        end
    end
end


function[y]=cellmean_one_weighted(x,weight,cores)

y=nan*zeros(length(x),1);
if strcmpi(cores(1:3),'par')
    parfor i=1:length(x)
        w=squared(weight{i}(:));
        y(i)=mean(w.*x{i}(:),1)./mean(w,1);
        if ~isfinite(y(i))
            y(i)=vmean(w.*x{i}(:),1)./vmean(w,1);
        end
    end
else
    for i=1:length(x)
        w=squared(weight{i}(:));
        y(i)=mean(w.*x{i}(:),1)./mean(w,1);
        if ~isfinite(y(i))
            y(i)=vmean(w.*x{i}(:),1)./vmean(w,1);
        end
    end
end
% y=vmean(col2mat(cell2col(x)),1);
% 
% if length(y)~=length(x)
%     disp('Fast version of CELLMEAN failed, perhaps because X contained NANs.')
%     disp('Be patient, implementing slow version.')
%     y=nan*zeros(length(x),1);
%     for i=1:length(x)
%         y(i)=vmean(vcolon(x{i}),1);
%     end
% end

