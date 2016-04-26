function[varargout]=cellstd(varargin)
%CELLSTD  Standard deviation of each element a cell array.
%
%   SIG=CELLSTD(X) where X is a cell array of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   returns the N x 1 array of standard deviation values SIG such that
%  
%      SIG(1)=STD(X1(:)),  SIG(2)=STD(X2(:)),...,  SIG(N)=STD(XN(:)).
%
%   In taking the standard deviation, non-finite values are ignored, as in
%   VSTD. The output SIG is a column vector of the same length as X.  
%
%   [SIG1,SIG2,...,SIGP]=CELLSTD(X1,X2,...,XP) also works for P different 
%   input arguments. 
%
%   CELLSTD(X1,X2,...XP);  with no output arguments overwrites the 
%   original input variables.
%   __________________________________________________________________   
%
%   Parallelization
%
%   CELLSTD(...,'parallel') parallelizes the computation using a PARFOR 
%   loop.  This requires that Matlab's Parallel Computing Toolbox be 
%   installed, and is useful for very large datasets.
%   __________________________________________________________________   
%
%   See also JCELL, VTOOLS.
%
%   Usage: sig=cellstd(x);
%          [sig1,sig2,sig3]=cellstd(x1,x2,x3);
%          cellstd(x1,x2,x3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details

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
    varargout{i}=cellstd_one(varargin{i},cores);
end

eval(to_overwrite(length(varargin)))

function[y]=cellstd_one(x,cores)
y=nan*zeros(length(x),1);
if strcmpi(cores(1:3),'par')
    parfor i=1:length(x)
        y(i)=std(x{i}(:),1,1);
        if ~isfinite(y(i))
            y(i)=vstd(x{i}(:),1);
        end
    end
else
    for i=1:length(x)
        y(i)=std(x{i}(:),1,1);
        if ~isfinite(y(i))
            y(i)=vstd(x{i}(:),1);
        end
    end
end

