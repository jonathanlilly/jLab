function[varargout]=cellmean(varargin)
%CELLMEAN  Mean value of each element a cell array.
%
%   M=CELLMEAN(X) where X is a cell array of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   returns the N x 1 array of mean values M with
%  
%      M(1)=MEAN(X1(:)),  M(2)=MEAN(X2(:)),...,  M(N)=MEAN(XN(:)).
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
%   Parallelization
%
%   CELLMEAN(...,'parallel') parallelizes the computation using a PARFOR 
%   loop.  This requires that Matlab's Parallel Computing Toolbox be 
%   installed, and is useful for very large datasets.
%   __________________________________________________________________   
%
%   See also JCELL, VTOOLS.
%
%   Usage: m=cellmean(x);
%          [m1,m2,m3]=cellmean(x1,x2,x3);
%          cellmean(x1,x2,x3);
%          cellmean(x1,x2,x3,'parallel');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2015 J.M. Lilly --- type 'help jlab_license' for details

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
    varargout{i}=cellmean_one(varargin{i},cores);
end

eval(to_overwrite(length(varargin)))


function[y]=cellmean_one(x,cores)
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

