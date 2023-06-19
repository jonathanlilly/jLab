function[varargout]=pm_index(varargin)
%PM_INDEX  Applies an index output by PM_SORT to sort data for POLYMAP.
%
%   PM_INDEX is called internally by POLYMAP.  However, for large problems
%   it may be preferable to call it externally, as documented in POLYMAP.
%
%   ZS=PM_INDEX(INDEXS,Z), where INDEXS is a cell array of indices as 
%   output by PM_SORT, applies INDEXS to the input field Z.
%
%   ZS then contains the elements of Z, arranged on the grid that was input 
%   to PM_SORT and by terms of increasing distance from the grid points.
%
%   If Z is empty, PM_SORT returns ZS as a cell array of empty arrays.
%
%   [ZS1,ZS2,...,ZSK]=PM_SORT(INDEXS,Z1,Z2,...,ZK) also works.
%
%   See also: POLYMAP and PM_SORT.
%
%   Usage: zs=pm_index(indexs,z);
%          [zs1,zs2,zs3]=pm_index(indexs,z1,z2,z3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2022 J.M. Lilly --- type 'help jlab_license' for details
 
indexs=varargin{1};
if ischar(varargin{end})
    parstr=varargin{end};
    varargin=varargin(1:end-1);
else
    parstr='ser';
end

varargin=varargin(2:end);
varargout=cell(length(varargin),1);

for k=1:length(varargin)
    if ~isempty(varargin{k})
        temp=indexs;
        if strcmpi(parstr(1:3),'par')
            parfor j=1:length(indexs)
                nonnani=~isnan(indexs{j});
                temp{j}(nonnani)=varargin{k}(indexs{j}(nonnani));
                %adjustment for missing data in complex-valued fields
                if ~isreal(varargin{k})
                    temp{j}(~nonnani)=temp{j}(~nonnani)+1i*nan;
                end
            end
        else
            for j=1:length(indexs)
                nonnani=~isnan(indexs{j});
                temp{j}(nonnani)=varargin{k}(indexs{j}(nonnani));
                %adjustment for missing data in complex-valued fields
                if ~isreal(varargin{k})
                    temp{j}(~nonnani)=temp{j}(~nonnani)+1i*nan;
                end
            end
        end
        varargout{k}=temp;
    else
        varargout{k}=cell(size(indexs));
    end
end