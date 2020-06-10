function[varargout]=polysmooth_sortfield(varargin)
%POLYSMOOTH_SORTFIELD  Returns sorted field values for a mapping problem.
%
%   ZS=POLYSMOOTH_SORTFIELD(INDEXS,Z), where INDEXS is a cell array of
%   indices as output by SPHERESORT or TWODSORT, applies INDEXS to the
%   input field Z, leading to the output ZS.
%
%   ZS contains the elements of Z, sorted in terms of increasing distance 
%   from a set of grid points, as specified by SPHERESORT or TWODSORT.
%
%   [ZS1,ZS2,...,ZSK]=POLYSMOOTH_SORTFIELD(INDEXS,Z1,Z2,...,ZK) also works.
%
%   This is mainly an auxillary function used with SPHERESORT or TWODSORT.
%   However, it is useful to call directly for the "One grid, many fields" 
%   case described in POLYSMOOTH.
%
%   Usage: zs=polysmooth_sortfield(indexs,z);
%          [zs1,zs2,zs3]=polysmooth_sortfield(indexs,z1,z2,z3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2018--2020 J.M. Lilly --- type 'help jlab_license' for details
 
indexs=varargin{1};
varargin=varargin(2:end);

varargout=cell(length(varargin),1);
for k=1:length(varargin)
    temp=indexs;
    for j=1:length(indexs)
        nonnani=~isnan(indexs{j});
        temp{j}(nonnani)=varargin{k}(indexs{j}(nonnani));
        %adjustment for missing data in complex-valued fields
        if ~isreal(varargin{k})
            temp{j}(~nonnani)=temp{j}(~nonnani)+1i*nan;   
        end
    end
    varargout{k}=temp;
end
