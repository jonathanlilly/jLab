function[varargout]=polysmooth_extract(varargin)
%POLYSMOOTH_EXTRACT
%
%   POLYSMOOTH_EXTRACT is an auxiliary function for use with POLYSMOOTH.
%
%
%
%   Usage: zs=polysmooth_extract(siz,extract,z);
%          [zs1,zs2,zs3]=polysmooth_extract(siz,extract,z1,z2,z3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2018 J.M. Lilly --- type 'help jlab_license' for details
 
siz=varargin{1};
extract=varargin{2};
varargin=varargin(3:end);

for i = 1:length(varargin)
    varargout{i}=nan*zeros(siz);
    if ~isreal(varargin{i})
            varargout{i}=varargout{i}+1*nan*zeros(siz);
    end
    varargout{i}(~isnan(extract))= varargin{i}(extract(~isnan(extract)));
end

%    zs=nan*zeros(size(xs));
%    zs(~isnan(extract))=z(extract(~isnan(extract)));