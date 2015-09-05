function[varargout]=deg180(varargin)
%DEG180  Converts degrees to the range [-180,180].
%
%   [TH1,TH2,...THN]=DEG180(TH1,TH2,...THN) converts the input
%   angles, which are measured in degrees, to the range [-180, 180].
%
%   See also DEG360.
%
%   'deg180 --t' runs a test.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007-2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    deg180_test,return
end
 
varargout=varargin;
for i=1:nargin;
    varargout{i}=deg360(varargin{i});
    bool=varargout{i}>180;
    varargout{i}(bool)=varargout{i}(bool)-360;
    %varargout{i}=angle(exp(sqrt(-1)*varargin{i}*2*pi/360))*360/2/pi;
    %varargout{i}(~isfinite(varargin{i}))=varargin{i}(~isfinite(varargin{i}));
%    varargout{i}=jrad2deg(jdeg2rad(varargin{i})); Same but slower
end

function[]=deg180_test
thi=[359 181 nan inf];
tho=[-1 -179 nan inf];
tol=1e-10;
reporttest('DEG180 simple',aresame(deg180(thi),tho,tol))
