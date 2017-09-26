function[varargout]=cms2kmd(varargin)
%CMS2KMD  Converts centimeters per second to kilometers per day.
%
%   Y=CMS2KMD(X)   <==>  Y=X*(3600*24/100/1000)
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2012--2015 J.M. Lilly --- type 'help jlab_license' for details


for i=1:nargin
    varargout{i}=(3600*24/100/1000).*varargin{i};
end
