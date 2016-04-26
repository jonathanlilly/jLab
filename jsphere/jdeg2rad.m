function[varargout]=jdeg2rad(varargin)
%JDEG2RAD  Converts degrees to radians.
%
%   [R1,R2,...,RN]=JDEG2RAD(D1,D2,...,DN) converts the input angles from
%   degrees to radians.  Output angles are in the range [-pi,pi), that
%   is, +/-180 degrees is defined to corresponds to -pi radians.
%
%   NANs and INFs in the input arguments are preserved.
%
%   See also JRAD2DEG, DEG180, DEG360, DEGUNWRAP.
% 
%   'jdeg2rad --t' runs a test.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    jdeg2rad_test,return
end

c=2*pi/360;

varargout=varargin;
for i=1:nargin;
    theta=angle(exp(sqrt(-1)*varargin{i}.*c));
    index=find(theta==pi);
    if ~isempty(index)
        theta(index)=-pi;
    end
    theta(~isfinite(theta))=varargin{i}(~isfinite(theta));      
    varargout{i}=theta;
end

 
function[]=jdeg2rad_test
th  =[0 90   180 -90   0    360+90 inf nan];
th2= [0 pi/2 -pi -pi/2 0 pi/2 inf nan];
reporttest('JDEG2RAD simple',aresame(jdeg2rad(th),th2,1e-10))
 
