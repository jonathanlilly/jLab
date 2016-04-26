function[varargout]=jrad2deg(varargin)
%JRAD2DEG  Converts radians to degrees.
%
%   [D1,D2,...,DN]=JRAD2DEG(R1,R2,...,RN) converts the input angles from
%   radians to degrees.  Output angles are in the range [-180,180), 
%   that is, +/- pi radians is defined to correspond to -180 degrees.
%
%   NANs and INFs in the input arguments are preserved.
%
%   See also JDEG2RAD, DEG180, DEG360, DEGUNWRAP.
% 
%   'jrad2deg --t' runs a test.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    jrad2deg_test,return
end
 
c=2*pi/360;

varargout=varargin;
for i=1:nargin;
    theta=angle(exp(sqrt(-1)*varargin{i}))./c;
    index=find(theta==180);
    if ~isempty(index)
        theta(index)=-180;
    end
    theta(~isfinite(theta))=varargin{i}(~isfinite(theta));      
    varargout{i}=theta;
end

function[]=jrad2deg_test
 
th= [0 pi/2 pi 3*pi/2 2*pi 5*pi/2 inf nan];
th2=[0 90   -180 -90   0    90 inf nan];
reporttest('JRAD2DEG simple',aresame(jrad2deg(th),th2,1e-10))

th=360*rand(100,1)-180;
reporttest('JRAD2DEG inverts JDEG2RAD',aresame(jrad2deg(jdeg2rad(th)),th,1e-10))
