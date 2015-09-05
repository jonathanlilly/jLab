function[varargout]=deg360(varargin)
%DEG360  Converts degrees to the range [0, 360].
%
%   [TH1,TH2,...THN]=DEG360(TH1,TH2,...THN) converts the input
%   angles, which are measured in degrees, to the range [0, 360].
%
%   More precisely, the output range is defined to include 0, but exclude
%   360, so DEG360(360) is defined as 360-EPS where EPS is 1e-10.
%
%   See also DEG180.
%
%   'deg360 --t' runs a test.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2013 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    deg360_test,return
end
 
varargout=varargin;
for i=1:nargin;
    temp=varargin{i};
    temp(temp==360)=(360-1e-10); %Correct for 360 
    varargout{i}=mod(temp,360);
    %temp(temp==-180|temp==)=temp
    
    %th=jrad2deg(jdeg2rad(th));
    %index=find(th<0);
    %if ~isempty(index)
    %    th(index)=th(index)+360;
    %end
    %varargout{i}=th;
    bool=~isfinite(varargin{i});
    varargout{i}(bool)=varargin{i}(bool);
end



function[]=deg360_test
thi=[-1 -179 nan inf];
tho=[359 181 nan inf];
tol=1e-10;
reporttest('DEG360 simple',aresame(deg360(thi),tho,tol))

th1=360*rand(100,1);
th2=deg360(deg180(th1));
reporttest('DEG360 inverts DEG180',aresame(th1,th2,tol))

 

%reporttest(' DEG360 ',aresame())
