function[varargout]=deg360(varargin)
%DEG360  Converts degrees to the range [0, 360].
%
%   [D1,D2,...,DN]=DEG360(D1,D2,...,DN) converts the input angles, which 
%   are measured in degrees, to the range [0, 360].
%
%   More precisely, the output range is defined to include 0, but exclude
%   360, so DEG360(360) is defined as 360-EPS where EPS is 1e-10.
%
%   DEG360 also works if D1,D2,...,DN are cell arrays of numerical arrays.
%
%   See also JDEG2RAD, JRAD2DEG, DEG360, DEGUNWRAP.
%
%   'deg360 --t' runs a test.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    deg360_test,return
end

varargout=varargin;
for i=1:nargin;
    if ~iscell(varargin{1})
        temp=varargin{i};
        temp(temp==360)=(360-1e-10); %Correct for 360
        varargout{i}=mod(temp,360);
        bool=~isfinite(varargin{i});
        varargout{i}(bool)=varargin{i}(bool);
    else
        for j=1:length(varargin{1})
            temp=varargin{i}{j};
            temp(temp==360)=(360-1e-10); %Correct for 360
            varargout{i}{j,1}=mod(temp,360);
            bool=~isfinite(varargin{i}{j});
            varargout{i}{j}(bool)=varargin{i}{j}(bool);
        end
    end
    %temp(temp==-180|temp==)=temp
    
    %th=jrad2deg(jdeg2rad(th));
    %index=find(th<0);
    %if ~isempty(index)
    %    th(index)=th(index)+360;
    %end
    %varargout{i}=th;

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
