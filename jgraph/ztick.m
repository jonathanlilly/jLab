function[]=ztick(varargin)
%ZTICK  Sets locations of z-axis tick marks.
%
%   ZTICK(DZ) or ZTICK DX where DX is a number uses DZ as the interval
%   between successive z-tick marks, with the endpoints being the axes 
%   limits, and applies this to the current axis.
%
%   ZTICK(Z) where Z is an array sets the 'ztick' property of the current 
%   axis to Z.
%
%   ZTICK with no input arguments attempts to make an educated guess for
%   setting 'nice' tickmark locations, usually with good results.
%
%   ZTICK(H) or ZTICK(H,Z) or ZTICK(H,DZ) applies the change to axes H. 
% 
%   See also XTICK, YTICK.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2017 J.M. Lilly --- type 'help jlab_license' for details  

h=gca;
dz=[];

if length(varargin)>1
    if ishandle(varargin{1})
        h=varargin{1};
        varargin=varargin(2:end);
    end
end

if length(varargin)>=1
   dz=varargin{1};
end

z=get(h,'zlim');
if isempty(dz)
    dzguess=[1/12 1/6 1/2 1 2 4 5 10 20 25 50 100];
    index=find(ceil((z(2)-z(1))./dzguess)<10,1,'first');
    dz=dzguess(index);
end   

if ischar(dz)
  dz=str2double(dz);
end

if length(dz)==1
  if z(2)>0&&z(1)<0
      zt=(0:dz:z(2));
      zt=[-fliplr(zt(2:end)) zt];
  else
      zt=(z(1):dz:z(2));
  end
  set(gca,'ztick',zt)
else
  set(gca,'ztick',dz);
end



