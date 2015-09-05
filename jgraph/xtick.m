function[]=xtick(varargin)
%XTICK  Sets locations of x-axis tick marks.
%
%   XTICK(DX) or XTICK DX where DX is a number uses DX as the interval
%   between successive x-tick marks, with the endpoints being the axes 
%   limits, and applies this to the current axis.
%
%   XTICK(X) where X is an array sets the 'xtick' property of the current
%   axis to X.
%
%   XTICK with no input arguments attempts to make an educated guess for
%   setting 'nice' tickmark locations, usually with good results.
%
%   XTICK(H) or XTICK(H,X) or XTICK(H,DX) applies the change to axes H. 
% 
%   See also YTICK, ZTICK.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2014 J.M. Lilly --- type 'help jlab_license' for details  

h=gca;
dx=[];

if length(varargin)>1
    if ishandle(varargin{1})
        h=varargin{1};
        varargin=varargin(2:end);
    end
end

if length(varargin)>=1
   dx=varargin{1};
end

x=get(h,'xlim');
if isempty(dx)
    dxguess=[1/12 1/6 1/2 1 2 4 5 10 20 25 50 100];
    index=find(ceil((x(2)-x(1))./dxguess)<10,1,'first');
    dx=dxguess(index);
end   

if ischar(dx)
  dx=str2double(dx);
end

if length(dx)==1
  if x(2)>0&&x(1)<0
      xt=(0:dx:x(2));
      xt=[-fliplr(xt(2:end)) xt];
  else
      xt=(x(1):dx:x(2));
  end
  set(gca,'xtick',xt)
else
  set(gca,'xtick',dx);
end




