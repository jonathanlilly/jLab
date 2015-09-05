function[]=ytick(varargin)
%YTICK  Sets locations of y-axis tick marks.
%
%   YTICK(DY) or YTICK DY where DY is a number uses DY as the interval
%   between successive y-tick marks, with the endpoints being the axes 
%   limits, and applies this to the current axis.
%
%   YTICK(Y) where Y is an array sets the 'ytick' property of the current
%   axis to Y.
%
%   YTICK with no input arguments attempts to make an educated guess for
%   setting 'nice' tickmark locations, usually with good results.
%
%   YTICK(H) or YTICK(H,Y) or YTICK(H,DY) applies the change to axes H.
% 
%   See also XTICK, ZTICK.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2014 J.M. Lilly --- type 'help jlab_license' for details   
    
h=gca;
dy=[];

if length(varargin)>1
    if ishandle(varargin{1})
        h=varargin{1};
        varargin=varargin(2:end);
    end
end

if length(varargin)>=1
   dy=varargin{1};
end

y=get(h,'ylim');
if isempty(dy)
    dyguess=[1/12 1/6 1/2 1 2 4 5 10 20 25 50 100];
    index=find(ceil((y(2)-y(1))./dyguess)<10,1,'first');
    dy=dyguess(index);
end   

if ischar(dy)
  dy=str2double(dy);
end

if length(dy)==1
  if y(2)>0&&y(1)<0
      yt=(0:dy:y(2));
      yt=[-fliplr(yt(2:end)) yt];
  else
      yt=(y(1):dy:y(2));
  end
  set(gca,'ytick',yt)
else
  set(gca,'ytick',dy);
end

