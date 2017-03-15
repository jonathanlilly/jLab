function[varargout]=linecolor(varargin)
%LINECOLOR  Set line colors based on a property value within a colormap.
%
%   LINECOLOR(H,C) set the lines with handles H to the colors C, determined
%   by looking up the values of C within the 'lansey' colormap.  
%
%   H and C should be arrays of the same size.
%
%   LINECOLOR(H,C,CMIN,CMAX) uses the values CMIN and CMAX as the lower and
%   uppermost values in the colormap, respectively.  The default behavior 
%   is to set CMIN and CMAX to the minimum and maximum values within C.
%
%   LINECOLOR(...,MAP) alternately uses the colormap with the name MAP.
%
%   Usage: linecolor(h,c);
%          linecolor(h,c,map);
%          linecolor(h,c,cmin,cmax);
%          linecolor(h,c,cmin,cmax,map);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2016 J.M. Lilly --- type 'help jlab_license' for details

h=varargin{1}(:);
c=varargin{2}(:);
c=vswap(c,nan,inf);

if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
else
    str='lansey';
end
varargin=varargin(3:end);

if length(varargin)==2
    cmin=varargin{1};
    cmax=varargin{2};
else
    cmin=minmin(c);
    cmax=maxmax(c);
end

map=colormap(str);
N=size(map,1);

m=frac(N-1,cmax-cmin);
b=1-m*cmin;
y=ceil(m*c+b);

y(y>N)=N;
y(y<1)=1;

set(h,'visible','off')
%get(h(1))
%y

for i=1:length(h),
    set(h(i),'color',map(y(i),:));
end

set(h,'visible','on')

