function[varargout]=linecolor(varargin)
%LINECOLOR  Set line colors based on a property value within a colormap.
%
%   LINECOLOR(H,C) set the lines with handles H to the colors C, determined
%   by looking up the values of C within the default figure colormap. 
%
%   H and C should be arrays of the same size.
%
%   LINECOLOR(H,C,CMIN,CMAX) uses the values CMIN and CMAX as the lower and
%   uppermost values in the colormap, respectively.  The default behavior 
%   is to set CMIN and CMAX to the minimum and maximum values within C.
%
%   LINECOLOR(...,MAP) alternately uses the colormap MAP, which can either
%   be a string or a colormap matrix.
%
%   Usage: linecolor(h,c);
%          linecolor(h,c,map);
%          linecolor(h,c,cmin,cmax);
%          linecolor(h,c,cmin,cmax,map);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2021 J.M. Lilly --- type 'help jlab_license' for details

h=varargin{1}(:);
c=varargin{2}(:);
c=vswap(c,nan,inf);

map0=colormap;

if nargin==3||nargin==5
    if ischar(varargin{end})
        str=varargin{end};
        map=colormap(str);
    else
        map=varargin{end};
    end
    varargin=varargin(1:end-1);
else
    map=get(0,'DefaultFigureColormap');
end

if length(varargin)==4
    cmin=varargin{3};
    cmax=varargin{4};
else
    cmin=minmin(c);
    cmax=maxmax(c);
end

N=size(map,1);

%[cmin,cmax]

m=frac(N-1,cmax-cmin);
b=1-m*cmin;
y=ceil(m*c+b);

y(y>N)=N;
y(y<1)=1;

set(h,'visible','off')
%get(h(1))
%y
%hold on,plot(y)

for i=1:length(h)
    set(h(i),'color',map(y(i),:));
end

colormap(map0);
set(h,'visible','on')

