function[h,hc]=jpcolor(varargin)
%JPCOLOR  Modified version of PCOLOR appropriate for cell-centered grids.
%
%   JPCOLOR(XMID,YMID,Z) makes a PCOLOR plot with XMID and YMID marking the
%   *centers* of the cells of Z.  X and Y need to be monotonic, but do not 
%   need to have uniform spacing.
%
%   This is unlike PCOLOR(X,Y,Z), where X and Y mark the cell *edges*.  
%
%   Similarly, unlike PCOLOR, JPCOLOR does not throw away the last row and
%   column of Z. It also specifies the axes to be centered on the cells.
%
%   JPCOLOR also automatically applies SQUEEZE to Z, which is useful for 
%   working with slices of multidimensional datasets.
%   
%   The default shading for JPCOLOR is FLAT rather than FACETED.  JPCOLOR
%   also sets the level of the PCOLOR image to the back of the plot, and 
%   sets the axis layer to the top so that tickmarks are not obscured.
%
%   The color axes are set to the 0.1% and 99.9% data quantiles using
%   COLORQUANT, rather than to the minimum and maximum data values. 
%
%   H=JPCOLOR(...) returns the handle H to the PCOLOR image.
%   _______________________________________________________________________
%
%   Adding a colorbar
%
%   JPCOLOR(...,LABEL), where LABEL is a string containing a colorbar
%   label, e.g. 'Temperature', adds a colorbar and labels it accordingly.
%
%   JPCOLOR(...,LABEL,LOC) puts the colorbar in location LOC.  LOC is one
%   of the eight strings 'N','S','E','W','NO','SO','EO','WO', corresponding
%   to the eight locations options discussed in COLORBAR.
%
%   JPCOLOR(...,[],LOC) adds a colorbar in location LOC with no label.
%
%   [H,HC]=JPCOLOR(...) also returns the handle HC to the colorbar.
%
%   'jpcolor --f' generates a figure showing the differences from PCOLOR.
%
%   Usage: jpcolor(z)
%          jpcolor(x,y,z);
%          jpcolor(x,y,z,'Temperature');
%          jpcolor(x,y,z,'Temperature','SO');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2020 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--f')
    type makefigs_jpcolor
    makefigs_jpcolor;
    return
end

label=[];
loc='EO';
if length(varargin)>1
    if ischar(varargin{end-1})||isempty(varargin{end-1}) 
    label=varargin{end-1};
    loc=varargin{end};
    varargin=varargin(1:end-2);
    end
end
if ischar(varargin{end})||isempty(varargin{end}) 
    label=varargin{end};
    varargin=varargin(1:end-1);
end

if length(varargin)==1
    x=varargin{1};
    z=x;
    x=[0:size(z,2)-1]+1/2;
    y=[0:size(z,1)-1]+1/2;
elseif length(varargin)>1
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
end    

if (length(find(size(x)>1))>1)||(length(find(size(y)>1))>1)
    error('Sorry, JPCOLOR does not support matrix-valued X and Y.')
end

x=x(:);
y=y(:);

dx=x(2)-x(1);
dy=y(2)-y(1);

x=interp1([1:length(x)]',x,[1:length(x)+1]','pchip','extrap');
y=interp1([1:length(y)]',y,[1:length(y)+1]','pchip','extrap');

%x=[x;x(end)+dx];
%x=x-dx/2;

%y=[y;y(end)+dy];
%y=y-dy/2;

z=squeeze(z);

z=[z z(:,end)];
z=[z; z(end,:)];
%z=vswap(z,0,nan);
z=vswap(z,-inf,nan);

%vsize(x,y,z);

%h=pcolor(x,y,z);shading flat
h=pcolor(x-dx/2,y-dy/2,z);shading flat
uistack(h,'bottom')
axis([min(x) max(x) min(y) max(y)])
boxon

colorquant

set(gca,'layer','top')
%This shows tickmarks

if ~isempty(label)
    hc=colorbar(loc);
    hc.Label.String=label;
end

if nargout==0, clear h hc, end
