function[varargout]=jpcolor(x,y,z)
%JPCOLOR  Modified version of PCOLOR appropriate for cell-centered grids.
%
%   JPCOLOR(XMID,YMID,Z) makes a PCOLOR plot with XMID and YMID marking the
%   *centers* of the cells of Z.  X and Y need to be monotonic, but do not 
%   need to have uniform spacing.
%
%   This is unlike PCOLOR(X,Y,Z), where X and Y mark the cell *edges*.  
%
%   Similarly, unlike PCOLOR, JPCOLOR does not throw away the last row and
%   column of Z.
%   
%   The default shading for JPCOLOR is FLAT rather than FACETED.  JPCOLOR
%   also sets the level of the PCOLOR image to the back of the plot.
%
%   'jpcolor --f' generates a figure showing the differences from PCOLOR.
%
%   Usage: jpcolor(x,y,z);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2016 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(x, '--f')
    type makefigs_jpcolor
    makefigs_jpcolor;
    return
end

if nargin==1
    z=x;
    x=[0:size(z,2)-1]+1/2;
    y=[0:size(z,1)-1]+1/2;
end

if (length(find(size(x)>1))>1)||(length(find(size(y)>1))>1)
    error('Sorry, JPCOLOR does not support matrix-valued X and Y.')
end

x=x(:);
y=y(:);

%dx=x(2)-x(1);
%dy=y(2)-y(1);

x=interp1([1:length(x)]',x,[1:length(x)+1]','pchip','extrap');
y=interp1([1:length(y)]',y,[1:length(y)+1]','pchip','extrap');

%x=[x;x(end)+dx];
%x=x-dx/2;

%y=[y;y(end)+dy];
%y=y-dy/2;

z=[z z(:,end)];
z=[z; z(end,:)];
%z=vswap(z,0,nan);
z=vswap(z,-inf,nan);

%vsize(x,y,z);

h=pcolor(x,y,z);shading flat
uistack(h,'bottom')
axis([min(x) max(x) min(y) max(y)])
boxon

