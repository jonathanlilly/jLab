function[z,tc]=trianglepath(varargin)
%TRIANGLEPATH  Moving instrument path composed of adjacent triangles.
%
%   [Z,TC]=TRIANGLEPATH(DX,LX,N) forms path consisting of a total of 2N 
%   equilateral triangles with all sides of length LX, and with DX giving 
%   the distance between successive points.  
%
%   The first cell begins at (0,0).  From there one moves southwest, east,
%   and northwest to form an upper triangle, then east to a new starting
%   point.  This process then repeats N times, with a final movement to the
%   southwest. The result is 2N triangles.
%
%   Z is the complex-valued path traced out in this way, such that PLOT(Z)
%   plots the triangles.  TC is a cell array with one triangle per cell.
%
%   This is used in MAKEFIGS_STOKES to draw a hypothetical sampling path 
%   consisting of adjacent closed cells formed by a single moving platform.
%
%   Usage: [z,tc]=trianglepath(dx,lx,N);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2024 J.M. Lilly --- type 'help jlab_license' for details
 
xo=0;
dx=varargin{1};
lx=varargin{2};
n=varargin{3};
if nargin>3
    xo=varargin{4};
end

if strcmp(dx, '--t')
    trianglepath_test,return
end
  
 
if n>0
    x=[(0:-dx:-lx)/2,(0:dx:lx)-lx/2,(0:-dx:-lx)/2+lx/2,(0:dx:lx)]';
    y=[(0:-dx:-lx),0*(0:dx:lx)-lx,(-lx:dx:0),0*(0:dx:lx)]';
    %vsize(x,y)
    z=[xo+x+1i*y;trianglepath(dx,lx,n-1,xo+lx)];
else
    z=[xo+(0:-dx:-lx)'/2+1i*(0:-dx:-lx)'];
end

if nargout>1
    %lower triangles
    x=[(0:-dx:-lx)/2,(0:dx:lx)-lx/2,(0:-dx:-lx)/2+lx/2]';
    y=[(0:-dx:-lx),0*(0:dx:lx)-lx,(-lx:dx:0)]';
    for i=1:2:2*n
        tc{i}=(x+1i*y)+(i-1)/2*lx;
    end

    %upper triangles
    x=[(0:-dx:-lx)/2+lx/2,(0:dx:lx),lx+(0:-dx:-lx)/2]';
    y=[(-lx:dx:0),0*(0:dx:lx),(0:-dx:-lx)]';
  
    for i=2:2:2*n
        tc{i}=(x+1i*y)+(i-2)/2*lx;
    end   
end