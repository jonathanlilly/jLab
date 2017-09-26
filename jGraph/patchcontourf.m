function[h]=patchcontourf(varargin)
%PATCHCONTOURF  Generate filled contours using patches, with specified colors.
%
%   Frequently one wishes to make a filled contour plot of one field over
%   another with two different colormap schemes.  This is not currently 
%   possible with CONTOURF because Matlab uses one colormap per axes.
%
%   For example, one may wish to plot sea surface temperature filled in 
%   color with the continents filled in black.
%
%   PATCHCONTOURS gets around this problem by geneating filled contours 
%   plots as patch objects, with out reference to the current colormap.  
%   Colors of filled contours may be specified explicitly. 
%
%   PATCHCONTOURF(Z,V) makes a filled contour plot of field Z at value V. 
%   The interiors of the contour are colored black by default.
%
%   PATCHCONTOURF(Z,V,C) alternately uses color C as both the fill color 
%   and the edge color.  C may be either a string, as in 'g' for green, 
%   or a color array of length 3, such as [1/2 1/2 1/2] for gray.
%
%   PATCHCONTOURF(X,Y,Z,V) and PATCHCONTOURF(X,Y,Z,V,C) use vectors X and Y
%   to specify the contour locations.  X has length SIZE(X,2) while Y has
%   length SIZE(X,1).  
%
%   H=PATCHCONTOURF(...) returns a handle to the patch object.
% 
%   PATCHCONTOURF(...,'m_map') alternately works with R. Pawlowicz's M_MAP
%   mapping package, available at http://www.eos.ubc.ca/~rich/map.html.
%
%   'patchcontourf --f' generates a sample figure.
%
%   Usage: patchcontourf(z,v);
%          patchcontourf(z,v,c);
%          patchcontourf(x,y,z,v,c);
%          h=patchcontourf(x,y,z,v,c);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--f')
    type makefigs_patchcontourf
    makefigs_patchcontourf;
    return
end


%Sort out string input arguments
str='k';
flag='mat';
if isstr(varargin{end})
    if length(varargin{end})>=3
        if strcmpi(varargin{end}(1:3),'m_m')||strcmpi(varargin{end}(1:3),'mat')
            flag=varargin{end};
            varargin=varargin(1:end-1);
        end
    end
end

if isstr(varargin{end})||length(varargin{end})==3
    str=varargin{end};
    varargin=varargin(1:end-1);
end

v=varargin{end};
varargin=varargin(1:end-1);

if length(varargin)==3
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
else
    x=[];
    y=[];
    z=varargin{1};
end

z1=zeros(size(z,1)+2,size(z,2)+2)+minmin(z);
z1(2:end-1,2:end-1)=z;

[ii,jj]=closedcurves(z1,v);

bool=ishold;
hold on

if isempty(x)
    for i=1:length(ii)
        if strcmpi(flag(1:3),'m_m')
            m_patch(ii{i},jj{i},str,'edgecolor',str);
        else
            patch(ii{i},jj{i},str,'edgecolor',str);
        end
    end
else
    dx1=x(2)-x(1);
    dxN=x(end-1)-x(end);
    x=[x(1)-dx1;x(:);x(end)-dxN];
    
    dy1=y(2)-y(1);
    dyN=y(end-1)-y(end);
    y=[y(1)-dy1;y;y(end)-dyN];
    
    for i=1:length(ii)
        x1=interp1(1:length(x),x,ii{i},'linear');
        y1=interp1(1:length(y),y,jj{i},'linear');
        if strcmpi(flag(1:3),'m_m')
            %The FLIPUDs keep M_PATCH from inverting the patch
            h=m_patch(flipud(x1),flipud(y1),str,'edgecolor',str);
        else
            h=patch(x1,y1,str,'edgecolor',str);
        end
    end
end

if ~ishold 
    hold off
end


