function[h]=seminfhaxby(N)
%SEMINFHAXBY  The seminf-Haxby colormap.
%   __________________________________________________________________
%
%   *|* seminfhaxby.png --- The seminf-haxby colormap versus lansey,
%   parula, and jet for a non-negative quantity.
%   Type 'jhelp seminfhaxby.png' to view this image. *|*
%   __________________________________________________________________
%
%   SEMINFHAXBY(M) returns an M-by-3 matrix containing the seminf-Haxby
%   colormap, see
%
%       http://soliton.vm.bytemark.co.uk/pub/cpt-city/jjg/misc/
%                                  tn/seminf-haxby.png.index.html
%
%   SEMINFHAXBY with no arguments returns a colormap having the same number 
%   of colors as the colormap of the current figure.
%
%   To make SEMINFHAXBY your default colormap, add to your startup.m file
%   the line "set(0,'DefaultFigureColormap',seminfhaxby)".
% 
%   'seminfhaxby --f' generates the figure shown above.
%
%   Usage: h=seminfhaxby(M);
%          colormap seminfhaxby
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2020 J.M. Lilly --- type 'help jlab_license' for details
 
if nargin>0
    if strcmpi(N, '--f')
       type makefigs_seminfhaxby
       makefigs_seminfhaxby;
       return
    end
end

if nargin < 1, N = size(get(gcf,'colormap'),1); end

h = cmap2linspecer(colorm(N));
h = cell2mat(h);

function[vOut] = cmap2linspecer(vIn) % changes the format from a double array to a cell array with the right format
vOut = cell(size(vIn,1),1);
for ii=1:size(vIn,1)
    vOut{ii} = vIn(ii,:);
end
 
function[cmap] = colorm(varargin)
n = varargin{1};
%frac=1;
cmapp= [0.95*[255,255,255];...%JML modifcation
208 216 251;...
186 197 247;...
143 161 241;...
97 122 236;...
0  39 224;...
25 101 240;...
12 129 248;...
24 175 255;...
49 190 255;...
67 202 255;...
96 225 240;...
105 235 225;...
123 235 200;...
138 236 174;...
172 245 168;...
205 255 162;...
223 245 141;...
240 236 120;...
247 215 103;...
255 189  86;...
255 160  68;...
244 116  74;...
238  79  77];

x = linspace(1,n,size(cmapp,1));
xi = 1:n;
cmap = zeros(n,3);
for ii=1:3
    cmap(:,ii) = pchip(x,cmapp(:,ii),xi);
end
cmap = cmap/255;


