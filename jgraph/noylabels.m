function[]=noylabels(str)
%NOYLABELS   Remove some or all y-axis tick mark labels.
%
%   NOYLABELS by itself removes all tick mark labels on the current y-axis.
%
%   NOYLABELS MAX removes only the label of the maximum y-tick mark.
%
%   NOYLABELS MIN removes only the label of the minimum y-tick mark.
%
%   NOYLABELS ALL is the same as the default.
%
%   See also NOXLABELS.
% 
%   Usage:  noylabels
%           noylabels xmin
%           noylabels ymin
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2008 J.M. Lilly --- type 'help jlab_license' for details            

if nargin==0
    str='all';
end

if strcmpi(str(1:3),'min')
    ylabels=get(gca,'yticklabel');
    ylabels(1,:)=' ';
    set(gca,'yticklabel',ylabels)
elseif strcmpi(str(1:3),'max')
    ylabels=get(gca,'yticklabel');
    ylabels(end,:)=' ';
    set(gca,'yticklabel',ylabels)
elseif strcmpi(str(1:3),'all')
    set(gca,'yticklabel',[])
end

