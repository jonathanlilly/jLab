function[]=noxlabels(str)
%NOXLABELS   Remove some or all x-axis tick mark labels.
%
%   NOXLABELS by itself removes all tick mark labels on the current x-axis.
%
%   NOXLABELS MAX removes only the label of the maximum x-tick mark.
%
%   NOXLABELS MIN removes only the label of the minimum x-tick mark.
%
%   NOXLABELS ALL is the same as the default.
%
%   See also NOYLABELS.
% 
%   Usage:  noxlabels
%           noxlabels xmin
%           noxlabels ymin
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2008 J.M. Lilly --- type 'help jlab_license' for details        

if nargin==0
    str='all';
end

if strcmpi(str(1:3),'min')
    xlabels=get(gca,'xticklabel');
    xlabels(1,:)=' ';
    set(gca,'xticklabel',xlabels)
elseif strcmpi(str(1:3),'max')
    xlabels=get(gca,'xticklabel');
    xlabels(end,:)=' ';
    set(gca,'xticklabel',xlabels)
elseif strcmpi(str(1:3),'all')
    set(gca,'xticklabel',[])
end
