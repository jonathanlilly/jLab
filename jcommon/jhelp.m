function[varargout]=jhelp(str)
%JHELP  Opens linked JLAB help files in Matlab's internal web browser.
%
%   JHELP with no input arguments opens the JLAB Contents.m file.
%   JHELP FILE opens the help information for JLAB function FILE.
%   
%   Usage: jhelp
%          jhelp polysmooth
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2015 J.M. Lilly --- type 'help jlab_license' for details

if nargin==0
    str='jLab';
end

jlabdir=whichdir('jlab_license.m');
web([jlabdir '/doc/' str '.html'])

