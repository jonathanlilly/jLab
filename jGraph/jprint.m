function[varargout]=jprint(varargin)
%JPRINT  Print to a specified directory and crop the resulting file.
%
%   JPRINT(DIRNAME,FILENAME) prints the current figure window to directory
%   DIRNAME in PNG format to file FILENAME, and crops the result using the
%   function CROP by Andy Bliss. 
%
%   JPRINT(DIRNAME,FILENAME,STR) instead prints with the file format STR,
%   with STR='png' being the default behavior.
%
%   Usage: jprint(dirname,filename)
%          jprint(dirname,filename,'jpeg')
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2019 J.M. Lilly --- type 'help jlab_license' for details

dir=varargin{1};
filename=varargin{2};
str='png';
if nargin>2
    str=varargin{3};
end

olddir=pwd;
cd(dir)
eval(['print -d' str ' ' filename])
fullfilename=findfiles(pwd,'*','include',filename);
bool=false(length(fullfilename),1);
for i=1:length(fullfilename)
    bool(i)=aresame(fullfilename{i}(1:length(filename)),filename);
end
fullfilename=fullfilename(bool);
eval(['crop ' fullfilename{1}])
cd(olddir)