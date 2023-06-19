function[varargout]=jprint(varargin)
%JPRINT  Print to a specified directory and crop the resulting file.
%
%   JPRINT(DIRNAME,FILENAME) prints the current figure window to directory
%   DIRNAME in PNG format to file FILENAME, and crops the result using the
%   function CROP by Andy Bliss. 
%
%   JPRINT(DIRNAME,FILENAME,FORMAT) prints to the file format FORMAT, with
%   with FORMAT='png' being the default behavior.
%
%   JPRINT(...,STR) where STR begins with a '-' as in STR='-r200', passes
%   this argument to the print function.
%
%   Usage: jprint(dirname,filename)
%          jprint(dirname,filename,'jpeg')
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2019--2023 J.M. Lilly --- type 'help jlab_license' for details

dir=varargin{1};
filename=varargin{2};
formatstr='png';
str='';
varargin=varargin(3:end);
for i=1:length(varargin)
    if strcmpi(varargin{i}(1),'-')
        str=varargin{i};
    else
        formatstr=varargin{i};
    end
end

olddir=pwd;
cd(dir)

eval(append('print ',str,' -d',formatstr,' ',filename))
fullfilename=findfiles(pwd,formatstr,'include',filename);

if iscell(fullfilename)
    fullfilename=fullfilename{1};
end

%fullfilename=[pwd '/' filename '.' formatstr]
%bool=false(length(fullfilename),1);
% for i=1:length(fullfilename)
%     fullfilename{i}(1:length(filename))
%     bool(i)=aresame(fullfilename{i}(1:length(filename)),filename);
% end
% fullfilename=fullfilename(bool);
if strcmpi(formatstr,'png')
    eval(['crop ' fullfilename])
end
cd(olddir)

