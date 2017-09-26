function[dir]=whichdir(name)
%   WHICHDIR  Returns directory name containing file in search path.
%
%   WHICHDIR NAME returns the directory containing a file called NAME,
%   i.e. the full path of the file excluding NAME and the trailing '/'. 
%
%   If more than one directory contains a file called NAME, then 
%   WHICHDIR returns a cell array of directories names.
%
%   If NAME does not include an extension, then it is interpreted as
%   the name of an m-file, i.e. having extension '.m' .
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2006--2015 J.M. Lilly --- type 'help jlab_license' for details
  
if isempty(strfind(name,'.'))
    name=[name '.m'];
end

dir=which(name,'-all');

for i=1:length(dir)
   dir{i}=dir{i}(1:end-length(name)-1);
end

if i==1
    dir=dir{1};
end
