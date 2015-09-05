function[fullpath]=findpath(dirname)
%FINDPATH  Returns the full pathname of a directory on the Matlab search path.
%
%   FINDPATH(DIRNAME) returns the full pathname of directory DIRNAME,
%   provided DIRNAME is on the Matlab search path.
%
%   For example, 'findpath jlab' on a Mac might return
%
%            '/Users/lilly/Desktop/DropBox/Matlab/jlab'.
%
%   FINDPATH returns the empty string if DIRNAME is not found.
%
%   Usage: findpath('dirname');
%          findpath dirname
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
  
pathstr=path;
b=strfind(pathstr,':');
a=[1 b(1:end-1)+1];

fullpath=[];
for i=1:length(a)
    str=pathstr(a(i):b(i));
    index=strfind(str,[dirname ':']);
    if ~isempty(index)
        fullpath=str(1:end-1);
    end
end

