%JLAB_ADDPATH  Adds JLAB and JDATA subdirectories to your Matlab search path. 
%
%   This script adds JLAB and JDATA subdirectories to your Matlab search path.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details

%First, find the path to JLAB
fullpath=which('jlab_license');
jlab_path=fullpath(1:end-15);

for k=1:2
    if k==2
        jlab_path=[jlab_path(1:end-4) 'jdata'];
        addpath(jlab_path)
    end
    dirlist=dir(jlab_path);
    for i=1:length(dirlist)
        if dirlist(i).isdir
            if ~strcmpi(dirlist(i).name(1),'.')
                %dirlist(i).name
                addpath([jlab_path '/' dirlist(i).name])
                dirlisti=dir([jlab_path '/' dirlist(i).name]);
                for j=1:length(dirlisti)
                    if dirlisti(j).isdir
                        if ~strcmpi(dirlisti(j).name(1),'.')
                            addpath([jlab_path '/' dirlist(i).name '/' dirlisti(j).name]);
                        end
                    end
                end
            end
        end
    end
end

clear i dirlist jlab_path fullpath