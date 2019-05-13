function[]=jlab_addpath
%JLAB_ADDPATH  Adds JLAB and JDATA subdirectories to your Matlab search path. 
%
%   This script adds JLAB subdirectories to your Matlab search path, as
%   well as the JDATA directory if it is found to exist. 
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015--2019 J.M. Lilly --- type 'help jlab_license' for details

%First, find the path to JLAB
fullpath=which('jlab_license');
jlab_path=fullpath(1:end-15);
bool=real(jlab_path)==real('\');
jlab_path(bool)='/';

%jlab_path

for k=1:2
    if k==1
        disp('Adding JLAB subdirectories to your Matlab search path.')
        if exist(jlab_path)
            jlab_addpath_subdirectories(jlab_path);
        end
    end
    if k==2
        jlab_path=[jlab_path(1:end-5) '/jdata'];
        if exist(jlab_path)
            disp('Looks like you have JDATA installed; adding this to your path also.')
            addpath(jlab_path)
            jlab_addpath_subdirectories(jlab_path);
        end
    end
end
clear i j k dirlist dirlisti jlab_path fullpath bool

function[]=jlab_addpath_subdirectories(path)
dirlist=dir(path);
for i=1:length(dirlist)
    if dirlist(i).isdir
        if ~strcmpi(dirlist(i).name(1),'.')
            %dirlist(i).name
            addpath([path '/' dirlist(i).name])
            dirlisti=dir([path '/' dirlist(i).name]);
            for j=1:length(dirlisti)
                if dirlisti(j).isdir
                    if ~strcmpi(dirlisti(j).name(1),'.')
                        %[path '/' dirlist(i).name '/' dirlisti(j).name]
                        addpath([path '/' dirlist(i).name '/' dirlisti(j).name]);
                    end
                end
            end
        end
    end
end