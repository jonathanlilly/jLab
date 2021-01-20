function[struct]=structindex(struct,index)
%STRUCTINDEX  Applies an index to all array-valued fields in a structure. 
%
%   STRUCTOUT=STRUCTINDEX(STRUCT,INDEX), where STRUCT is a structure
%   containing a set of fields that are non-string arrays all having the 
%   same length N>1, applies INDEX to all of them and returns the result in
%   STRUCTOUT.  All other fields of STRUCT are left unchanged. 
%
%   STRUCTINDEX(STRUCT,BOOL) where BOOL is a boolean array of length N also 
%   works.  
%
%   STRUCTINDEX is useful for dealing with the ridge structures output by 
%   EDDYRIDGES.
%
%   See also EDDYRIDGES.
%
%   Usage: struct=structindex(struct,index);
%          struct=structindex(struct,bool);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2019 J.M. Lilly --- type 'help jlab_license' for details
 
names=fieldnames(struct);
for i=1:length(names)
    temp=getfield(struct,names{i});
    if length(temp)~=1&&~ischar(temp)
        %names{i}
        temp=temp(index);
        struct=setfield(struct,names{i},temp);
    end
end