function[varargout]=findfiles(varargin)
%FINDFILES  Returns all files in a directory with a specified extension.
%
%   FILES=FINDFILES(DIRNAME,EXT) where DIRNAME is a full directory path, 
%   returns a cell arrays FILES of all files having extension EXT in 
%   directory DIRNAME.
%
%   For example, FILES=FINDFILES('/Users/lilly/Home','m') returns a cell
%   cell array of all m-files in the Home directory.
%
%   FINDFILES(DIRNAME,'*') returns files having any extension.  
%
%   FINDFILES(DIRNAME,'') returns files having no extension.  This can
%   be used to return directory names.
%   ___________________________________________________________________
%
%   Recursive searching
%
%   FINDFILES(DIRNAME,...,'recursive') recursively searches though
%   directory DIRNAME as well as any subdirectories it contains. 
%   
%   [FILES,PATHS]=FINDFILES(DIRNAME,...,'recursive') also outputs a cell 
%   array PATHS specifying the full path name to each file in FILES. 
%
%   FINDFILES(DIRNAME,...,'recursive','ignore_private') will ignore any 
%   directories named 'private' as a part of this recursive search. 
%   ___________________________________________________________________
%
%   Additional options
%
%   FINDFILES(...,'include',STR1) additionally returns only those file
%   names which inlude the text STR1.
%
%   FINDFILES(...,'exclude',STR2) excludes those files names which 
%   include the text STR2.
%  
%   The include and exclude options must be input after the first two
%   variables.  Both options may be input simultaneously.
%   ___________________________________________________________________
%
%   Usage: files=findfiles(pwd,'m');
%          files=findfiles(dirname,ext);    
%          files=findfiles(dirname,ext,'include','help');
%          [files,dirs]=findfiles(dirname,ext,'include','help','recursive');
%          files=findfiles(dirname,ext,'include','help','exclude','jlab');
%
%   'findfiles --t' runs a test.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2019 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    findfiles_test,return
end
 
dirname=varargin{1};
ext=varargin{2};

if length(ext)>1
    if strcmpi(ext(1),'.'),ext=ext(2:end);end
end
if strcmpi(dirname(end),'/'),dirname=dirname(1:end-1);end

vars=varargin(3:end);
excludestr=[];
includestr=[];
ignorestr='searchprivate';
recursionstr='norecursion';
while ~isempty(vars)
    str=vars{1};
    if strcmpi(str(1:3),'exc')
        excludestr=vars{2};
        vars=vars(3:end);
    elseif strcmpi(str(1:3),'inc')
        includestr=vars{2};
        vars=vars(3:end);
    elseif strcmpi(str(1:3),'rec')
        recursionstr=vars{1};
        vars=vars(2:end);
    elseif strcmpi(str(1:3),'ign')
        ignorestr=vars{1};
        vars=vars(2:end);
    end
end

dirstruct=dir(dirname);
files=cell(length(dirstruct),1);
isdir=false(length(dirstruct),1);
for i=1:length(dirstruct)
    files{i}=dirstruct(i).name;
    isdir(i)=dirstruct(i).isdir;
end
dirs=files(isdir);


bool=true(size(dirs));
for i=1:length(dirs)
    if strcmpi(dirs{i},'.')||strcmpi(dirs{i},'..')
        bool(i)=false;
    end
    if ignorestr
        if strcmpi(dirs{i},'private')
            bool(i)=false;
        end
    end
end
dirs=dirs(bool);

bool=false(length(files),1);
N=length(ext);
if N==0
    for i=1:length(files)
        bool(i,1)=isempty(strfind(files{i},'.'));
    end
elseif N==1 && strcmpi(ext,'*')
    for i=1:length(files)
        bool(i,1)=~isempty(strfind(files{i},'.'))&~strcmpi(files{i}(end),'.');
    end
else   
    for i=1:length(files)
        if length(files{i})>N
            bool(i,1)=strcmpi(files{i}(end-N:end),['.' ext]);
        end
    end
end
if N~=0
    for i=1:length(files)
        if ~isempty(includestr)
            bool(i,1)=bool(i,1)&contains(files{i}(1:end-N-1),includestr);
        end
        if ~isempty(excludestr)
            bool(i,1)=bool(i,1)&~contains(files{i}(1:end-N-1),excludestr);
        end
    end
end

index=find(bool);
if ~isempty(index)
    files=files(index);
else 
    files=[];
end

dirs_out=files;
for i=1:length(files)
    dirs_out{i}=dirname;
end

if strcmpi(recursionstr(1:3),'rec')
    for i=1:length(dirs)
        %[dirname '/' dirs{i}]
        [filesi,dirsi]=findfiles([dirname '/' dirs{i}],ext,'include',includestr,'exclude',excludestr,recursionstr);
        files=[files;filesi];
        dirs_out=[dirs_out;dirsi];
    end
end
                
varargout{1}=files;
varargout{2}=dirs_out;

function[]=findfiles_test

dirname=whichdir('jlab_license');
if iscell(dirname)
    dirname=dirname{1};
end
files=findfiles(dirname,'m');
bool=0;
for i=1:length(files)
    if strcmpi('jlab_license.m',files{i})
        bool=1;
    end
end
reporttest('FINDFILES found jlab_license', bool)
%files=findfiles(dirname,'m','include','jlab_license');
%reporttest('FINDFILES with include flag', aresame(files{1},'jlab_license.m'))

