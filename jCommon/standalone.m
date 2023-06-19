function[mfiles]=standalone(varargin)
%STANDALONE  Create stand-alone version of an m-file, including dependencies.
%
%   DEPENDS=STANDALONE(FILENAME,SOURCEDIR) looks in m-file FILENAME for all
%   calls to m-files in directory SOURCEDIR and subdirectories, and returns 
%   in cell array DEPENDS a list of m-files called by FILENAME, or called
%   by any m-files called by those files, etc., excluding FILENAME itself.
%
%   Thus, in order to create a stand-alone version of FILENAME, one needs
%   to include all the files in DEPENDS from SOURCEDIR and subdirectories.
%   These files are called the 'dependencies' of FILENAME.
%   
%   Note that SOURCEDIR must be the full pathname of a directory, e.g.
%   SOURCEDIR='/Users/lilly/Desktop/Dropbox/Matlab/jlab'.
%
%   FILENAME may either be a string, or a cell array of strings containing
%   multiple filenames.  The trailing '.m' extension is not needed.   
%   ___________________________________________________________________
%
%   Copying files
%
%   STANDALONE(FILENAME,SOURCEDIR,TARGETDIR) on Mac or Unix systems copies
%   all the files in FILENAME to directory TARGETDIR, which is again the 
%   full pathname of a directory, and copies all files in DEPENDS to a
%   subdirectory within TARGETDIR named 'private'.  
%
%   The directories are created if they do not already exist. 
%
%   The subdirectory 'private' has a special meaning in Matlab.  Functions
%   in this folder can only be called by functions in the enclosing 
%   directory, or by scripts called by these functions. Thus, the 
%   dependency files cannot be accidentally called by other functions.  
%   
%   WARNING: This routine will automatically copy potentially large numbers
%   of files by directly accessing the 'cp' command.  If you make an error 
%   in the input directories, you could end up with a big mess, or even
%   overwrite things you don't want to overwrite.  Please use with caution.
%
%   For example, to make a standalone version of this function, type
%
%        standalone('standalone',....,
%               '/Users/lilly/Desktop/DropBox/Matlab/jlab',...
%               '/Users/lilly/Desktop/DropBox/Matlab/standalone')
%
%   which copies 'standalone.m' into a directory 'standalone', and its
%   dependencies into the subdirectory 'standalone/private'. 
%   ___________________________________________________________________
%
%   Usage: depends=standalone(filename,sourcedir);
%          standalone(filename,sourcedir,targetdir)
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015--2017 J.M. Lilly --- type 'help jlab_license' for details
 
brecursion=true;
if ~ischar(varargin{end})
    brecursion=varargin{end};
    varargin=varargin(1:end-1);
end

filename=varargin{1};
sourcedir=varargin{2};
    
targetdir=[];
if length(varargin)>2
    targetdir=varargin{3};
end

if ~isempty(targetdir)
    if ~ismac||~isunix
        error('Sorry, STANDALONE can only move files on Mac or Unix systems.')
    end
end

%Make sure the input filename has the right format
if iscell(filename)
    filename=filename(:);
end
if ~iscell(filename)
    filenames{1}=filename;
else
    filenames=filename;
end
%Remove trailing m-file extension
for i=1:length(filenames)
    if strcmpi(filenames{i}(end-1:end),'.m')
        filenames{i}=filenames{i}(1:end-2);
    end
end
%Remove any accidentally duplicated filenames
for i=1:length(filenames)
    for j=i+1:length(filenames)
        if strcmpi(filenames{i},filenames{j})
            filenames{j}=[];
        end
    end
end
filenames=cellprune(filenames,'quiet');

%Read contents of the source file(s)
text=[];
for i=1:length(filenames)
    fullsourcedir=whichdir(filenames{i});
    if iscell(fullsourcedir)
        fullsourcedir=fullsourcedir{1};
    end
    fullfilename=[fullsourcedir '/' filenames{i} '.m'];
    fid=fopen(fullfilename,'r');
    text=[text char(fread(fid)')];
    fclose(fid);
end

%Find all m-files in sourcedirectory
mfiles=findfiles(sourcedir,'m','recursive');
%Remove the m-extension
for i=1:length(mfiles)
    mfiles{i}=mfiles{i}(1:end-2);
end

%Prune list of file names
bool=true(size(mfiles));
for i=1:length(mfiles)
    %Remove names of this or other input filenames
    for j=1:length(filenames)
        if strcmpi(mfiles{i},filenames{j})
            bool(i)=false;
        end
    end
end
mfiles=mfiles(bool);

%Find all calls
bool=false(size(mfiles));
for i=1:length(mfiles)
    if ~isempty(strfind(text,mfiles{i}))
       %True if any usage is a non-commented call
       bool(i)=~iscommented(text,mfiles{i});
    end
end
mfiles=mfiles(bool);
%return   %XXX  

%Remove duplicates from multiple calls
mfiles=remove_repeated(mfiles);
mfiles=remove_duplicates(filenames,mfiles);

%These are now the m-files used by filenames, that were not previously 
%included in a higher-level recursion, and not including files in filenames

%And recurse, one depth level at a time
if brecursion
    %depth=0;
    bdone=false;
    mfiles_previous=mfiles;  %Files from previous depth level
    while ~bdone
        %depth=depth+1
        %mfiles_previous,sourcedir
        mfiles_next=standalone(mfiles_previous,sourcedir,false);
        %Files from next depth level
        
        %Remove duplicates
        mfiles_next=remove_duplicates(mfiles,mfiles_next);
        mfiles_next=remove_duplicates(filenames,mfiles_next);
        
        if isempty(mfiles_next)
            bdone=true;
        else
            %Accumulate into a list of all files called at any level
            mfiles=vertcat(mfiles,mfiles_next);
            mfiles_previous=mfiles_next;
        end
    end
end
%Variable mfiles is now a list of called m-files at any recursion depth

if ~isempty(targetdir)
    movefiles(filenames,targetdir)
    movefiles(mfiles,[targetdir '/private'])
end

if nargout==0
    clear mfiles
end

function[x]=remove_repeated(x)
%Removes repeated entries in x
bool=true(size(x));
for i=1:length(x)
    for j=i+1:length(x)
        if strcmpi(x{i},x{j})
            bool(j)=false;
        end
    end
end
x=x(bool);

function[y]=remove_duplicates(x,y)
%Removes all entries in y that already exist in x
bool=true(size(y));
for i=1:length(x)
    for j=1:length(y)
        if strcmpi(x{i},y{j})
            bool(j)=false;
        end
    end
end
y=y(bool);

function[]=movefiles(mfiles,targetdir)
%Moves a list of files into targetdir
currentdir=pwd;
try
    cd(targetdir)
    cd currentdir
catch
    disp(['Directory ' targetdir ' does not exist; creating.'])
    eval(['!mkdir ' targetdir])
end

disp(['STANDALONE copying ' int2str(length(mfiles)) ' files.'])
for i=1:length(mfiles)
    fullsourcedir=whichdir([mfiles{i} '.m']);
    if iscell(fullsourcedir)
        fullsourcedir=fullsourcedir{1};
    end
    fullfilename=[fullsourcedir '/' mfiles{i} '.m'];
    eval(['!cp ' fullfilename ' ' targetdir])
end

function[bool]=iscommented(text,str)
%Returns true if all occurrences of STR within TEXT are commented out
index=strfind(text,str);
bool=false(size(index));
for i=1:length(index)
     textuntil=text(1:index(i));
     lastcomment=find(real(textuntil)==real('%'),1,'last');
     lastendline=find(real(textuntil)==10,1,'last');
     if isempty(lastendline)&&~isempty(lastcomment)
         bool(i)=true;
     elseif ~isempty(lastendline)&&~isempty(lastcomment)
         bool(i)=(lastendline<lastcomment);
     end
end
bool=prod(bool);


