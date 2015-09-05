function[mfiles]=standalone(varargin)
%STANDALONE  Create stand-alone version of an m-file, including dependencies.
%
%   FILES=STANDALONE(FILENAME,SOURCEDIR) looks in m-file FILENAME for all
%   calls to m-files in directory SOURCEDIR and subdirectories, and returns 
%   a list of m-files called by FILENAME, or called by any m-files called 
%   by those files, and so forth. Commented-out calls are not included. 
%
%   This list of called m-files is returned in the cell array FILES.  The
%   first element of FILES is the original m-file FILENAME itself. 
%
%   Thus, in order to create a stand-alone version of FILENAME, one needs
%   to include all the files in FILES from SOURCEDIR and subdirectories.
%   
%   Note that SOURCEDIR must be the full pathname of a directory, e.g.
%   SOURCEDIR='/Users/lilly/Desktop/Dropbox/Matlab/jlab'.
%   
%   STANDALONE(FILENAME,SOURCEDIR,TARGETDIR) on Mac or Unix systems moves
%   all the files to directory TARGETDIR, which is again the full pathname 
%   of a directory.  TARGETDIR is created if does not already exist. 
%   
%   WARNING: This routine will automatically copy potentially large numbers
%   of files by directly accessing the 'cp' command.  If you make an error 
%   in the input directories, you could end up with a big mess, or even
%   overwrite things you don't want to overwrite.  Please use with caution.
%
%   For example, to make standalone version of STANDALONE, one would type
%
%        standalone('standalone.m',....,
%               '/Users/lilly/Desktop/DropBox/Matlab/jlab',...
%               '/Users/lilly/Desktop/DropBox/Matlab/standalone')
%
%   which moves STANDALONE and all its dependencies into a directory of the
%   same name. 
%
%   Usage: files=standalone(filename,sourcedir);
%          standalone(filename,sourcedir,targetdir)
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
 

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


%Read contents of the source file
fullsourcedir=whichdir(filename);
fullfilename=[fullsourcedir '/' filename];
fid=fopen(fullfilename,'r');
text=char(fread(fid)');
fclose(fid);

%Find all m-files in sourcedirectory
mfiles=findfiles(sourcedir,'m','recursive');

%Remove self name
bool=true(size(mfiles));
for i=1:length(mfiles)
    if strcmpi(mfiles{i},filename)
        bool(i)=false;
    end
end
mfiles=mfiles(bool);

        
%Remove the m-extension
for i=1:length(mfiles)
    mfiles{i}=mfiles{i}(1:end-2);
end

%Find all calls
bool=false(size(mfiles));
for i=1:length(mfiles)
    if ~isempty(strfind(text,mfiles{i}))
       %True if any usage is a non-commented call
       bool(i)=~iscommented(text,mfiles{i});
    end
end
mfiles=mfiles(bool);  %These are now the m-files used by filename

%And recurse
for i=1:length(mfiles)
      mfilesi=standalone([mfiles{i} '.m'],sourcedir);
      mfiles=[mfiles;mfilesi];
end

if ~isempty(targetdir)
    mfiles=[filename(1:end-2);mfiles];
    
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
        fullfilename=[fullsourcedir '/' mfiles{i} '.m'];
        eval(['!cp ' fullfilename ' ' targetdir])
    end
end

if nargout==0
    clear mfiles
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


