function[]=jlab_runtests(str)
%JLAB_RUNTESTS  Runs a test suite for JLAB package.
%
%   JLAB_RUNTESTS runs automated tests for the JLAB package.
%
%   'jlab_runtests' runs all automated tests.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2015 J.M. Lilly --- type 'help jlab_license' for details      


% 'jlab_runtests figures' makes all sample figures.  Note that this is
% going to dump about 50 figures onto your desktop. 
 
if nargin==0
    str='tests';
end


disp(['This may take a few minutes.'])

a=now;
global BOOL_JLAB_RUNTEST
global FUNCTION_NUMBER
global FUNCTION_NUMBER_ARRAY

BOOL_JLAB_RUNTEST=[];
FUNCTION_NUMBER_ARRAY=[];

dirname=whichdir('jlab_license');
if iscell(dirname)
    dirname=dirname{1};
end
names=findfiles(dirname,'m');
for i=1:length(names)
    names{i}=names{i}(1:end-2);
    dirnames{i}=dirname;
end

%Recursively add m-files from directories in jlab directory
dirlist=dir(dirname);
n=length(names);
for i=1:length(dirlist)
    if dirlist(i).isdir
        if ~strcmpi(dirlist(i).name(1),'.')
            namesi=findfiles([dirname '/' dirlist(i).name],'m');
            for j=1:length(namesi)
                if ~strcmp(namesi{j}(1),'C')
                    n=n+1;
                    names{n}=namesi{j}(1:end-2);
                    dirnames{n}=[dirname '/' dirlist(i).name];
                end
            end
        end
    end
end

%Alphabetize
namecode=zeros(length(names),1);
for i=1:length(names)
    namecode(i)=real(lower(names{i}(1)));
end
[namecode,sorter]=sort(namecode);
names=names(sorter);
dirnames=dirnames(sorter);

testfailures=zeros(length(names),1);
for i=1:length(names)
    FUNCTION_NUMBER=i;
    if ~strcmpi(names{i},'jlab_runtests')&&~strcmpi(names{i},'Contents') %No recursion please
        fid=fopen([dirnames{i} '/' names{i} '.m'], 'r');   
        filestr=char(fread(fid)');
        fclose(fid);
        if strcmpi(str(1:3),'tes')||strcmpi(str(1:3),'bot')
            if ~isempty(strfind(filestr,'''--t'''))
                try
                    eval([names{i} '(''--t'');']);
                catch    
                    testfailures(i)=1;
                end
            end
        end
        if strcmpi(str(1:3),'fig')||strcmpi(str(1:3),'bot')
            if ~isempty(strfind(filestr,'''--f'''))            
                try
                    disp(['Generating figures for ' names{i} '...'])
                    eval([names{i} '(''--f'');']);
                catch    
                end
            end
        end
    end
end
%close all


if strcmpi(str(1:3),'tes')||strcmpi(str(1:3),'bot') 
    disp('---------------------------------------------')
    disp('Please email the following to eponym@jmlilly.net after your first run.')
    disp(['JLAB_RUNTESTS --- ' int2str(sum(BOOL_JLAB_RUNTEST)) ' of '  int2str(length(BOOL_JLAB_RUNTEST)) ' tests passed.'])
    if sum(~BOOL_JLAB_RUNTEST)>0
          disp('    Tests in the following routines failed:')
          for i=1:length(names)
              Nfailed=length(find((FUNCTION_NUMBER_ARRAY==i)&(BOOL_JLAB_RUNTEST==0)));
              if Nfailed>0
                  disp(['          ' names{i} ', ' int2str(Nfailed) ' test(s) failed.'])
              end
          end
    end
    if exist('jdata')~=7
        disp('JDATA directory not found.  That''s only a problem if you thought you installed it.')
        disp('Tests dependent upon the JDATA data will not execute.')
    end
    if vsum(testfailures,1)>0
        disp('    Tests in the following routines did not execute:')
        for i=1:length(names)
            if testfailures(i)==1
                disp(['          ' names{i} ' tests did not execute.' ])
            end
        end
    end
end
b=now;
disp(['JLAB_RUNTESTS took ' num2str((b-a)*24*60) ' minutes running on a ' computer ' running Matlab ' version '.'])
ver
