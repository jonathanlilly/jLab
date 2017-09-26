function[str]=commentlines(dirname)
%COMMENTLINES  Returns the comment lines from m-files. 
%
%   COMMENTLINES is used to facilitate making 'Contents.m' files.
%
%   COMMENTLINES returns a string matrix which contains all the comment
%   lines, i.e the first line beginning with '%', from all m-files in the
%   current directory.  
%  
%   COMMENTLINES DIRNAME applies to directory DIRNAME. 
%   
%   COMMENTLINES with no input argument appies itself to the current
%   working directory.  
%   
%   COMMENTLINES FILENAME applies to just m-file FILENAME.
%  
%   Examples:  commentlines plot
%              commentlines elmat
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002--2015 J.M. Lilly --- type 'help jlab_license' for details

bmname=0;
dirname=deblank(dirname);

if nargin~=0
  if exist(dirname)==2 || exist(dirname)==5
     x=dirname;
     if ~ismname(x)
         x=[x '.m'];
     end    
     bmname=1;
  elseif exist(dirname)==7
     olddir=pwd;
     evalin('base',['cd ' findpath(dirname)])
  end
else
  dirname=pwd;
end

str=[];
%/********************************************************
%Make a list of all m-file names
if ~bmname
  x=dir;
  for i=1:length(x)
    y{i}=x(i).name;
  end
  x=y;
  
  %remove elements which are not .m 
  x=x(ismname(x));
  x=strs2mat(x);
end
%\********************************************************


if isempty(x)
  disp(['No m-files found in directory ' dirname '.'])
else
  for i=1:size(x,1)
    fid=fopen(deblank(x(i,:)));
    firstline=fgetl(fid);
    secondline=fgetl(fid);
    if length(firstline)>=8
       bfunction=strcmpi(firstline(1:8),'function');
    else
       bfunction=0;
    end
    
    if ~isempty(secondline)
       if aresame(secondline,-1)|| ~bfunction
           secondline=firstline;       
       end
    else
        secondline=firstline;
    end
    
  
    if ~isempty(secondline)
      %remove both leading and trailing blanks
      a=find(real(secondline)~=real('%') & ~isblank(secondline),1,'first');
      temp=deblank(secondline(a:end));
      str{i}=temp;
    else 
      str{i}=[];
    end
    
    fclose(fid);
  end
  
  str=packstrs(str);
    
  for i=1:length(str)
    
      str1=str{i};
      a=find(real(str1)~=real('%'),1,'first');
      b=find(isblank(str1)|istab(str1),1,'first')-1;
      
      name=lower(safeindex(a,b,str1));
      str1=safeindex(b+1,length(str1),str1);
       
      if ~isempty(str1)
         a=find(~isblank(str1) & ~istab(str1) & real(str1)~=real('%'),1,'first');
         str1=safeindex(a,length(str1),str1);
      end
      
      if ~isempty(name)
        nfl=length(name);
        if nfl<13
            str1=['%   ' name blanks(13-nfl) ' - ' str1];
        else
            str1=['%   ' name ' - ' str1];
        end
        str{i}=str1;
      else 
        str{i}=[];
      end  
  end
  
  str=packstrs(str);
  str=strs2mat(str);
  
  if ~isempty(str)
      firstcol=real(str(:,1));
      [firstcol,index]=sort(firstcol);
      str=str(index,:);
  end
end %m-files found

%Upcase initial word
% for i=1:size(str,1)
%    index=strfind(str(i,:),' ');
%    index=index(find(index>4));
%    if ~isempty(index);
%       str(i,1:index-1)=upper(str(i,1:index-1));
%    end
% end

if nargin~=0 && ~bmname
   evalin('base',['cd ' olddir])
end


function[y]=safeindex(a,b,x)
% SAFEINDEX   Indexes an array; returns empty if index is empty
%
%   Y=SAFEINDEX(INDEX,X) <==>  Y=X(INDEX), or [] if empty INDEX
%   Y=SAFEINDEX(A,B,X) <==>  Y=X(A:B) if B>=A, [] otherwise
 
y=[];
 
if nargin==2
  index=a;
  x=b;
  if ~isempty(a)
    y=x(index);
  end
elseif b-a>=0
   index=a:b;
   y=x(index);
end

function[mat]=strs2mat(x,flag)
%STRS2MAT  Converts a cell array of strings into a string matrix.
%
%   MAT=STRS2MAT(X) where X is a cell array of strings will return a
%   matrix MAT, with rows containing the elements of the cell array.
%
%   The text in MAT is formatted to be flush on the left.
%
%   See also FLUSHLEFT, FLUSHRIGHT
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2003--2006 J.M. Lilly --- type 'help jlab_license' for details      
  
%Flag prevents recursion from flushleft
if nargin==1
  flag=0;
end

if ~iscell(x)
  xtemp=x;
  clear x
  x{1}=xtemp;
end

%x=packstrs(x);

if isempty(x)
    mat=[];
else
    M=length(x);
    for i=1:M
        n(i)=length(x{i});
    end
    N=max(n);
    
    mat=char(32*ones(M,N)); %32 <==> ' '
    
    for i=1:M
        if n(i)>0
            mat(i,1:n(i))=x{i};
        end
    end
    
    if ~flag
        mat=flushleft(mat);
    end
end


function[x]=packstrs(x)
%PACKSTRS Removes empty entries from a cell array of strings
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2008 J.M. Lilly --- type 'help jlab_license' for details        
  
x=deblankstrs(x);
n=cellength(x);
x=x(n>=1);

function[x]=deblankstrs(x)
%DEBLANKSTRS Deblanks all elements in a cell array of strings
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details  
  
for i=1:length(x)
  x{i}=deblank(x{i});
end

function[b]=isblank(x)
%ISBLANK  Tests whether elements of a string are blanks.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details  
 
x=real(x);
blk=real(' ');
b=(x==blk);

function[b]=ismname(x)
%ISMNAME  Tests whether a string is an m-file name; cells OK  
%
%   ISMNAME STR tests whether a string STR has the ending ".m".
%
%   ISMNAME(X) where X is a cell array of strings returns a boolean
%   array testing whether each string has the ending ".m".
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details  
  
if iscell(x)
  for i=1:length(x)
   b(i)=ismname1(x{i});
  end
else
   b=ismname1(x);
end

function[b]=ismname1(x)
 
x=deblank(x);
if length(x)>=2
    b=strcmpi(x(end-1:end),'.m');
else 
    b=false;
end

if strcmpi(x(1),'#') || strcmpi(x(1),'.')  %Exclude funny names
   b=false;
end

function[b]=istab(x)
% ISTAB  Tests whether elements of a string are tab markers.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details        
 
  
x=real(x);
blk=9;
b=(x==blk);

function[x]=flushleft(x)
%FLUSHLEFT   Makes a blank-padded string matrix flush on the left
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002--2015 J.M. Lilly --- type 'help jlab_license' for details  

for i=1:size(x,1)
  index=find(real(x(i,:))~=real(' '),1,'first');
  if ~isempty(index)
    x(i,:)=vshift(x(i,:),index-1,2);
  end
end


