function[]=use(x)
%USE  Copies structure fields into named variables in workspace.
%
%   USE STRUCT, where STRUCT is a structure, copies all fields of the form
%   STRUCT.X into variables named X.
%  
%   This is useful for handling multiple datasets with the same variable
%   names.  The structures can be then kept in memory and 'mapped' into
%   variables as needed.
%
%   See also MAKE, MATSAVE.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2015 J.M. Lilly --- type 'help jlab_license' for details    


if strcmpi(x, '--t')
    use_test;return
end


%'if ~exist (' x ')==1, load 'x'  

clear str
str{1}    =['if ~isempty(' x '),'];
str{end+1}=['  ZZFNAMES=fieldnames(' x ');' ];
str{end+1}='  for ZZi=1:length(ZZFNAMES),';
str{end+1}=[' 	  eval([ZZFNAMES{ZZi}, ''=getfield(' x ',ZZFNAMES{ZZi});'']);'];
str{end+1}='  end;';
str{end+1}='else;';
str{end+1}='  disp([''Contains no data.'']);'; 
str{end+1}='end;';
str{end+1}='clear ZZi ZZFNAMES';

str=strs2sray(str);
evalin('caller',str)

function[]=use_test

load solomon
use solomon
reporttest('USE',aresame(x,solomon.x))

function[row]=strs2sray(x)
%STRS2SRAY  Converts a cell array of strings into a string array /w returns
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details      
  
  
if ~iscell(x)
  xtemp=x;
  clear x
  x{1}=xtemp;
end

M=length(x);
for i=1:M
    n(i)=length(x{i});
end
N=max(n);

row=[];

for i=1:M
   row=[row,x{i},char(10)]; 
end


