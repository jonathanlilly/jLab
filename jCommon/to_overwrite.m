function[str]=to_overwrite(N)
%TO_OVERWRITE Returns a string to overwrite original arguments.
%
%   STR=TO_OVERWRITE(N), when called from within an m-file which has
%   VARARGIN for the input variable, returns a string which upon
%   EVAL(STR) will cause the first N input variables in the caller
%   workspace with the values contained in the first N elements of
%   VARARGOUT.
%
%   See also TO_GRAB_FROM_CALLER.
%
%   Usage: eval(to_overwrite(N))
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2006 J.M. Lilly --- type 'help jlab_license' for details      

str{1}=     'if nargout==0';
str{end+1}= '   global ZZoutput';
str{end+1}= '   evalin(''caller'',[''global ZZoutput''])';
str{end+1}=['   for i=1:' int2str(N)];
str{end+1}= '     if ~isempty(inputname(i))';
str{end+1}= '       ZZoutput=varargout{i};';
str{end+1}= '       assignin(''caller'',inputname(i), ZZoutput)';
str{end+1}= '     end';
str{end+1}= '   end';
str{end+1}='   evalin(''caller'',[''clear ZZoutput''])';
str{end+1}='end';

str=strs2row(str);


function[row]=strs2row(x)
%STRS2ROW  Converts a cell array of strings into a row array

M=length(x);
for i=1:M
    n(i)=length(x{i});
end
N=max(n);

row=[];

for i=1:M
    row=[row,char(10),x{i}]; 
end
