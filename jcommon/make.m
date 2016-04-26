function[evalme]=make(varargin)
%MAKE  Create a structure containing named variables as fields.
%
%   MAKE NAME V1 V2 ... VN where the VI are names of variables, creates
%   a structure NAME with fields
%
%       NAME.V1, NAME.V2, ... NAME.VN
%  
%   The original variables are not deleted.  If structure NAME has any
%   existing fields, the new fields are appended after the existing fields.
%
%   MAKE(CELL), where CELL is a cell array of strings such that 
%   CELL{1}=NAME, CELL{2}=V1, ... CELL{N}=VN, also works.
%  
%   This is useful for handling multiple datasets with the same variable 
%   names.  The structures can be then kept in memory and 'mapped' into 
%   variables as needed using USE. 
%
%   See also USE, MATSAVE.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2016 J.M. Lilly --- type 'help jlab_license' for details        
 

nargs=length(varargin);

if iscell(varargin{1}) && nargin==1
    args=varargin{1};
else
    args=varargin;
end

name=args{1};
nargs=length(args);

i=1;
evalme='';
while i<nargs
    i=i+1;
    nameB=args{i};
    nameA=[name '.' nameB];
    evalme=[evalme,nameA,'=',nameB,';'];
    evalme=[evalme,char(10)];
end

if nargout==0
   evalin('caller',evalme)
   clear evalme
end

%function[]=make_test
%Do later




