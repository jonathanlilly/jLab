function[]=matsave(varargin)
%MATSAVE  Create and save structure of variables as a mat-file. 
%
%   MATSAVE NAME X1 X2 ... XN  where X1 ... XN are names of variables in 
%   the current workspace, first creates a structure
%
%       NAME.X1, NAME.X2, NAME.X3
%
%   then saves it as a .mat file called NAME, and finally loads this 
%   mat-file into the current workspace.
%
%   By default, the mat-file is saved in a directory called 'Matfiles'.
%   This directory can be anywhere, but its location must be added to the 
%   Matlab search path using ADDPATH, e.g. 'addpath ~/Matfiles'.  
%
%   It is preferable to put the ADDPATH statement in your STARTUP file so 
%   you do not have to re-enter it every session. 
%
%   If a directory called 'Matfiles' is not found, MATSAVE saves the 
%   mat-file to the current directory.  
%
%   By default NAME is a string rather than a variable.  However, NAME can 
%   also be a variable whose value is a string, through the use of the '-n'
%   flag.  Then
%
%        MATSAVE -N VAR X1 ... XN
%
%   with VAR='sample' will create a structure named 'sample'.
%
%   The version of mat-file format to save in may also be specified through
%   the of the -v flag, such that
%
%        MATSAVE -V6 NAME X1 ... XN
%
%   will create a structure NAME and save it as a Version 6 mat-file; see 
%   SAVE for allowable version options.  The default is -v7.3.
%
%   See also MAKE, USE.
%
%   Usage:  matsave name x1 x2 x3
%           matsave -n var x1 x2 x3
%           matsave -v6 -n var x1 x2 x3
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2018 J.M. Lilly --- type 'help jlab_license' for details


version='-v7.3';
vars=varargin;

nflagsin=0;
fflagin=0;

i1=vars{1};
if strcmpi(i1(1),'-')
   if strcmpi(i1(2),'n')
       %-n flag for filename was input
       fflagin=1;
       vars=vars(2:end);
       nflagsin=1;
   elseif strcmpi(i1(2),'v')
       %-v flag for version was input
       version=i1;
       vars=vars(2:end);
       nflagsin=1;
              
       i1=vars{1};
       if strcmpi(i1(2),'n')
          %-n flag for filename was input
          fflagin=1;
          vars=vars(2:end);
          nflagsin=2;
       end
   end
end

eval(to_grab_from_caller(nflagsin+1:nargin))
if fflagin
     vars{1}=eval(vars{1});
end

matsave_name=vars{1};

evalme='make ';
for i=1:length(vars)
    evalme=[evalme vars{i} ' '];
end
evalin('caller',evalme)

dirname=findpath('/Matfiles');

if isempty(dirname)
    dirname=pwd;
end

matsave_fullname=[dirname '/' matsave_name];
evalin('caller',['save(''' matsave_fullname ''',''' version ''',''' matsave_name ''')'])

%evalin('caller',['save( 'matsave_fullname ',' version ',' matsave_name ')'])
%evalin('caller','save(matsave_fullname,version,matsave_name)')
%evalin('caller',['load ' matsave_fullname])


