function[]=reporttest(str,bool)
%REPORTTEST  Reports the result of an m-file function auto-test.
%
%   Called by JLAB_RUNTESTS.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2003--2014 J.M. Lilly --- type 'help jlab_license' for details


global BOOL_JLAB_RUNTEST
global FUNCTION_NUMBER
global FUNCTION_NUMBER_ARRAY

FUNCTION_NUMBER_ARRAY=[FUNCTION_NUMBER_ARRAY;FUNCTION_NUMBER];

if bool
    disp([str ' test: passed'])
    BOOL_JLAB_RUNTEST=[BOOL_JLAB_RUNTEST;1];
else  
    disp([str ' test: FAILED'])
    BOOL_JLAB_RUNTEST=[BOOL_JLAB_RUNTEST;0];
end
