function[x]=tmat()
%TMAT  2x2 complex grouping matrix.  TMAT = [1  i; 1 -i] / SQRT(2)
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details        

x=ones(2,2);
x(1,2)=sqrt(-1);
x(2,2)=-sqrt(-1);
x=frac(1,sqrt(2)).*x;
