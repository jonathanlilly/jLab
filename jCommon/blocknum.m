function[num,a,b]=blocknum(x,delta)
%BLOCKNUM  Numbers the contiguous blocks of an array.
%
%   Suppose X is a column vector which contains blocks of identical
%   values, e.g. X=[0 0 0 1 1 3 2 2 2]';
%
%   N=BLOCKNUM(X) gives each contiguous block of identical values a
%   unique number, in order beginning with 1.  Each elements of N
%   specifies the number of the block to which the corresponding
%   element of X belongs.
%
%   In the above example, N=[1 1 1 2 2 3 4 4 4]';
%  
%   [N,A,B]=BLOCKNUM(X) also returns arrays A and B which are indices
%   into the first and last point, respectively, of each block.  In
%   the above example, A=[1 4 6 7]' and B=[3 5 6 9]';
%
%   [...]=BLOCKNUM(X,D) defines the junction between two blocks as
%   locations where ABS(DIFF(X))>D.  Thus D=1 is a 'rate of change'
%   definition, and X=[1 2 3 5 6 10 16 17 18]'; will yield the same 
%   result for N as in the previous example.
%
%   See also BLOCKLEN.
%
%   Usage: num=blocknum(x);
%          [num,a,b]=blocknum(x);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2000--2014 J.M. Lilly --- type 'help jlab_license' for details
  
if strcmpi(x,'--t')
  blocknum_test;return;
end

if nargin==1
  delta=0;
end

dx=diff(x);
index=find(abs(dx)>delta);
num=zeros(size(x));
if ~isempty(index)
  num(index+1)=1;
end
num(1)=1;
num=cumsum(num);


if nargout>=2
   a=[1;find(diff(num)~=0)+1];
end
if nargout>=3
  b=[find(diff(num)~=0);length(num)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[]=blocknum_test
x= [0 0 0 1 1 3 2 2 2]';
y= [1 1 1 2 2 3 4 4 4]';

reporttest('BLOCKNUM with D==0',all(y==blocknum(x)))

x=[1 2 3 5 6 10 16 17 18]';
reporttest('BLOCKNUM with D==1',all(y==blocknum(x,1)))

