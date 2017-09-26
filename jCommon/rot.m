function[x]=rot(x)
% ROT Complex-valued rotation:  ROT(X)=EXP(SQRT(-1)*X)
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2015 J.M. Lilly --- type 'help jlab_license' for details      

if strcmpi(x,'--t')
  rot_test;clear x;return
end

x=rot_mainangle(x);
index1=find(x==pi/2);
index2=find(x==-pi/2);
index3=find(x==pi|x==-pi);
index4=find(x==0|x==2*pi|x==-2*pi);

x=exp(sqrt(-1).*x);

%It is unfortunately necessary to specify some special cases
%This is because in Matlab, exp(sqrt(-1))==sqrt(-1) is false

if ~isempty(index1)
   x(index1)=sqrt(-1);
end
if ~isempty(index2)
   x(index2)=-sqrt(-1);
end
if ~isempty(index3)
   x(index3)=-1;
end
if ~isempty(index4)
   x(index4)=1;
end


function[x]=rot_mainangle(x)
%Put X between -pi and pi
x=mod(x,2*pi);
index=find(x>pi);
if ~isempty(index)
   x(index)=x(index)-2*pi;
end

function[]=rot_test
x=[pi/2      -pi/2     pi -pi 2*pi -2*pi 0];
y1=[sqrt(-1) -sqrt(-1) -1 -1  1     1    1];
reporttest('ROT handles special cases at n*pi/2',aresame(rot(x),y1))
x=(-3*pi:.1:3*pi);
reporttest('ROT inverts angle',aresame(angle(rot(x)),rot_mainangle(x),1e-11))
%plot(angle(rot(x)),'b'),hold on
%plot(rot_mainangle(x),'r--')
%plot(abs(rot_mainangle(x)-angle(rot(x))))
