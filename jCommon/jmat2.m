function[x]=jmat2(theta)
%JMAT2 2x2 rotation matrix through specified angle.
%
%   J=JMAT2(PHI) creates a rotation matrix 
%
%       J=[COS(PHI) -SIN(PHI);
%          SIN(PHI)  COS(PHI)]
%
%   such that J*X rotates the column-vector X by PHI radians
%   counterclockwise.
%
%   If LENGTH(PHI)>1, then J will have dimension  2 x 2 x SIZE(PHI);
%
%   See also VECTMULT, JMAT3, and TMAT.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details        



x=zeros(2,2,size(theta,1),size(theta,2),size(theta,3),size(theta,4));
x(1,1,:,:,:,:)=cos(theta);
x(1,2,:,:,:,:)=-sin(theta);
x(2,1,:,:,:,:)=sin(theta);
x(2,2,:,:,:,:)=cos(theta);

% x=zeros(2,2,length(theta(:)));
% x(1,1,:)=cos(theta(:));
% x(1,2,:)=-sin(theta(:));
% x(2,1,:)=sin(theta(:));
% x(2,2,:)=cos(theta(:));
% 
% 
% if size(theta,1)~=length(theta(:));
%     x=reshape(x,[2,2,size(theta)]);
% end
