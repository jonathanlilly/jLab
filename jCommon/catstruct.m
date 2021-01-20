function[struct]=catstruct(varargin)
%CATSTRUCT  Concatenates the array elements of a cell array of structures.
%
%   Let C be a cell array of M structures having N identical fields, where 
%   all of the fields are themselves arrays,
%         
%      C{1}.F1=A11, C{1}.F2=A12, ... , C{1}.FN=A1N
%      C{2}.F1=A21, C{2}.F2=A22, ... , C{2}.FN=A2N
%         :        :        :        :        :
%      C{M}.F1=AM1, C{2}.F2=AM2, ... , C{M}.FN=AMN
%   
%    and where all the AMN for any particular field (F1,say), have the same
%    number of columns, that is, SIZE(A11,2)=SIZE(A21,2)= ... =SIZE(AM1,2).
%
%    S=CATSTRUCT(C) then returns a structure in which all fields of the
%    same type have been concatenated:
%
%      SNEW.F1=[A11;   SNEW.F2=[A12;    ...    SNEW.FN=[A1N;
%               A21;            A22;    ...             A2N;
%                :               :      ...              :
%               AM1]            AM2]    ...             AMN];
%
%   Note S is a single structure, whereas C is a cell array of structures.
%
%   As a simple example, with C{1}.F=[1;2] and C{2}.F=[3;4], S=catstruct(C)
%   returns S with a field F such that S.F=[1;2;3;4].
%
%   One use of CATSTRUCT is in concatentating output from multiple cores
%   when using SPMD in the parallel computing toolbox. 
%
%   SNEW=CATSTRUCT(S,INDEX) only concatenates those fields listed in INDEX.
%   Here index is a list of numbers indicating the fields by the order in 
%   which they appear, e.g., INDEX=1:10 concatenates the first 10 fields.
%
%   Usage:  snew=catstruct(s);
%           snew=catstruct(s,index);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2003--2021 J.M. Lilly --- type 'help jlab_license' for details


structo=varargin{1};

fx=fields(structo{1});
if nargin==1
    index=1:length(fx);
else
    index=varargin{2};
end
for j=1:length(fx)
    struct.(fx{j})=[];
end

%   CATSTRUCT also works if the AMN are all cell arrays rather than 
%   numerical arrays, with the output fields becoming concatenated cell
%   cell arrays.

for i=1:length(structo)
     structi=structo{i};
     for j=1:length(index)
%         struct.(fx{j})=[struct.(fx{j});structi.(fx{j})];
         struct.(fx{index(j)})=[struct.(fx{index(j)});structi.(fx{index(j)})];
     end
end
        
% z=varargin{1};
% for i=2:nargin
%     z=catstruct1(z,varargin{i});
% end
% 
% function[z]=catstruct1(x,y)
% 
% fx=fields(x);
% fy=fields(y);
% 
% for k=1:length(fx)
%   if ~aresame(fx{k},fy{k})
%     error('X and Y must have identical fields')
%   end
% end
% 
% z=[];
% for i=1:length(fx);
%    vx=x.(fx{i});
%    vy=y.(fx{i});
%    %vx=getfield(x,fx{i});
%    %vy=getfield(y,fx{i});
%    N=max(size(vx,1),size(vy,1));
%  
%    iscell(vx)
%    if all(isreal(vx))  
%       vz=nan*zeros(N,size(vx,2)+size(vy,2));
%    else
%       vz=(nan+sqrt(-1)*nan)*zeros(N,size(vx,2)+size(vy,2));
%    end
%    vz(1:size(vx,1),1:size(vx,2))=vx;
%    vz(1:size(vy,1),(1:size(vy,2))+size(vx,2))=vy;
%    %z=setfield(z,fx{i},vz);
%    z.(fx{i})=vz;
% end
% 
%      
