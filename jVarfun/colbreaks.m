function[varargout]=colbreaks(varargin)
%COLBREAKS  Insert NANs into discontinuties in a vector.
%
%   O1=COLBREAKS(I1), where I1 is a column vector, inserts NANs into I1 at
%   its discontinuities (and puts one at the end), for example:
%		
%	  COLBREAKS([2;2;2;3;3])=[2;2;2;NAN;3;3;NAN];
%
%   [O1,O2,O3,...]=COLBREAKS(I1,I2,...), where I1,I2,... are all column 
%   vectors of the same length, uses I1 as a reference for the other 
%   vectors and inserts NANs into all the vectors at locations where I1 is 
%   discontinous.
%	
%   For instance, using station number for I1 will sort a column of CTD 
%   data into a NAN-padded matrix.
%
%   Note that for complex-valued input arrays, a complex-valued NAN, 
%   NAN+SQRT(-1)*NAN, is used to indicate the break locations. 
%
%   MAT2=COLBREAKS(MAT1), where MAT1 and MAT2 are matrices of the same 
%   size, also works.  In this case the first column of MAT1 is used as the
%   reference vector.
%
%   COLBREAKS, COL2MAT, and MAT2COL together form a system for moving data
%   with segments of nonuniform length rapidly back and forth between a 
%   column format and a padded-matrix format. CTD or float data, for 
%   instance, can be stored in the (usually much smaller) column format and
%   converted into the matrix format upon loading.
%  
%   COLBREAKS(C1,C2,...); with no output arguments overwrites the original
%   input variables.  
%
%   See also COL2MAT, MAT2COL, COL2CELL, CELL2COL, ORBITBREAKS.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2000--2013 J.M. Lilly --- type 'help jlab_license' for details

if strcmpi(varargin{1},'--t'),return,end
%Look ma, no loops!

bmat=0;
nargs=nargin;
ref=varargin{1};
if nargin==1 && size(ref,2)>1
    temp=ref;
    ref=temp(:,1);
    c=cell(size(temp,2),1);
    for i=1:size(temp,2)-1
        c{i}=temp(:,i+1);
    end
    nargs=size(temp,2);
    clear temp
    bmat=1;
else
    c=varargin(2:end);
end



ii=0;
index=find(diff(ref)~=0)+1;
bool=false(size(ref));
bool(index)=ones(size(index));
index=(1:length(ref))'+cumsum(bool);

while ii<nargs - 1
    ii=ii+1;
    col=c{ii};
    
    if isreal(col)
        colout=nan*ones(max(index),1);
        colout(index)=col;
        colout=[colout;nan];
    else
        colout=(nan+sqrt(-1)*nan)*ones(max(index),1);
        colout(index)=col;
        colout=[colout;(nan+sqrt(-1)*nan)];
    end
    
    varargout{ii+1}=colout;
end

if isreal(ref)
    refout=nan*ones(max(index),1);
    refout(index)=ref;
    refout=[refout;nan];
else
    refout=(nan+sqrt(-1)*nan)*ones(max(index),1);
    refout(index)=ref;
    refout=[refout;(nan+sqrt(-1)*nan)];
end

varargout{1}=refout;
if bmat
    for i=1:nargs-1
        varargout{1}(:,i+1)=varargout{i+1};
    end
    varargout=varargout(1);
end
if nargout==0
    eval(to_overwrite(nargs));
end



