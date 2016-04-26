function[varargout]=col2mat(varargin)
%COL2MAT  Expands 'column-appended' data into a matrix.
%
%   M=COL2MAT(C), where C is a column vector of data segments separated by
%   NANs, returns a matrix M in which each segment has its own column.  
%
%   M will have enough rows to fit the longest data segment. Segments are
%   read into columns of M from the top down, leaving empty spaces at the 
%   bottom which are filled in with NANs.
%
%   If the data is complex, then gaps are filled with NAN+SQRT(-1)*NAN;
%
%   [M1,M2,...,MN]=COL2MAT(C1,C2,...,CN), where the CN are input column
%   vectors, also works, as does [M1,M2,...,MN]=COL2MAT(MAT), where MAT is 
%   a matrix of column vectors.  In both cases the locations of NaNs in the 
%   first column is used as the reference for the others. 
%
%   COL2MAT(C1,C2,...); with no output arguments overwrites the original 
%   input variables.
%
%   COL2MAT, MAT2COL, and COLBREAKS together form a system for moving data
%   with segments of nonuniform length rapidly back and forth between a 
%   column format and a padded-matrix format. 
%
%   Note that while COL2CELL and CELL2COL work for arrays having multiple
%   columns, COL2MAT only works with column vectors. 
%
%   See also MAT2COL, COLBREAKS, COL2CELL, CELL2COL.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2000--2016 J.M. Lilly --- type 'help jlab_license' for details
  
%no loops!

if ischar(varargin{1})
    if strcmpi(varargin{1}(1:3),'--t')
        return
    end
end

if isempty(varargin{1})
   for i=1:nargout
       varargout{i}=[];
   end
   return
end

if nargin>1
	for i=1:nargin
        	eval(['data(:,' int2str(i) ')=varargin{' int2str(i) '};'])
	end
else 
	data=varargin{1};
end
  
for ii=1:size(data,2)
	col=data(:,ii);
	if ii==1;
        %Account for potential missing NaNs at the end
        if ~isnan(col(end))
            col(end+1)=nan;
            data(end+1,:)=nan;
        end
		%use index into NANs and index derivative
		%to determine size of new matix
		nani=find(isnan(col));
		dnani=diff([0;nani]);
		nrows=max(dnani);
		ncols=length(nani);

		%determine index into top of each column ==a
		%and index into first NAN in each column ==b
		a=1+nrows*(0:ncols-1)';
        %size(dnani)
        %size(a)
		b=dnani+a;
		index=zeros(nrows*ncols,1);
		
		b=b(b<length(index));
		%mark all the numbers between each a and each b
		index(a)=1;	
		index(b)=index(b)-1;%this matters (a may equal b)		
		index=find(cumsum(index));
		if length(index)>length(col)
		   index=index(1:length(col));
		end
	end
	mat=nan*ones(nrows,ncols);
    mat(index)=col;
    if ~isreal(mat)
        vswap(mat,nan,nan+sqrt(-1)*nan);
    end
	varargout{ii}=mat;
end



if nargout>nargin
	varargout{int2str(size(data,2)+1)}=index;
end

if nargout==0
  eval(to_overwrite(nargin));
end


