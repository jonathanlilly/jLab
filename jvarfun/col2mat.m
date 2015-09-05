function[varargout]=col2mat(varargin)
%COL2MAT  Expands 'column-appended' data into a matrix.
%
%   M=COL2MAT(C), where C is a column vector of data segments
%   separated by NANs, expands the data into a matrix such that each
%   segment has its own column.  For instance, the intervals may
%   represent CTD casts of various depths, and the NANs mark
%   transitions between casts.
%
%   M will have enough rows to fit the longest data segment. Segments
%   are read into columns of M from the top down, leaving empty spaces
%   at the bottom which are filled in with NANs.
%
%   If the data is complex, then gaps are filled with NAN+SQRT(-1)*NAN;
%
%   [M1,M2,....]=COL2MAT(C1,C2,...), where C1,C2,... are input column
%   vectors, also works, as does [M1,M2,....]=COL2MAT(MAT), where MAT
%   is a matrix of column vectors.  In both cases the first column is
%   used as the reference for the others, so that only the locations
%   of NANs in the first column matters.
%
%   COL2MAT, MAT2COL, and COLBREAKS together form a system for moving
%   data with segments of nonuniform length rapidly back and forth
%   between a column format and a padded-matrix format.  CTD or float
%   data, for instance, can be stored in the (usually much smaller)
%   column format and converted into the matrix format upon loading.
%
%   COL2MAT(C1,C2,...); with no output arguments overwrites the
%   original input variables.
%
%   See also MAT2COL, COLBREAKS, COL2CELL, CELL2COL.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2000--2013 J.M. Lilly --- type 'help jlab_license' for details
  
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


