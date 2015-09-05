function[varargout]=cellchunk(varargin)
%CELLCHUNK  Converts cell array data into uniform length 'chunks'.
%
%   Y=CELLCHUNK(X,L), where X is a cell array of variable-length numerical
%   arrays, extracts non-overlapping 'chunks' of the uniform length L and
%   returns these in Y, a matrix having L rows.
%   
%   Successive chunks of length L data in X are become successive columns
%   in Y, regardless of the cell to which they initially belonged.  
%   Residuals of length less then L at the end of each cell are discarded.
%
%   For example:
%
%       x{1}=[1 2 3 4]'; x{2}=[5 6 7]'; x{3}=[8 9]'; 
%
%       cellchunk(x,2) =  [1 3 5 8]
%                         [2 4 6 9]
%
%   This is useful, for example, in examining spectra from uniform-length 
%   intervals in Lagragnian data such as FLOATS.MAT or DRIFTERS.MAT.
%
%   If there is a remainder R in how many times L goes into the length of
%   the time series, FLOOR(R/2) points will be thrown away from the
%   beginning of the time series and CEIL(R/2) points from the end. 
%
%   [Y1,Y2,...,YN]=CELLCHUNK(X1,X2,...,XN,L) also works, where the XN are
%   all cell array of the same size.
% 
%   CELLCHUNK also works when the XN are simply column arrays.
%
%   CELLCHUNK with no output arguments overwrite the original named input
%   variables. 
%   __________________________________________________________________
%   
%   Overlap
%
%   CELLCHUNK(...,L,'overlap') instead outputs chunks of length L with a 
%   50% overlap.  That is, successive columns of the output overlap by L/2.
%
%   In this case, remainders are discarded from the ends of the time 
%   series, instead of being split between the beginning and the end as in
%   the non-overlapping case.  
%   __________________________________________________________________   
%
%   Parallelization
%
%   CELLCHUNK(...,'parallel') parallelizes the computation using a PARFOR 
%   loop.  This requires that Matlab's Parallel Computing Toolbox be 
%   installed, and is useful for very large datasets.
%   __________________________________________________________________
%
%   See also TRAJCHUNK.
%
%   'cellchunk --t' runs a test.
%
%   Usage: y=cellchunk(x,L);
%          [y1,y2,y3]=cellchunk(x1,x2,x3,L);
%          [y1,y2,y3]=cellchunk(x1,x2,x3,L,'overlap');
%          [y1,y2,y3]=cellchunk(x1,x2,x3,L,'overlap','parallel');
%          cellchunk(x1,x2,x3,L);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    cellchunk_test,return
end


str='nooverlap';
cores='serial';

for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'ser')|| strcmpi(varargin{end}(1:3),'par')
            cores=varargin{end};
        else
            str=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end


L=varargin{end};
varargin=varargin(1:end-1);
na=length(varargin);
for j=1:na
    if strcmpi(str(1:3),'ove')
        varargout{j}=cellchunk_overlap(L,varargin{j},cores);
    else
        varargout{j}=cellchunk_nooverlap(L,varargin{j},cores);
    end
end

eval(to_overwrite(na));

function[x1,x2]=cellchunk_nooverlap(L,x,cores)

%Need to use a boolean here because parfor does not see the nargout
if nargout==2
    bool=true;
else
    bool=false;
end

if ~iscell(x);
    temp{1}=x;
    x=temp;
end

if ~strcmpi(cores(1:3),'par')
    for i=1:length(x)
        N=length(x{i});
        M=mod(N,L);
        index=1+floor(M/2):N-ceil(M/2);
        x1{i}=vindex(x{i},index,1);
        if bool
            index=1+floor((M+L)/2):N-ceil((M+L)/2);
            x2{i}=vindex(x{i},index,1);
        end
    end
else  %Same exact thing but with parfor loop
    parfor i=1:length(x)
        N=length(x{i});
        M=mod(N,L);
        index=1+floor(M/2):N-ceil(M/2);
        x1{i}=vindex(x{i},index,1);
        if bool
            index=1+floor((M+L)/2):N-ceil((M+L)/2);
            x2{i}=vindex(x{i},index,1);
        end
    end
end


x1=nonnan(cell2col(x1));
x1=reshape(x1,L,length(x1)/L);
if nargout==2
    x2=nonnan(cell2col(x2));
    x2=reshape(x2,L,length(x2)/L);
end

function[x]=cellchunk_overlap(L,x,cores)

if ~iscell(x);
    temp{1}=x;
    x=temp;
end

xcol=cell2col(x);
num=[1:length(xcol)]';
num(isnan(xcol))=nan;
col2cell(num);

[num1,num2]=cellchunk_nooverlap(L,num,cores);
[y1,y2]=cellchunk_nooverlap(L,x,cores);

x=[y1 y2];
num=[num1(1,:) num2(1,:)];
[sorted,sorter]=sort(num);
x=x(:,sorter);

function[]=cellchunk_test
 
load ebasnfloats
use ebasnfloats

[numc,latc,lonc]=cellchunk(num(3),lat(3),lon(3),200);
dnum=numc(1,2:end)-numc(end,1:end-1);%Should differ by 1 at L
%figure,plot(numc)
reporttest('CELLCHUNK no overlap',allall(dnum==1))


[numc,latc,lonc]=cellchunk(num(3),lat(3),lon(3),200,'overlap');
numc(1,2:end)-numc(100,1:end-1);%Should differ by 1 at L/2
%figure,plot(numc)
reporttest('CELLCHUNK with 50% overlap',allall(dnum==1))

