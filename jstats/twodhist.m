function[mat,xmid,ymid]=twodhist(varargin)
%TWODHIST  Two-dimensional histogram.
%   __________________________________________________________________
%
%   *|* twodhist.png --- Histogram of the global drifter dataset.  
%   Type 'jhelp twodhist' to view this image. *|*
%   __________________________________________________________________
%
%   MAT=TWODHIST(X,Y,XBIN,YBIN) where X and Y are arrays of the same
%   length, creates a two-dimensional histogram MAT with bin edges
%   specified by XBIN and YBIN. 
%
%   If XBIN and YBIN are length N and M, respectively, then MAT is of
%   size M-1 x N-1.  XBIN and YBIN must be monotonically increasing. 
%
%   [MAT,XMID,YMID]=TWODHIST(...) optionally returns the midpoints XMID
%   and YMID of the bins.
%
%   TWODHIST, TWODSTATS, and TWODMED are three related functions for 
%   computing statistics as a function two variables using very fast
%   algorithms that avoid any loops through efficient use of indexing.
%
%   X and Y can also be cell arrays of numerical arrays, in which case 
%   all data values are concatented prior to finding the histogram.
%   __________________________________________________________________
%
%   Automatic bin calculation
%
%   TWODHIST can compute appropriate bins internally.
%
%   [MAT,XMID,YMID]=TWODIST(X,Y,N) uses N bins in the X and Y directions,
%   linearly spaced between the minimum and maximum values, and returns the
%   bin midpoints in XMID and YMID.  MAT is N-1 x N-1.
%
%   [MAT,XMID,YMID]=TWODIST(X,Y,[XMIN XMAX],[YMIN YMAX],N) similarly uses N
%   bins, linearly spaced between the designated X and Y values.  
%   __________________________________________________________________
% 
%   Algorithms
%
%   By default, TWODHIST now works with an internal call to Matlab's
%   HISTCOUNTS2 function, available as of Matlab 2015b.  This is much 
%   faster than the previous algorithm.
%
%   If HISTCOUNTS2 is not available, TWODHIST uses loopless algorithm that
%   is in turn much faster than an explicit loop.  TWODHIST(...,'jlab') 
%   uses this algorithm, while TWODHIST(...,'slow') uses the explicit loop.
%   These options are mostly used for testing purposes.
%   __________________________________________________________________
%
%   Parallelization
%
%   TWODHIST(...,'parallel') parallelizes the computation using the fast 
%   algorithm together with SPMD.  This requires that Matlab's Parallel
%   Computing Toolbox be installed.  While TWODHIST is already very fast,
%   parallelization may be useful for extremely large datasets.
%   __________________________________________________________________
%
%   See also TWODMED, TWODSTATS.
%
%   'twodhist --f' generates the sample figure shown above.
%   'twodhist --t' runs some tests.
%
%   Usage: [mat,xmid,ymid]=twodhist(x,y,N);
%          [mat,xmid,ymid]=twodhist(x,y,[xmin xmax],[ymin ymax],N);
%          mat=twodhist(x,y,xbin,ybin);
%          [mat,xmid,ymid]=twodhist(x,y,xbin,ybin,'parallel');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details    

  
%   Additional output
%   
%   [MAT,XMID,YMID,INDEX]=TWODHIST(X,Y,...) also returns INDEX, an array of
%   the same size as X and Y giving the index into the matrix MAT
%   corresponding to each (X,Y) data point. 
%   __________________________________________________________________


if strcmpi(varargin,'--t')
   twodhist_test;return
elseif strcmpi(varargin,'--f')
   type makefigs_twodhist
   makefigs_twodhist;
   return
end

parstr='serial';
str='histcounts2';

xdata=varargin{1};
ydata=varargin{2};

if iscell(xdata)
    [xdata,ydata]=cell2col(xdata,ydata);
end
if ~isreal(xdata)||~isreal(ydata)
    error('X and Y must be real-valued.');
end
if ~aresame(size(xdata),size(ydata))
     error('X and Y should have the same size.')
end
vcolon(xdata,ydata);

bool=isfinite(xdata)&isfinite(ydata);
xdata=xdata(bool);
ydata=ydata(bool);

for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'ser')||strcmpi(varargin{end}(1:3),'par')
            parstr=varargin{end};
        else
            str=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end

if strcmpi(str(1:3),'his')&&(exist('histcounts2')~=2)
    str='jlab';
    disp('HISTCOUNTS2 not found; reverting to former JLAB algorithm.')
end    

if strcmpi(parstr(1:3),'par')
    if exist('parpool')~=2
        disp('Sorry, parallel algorithm requires the Parallel Computing Toolbox.')
        disp('Defaulting to the serial algorithm.')
        str='serial';
    end
end

if length(varargin)==4
    xbin=varargin{3};
    ybin=varargin{4};
elseif length(varargin)==3
    N=varargin{3};
    xbin=linspace(minmin(xdata),maxmax(xdata),N);
    ybin=linspace(minmin(ydata),maxmax(ydata),N);
elseif length(varargin)==5
    N=varargin{5};
    xbin=linspace(varargin{3}(1),varargin{3}(2),N);
    ybin=linspace(varargin{4}(1),varargin{4}(2),N);
end
    
xbin=xbin(:);
ybin=ybin(:);
if any(diff(xbin)<0)
  error('XBIN must be monotonically increasing')
end
if any(diff(ybin)<0)
  error('YBIN must be monotonically increasing')
end

%Exclude points which are obviously outside of the domain
bool=xdata<xbin(end)&xdata>xbin(1)&ydata<ybin(end)&ydata>ybin(1);
xdata=xdata(bool);
ydata=ydata(bool);

if ~isempty(xdata)
    if strcmpi(parstr(1:3),'ser')
        mat=twodhist_one(xdata,ydata,xbin,ybin,str);
    elseif strcmpi(parstr(1:3),'par')
        %Parallel Algorithm
        disp('TWODHIST employing parallel algorithm.')
        pool=gcp;
        Nworkers=pool.NumWorkers;
        N=length(xdata);
        M=floor(N/Nworkers);
        spmd 
            %Determine the data to send to each worker
            if labindex<Nworkers
                spmdindex=(labindex-1)*M+1:labindex*M;
            else
                %A few leftover time series go to the last worker
                spmdindex=(labindex-1)*M+1:N;
            end
            mati=twodhist_one(xdata(spmdindex),ydata(spmdindex),xbin,ybin,str); 
        end
        mat=zeros(size(mati{1}));
        for i=1:length(mati)
            mat=mat+mati{i};
        end
    end
else
    disp('Warning: No valid data in specified region.')
    mat=0*oprod(ybin(1:end-1),xbin(1:end-1));
end

if nargout>1
  xmid=(xbin+vshift(xbin,1,1))./2;
  xmid=xmid(1:end-1);
end
if nargout>2
  ymid=(ybin+vshift(ybin,1,1))./2;
  ymid=ymid(1:end-1);
end

function[mat]=twodhist_one(xdata,ydata,xbin,ybin,str)

if strcmpi(str(1:3),'his')
    mat=histcounts2(xdata,ydata,xbin,ybin)';
elseif strcmpi(str(1:3),'jla')
    mat=twodhist_jlab(xdata,ydata,xbin,ybin);
elseif strcmpi(str(1:3),'slo')
    mat=twodhist_slow(xdata,ydata,xbin,ybin);
end

function[mat]=twodhist_jlab(xdata,ydata,xbin,ybin)

[xnum,xi,xmid]=bindata(xbin,xdata);
[ynum,yi,ymid]=bindata(ybin,ydata);

mat=zeros([length(ybin)-1,length(xbin)-1]);
index=nan*zeros(size(xdata));

nani=(~isnan(xnum)&~isnan(ynum));

if sum(nani(:))>0
    index(nani)=sub2ind([length(ybin)-1,length(xbin)-1],ynum(nani),xnum(nani));
    [indexsorted,sorter]=sort(index(nani));
    
    [L,ia]=blocklen(indexsorted);
    mat(indexsorted(ia))=L(ia);
end



function[mat]=twodhist_slow(xdata,ydata,xbin,ybin)
mat=zeros(length(ybin),length(xbin));
[xbinb,ybinb]=vshift(xbin,ybin,1,1);
for i=1:length(xbin)
    for j=1:length(ybin)
        mat(j,i)=length(find(xdata>xbin(i)&xdata<=xbinb(i)&...
            ydata>ybin(j)&ydata<=ybinb(j)));
    end
end
mat=mat(1:end-1,:);
mat=mat(:,1:end-1);

function[]=twodhist_test
L=10000;
xdata=3*abs(rand(L,1));
ydata=3*abs(rand(L,1));
xbin=(0:.1:2);
ybin=(0:.2:2);
tic;
mat1=twodhist(xdata,ydata,xbin,ybin);
dt1=toc;
tic
mat2=twodhist(xdata,ydata,xbin,ybin,'slow');
dt2=toc;
bool=aresame(mat1,mat2,1e-10);
reporttest('TWODHIST fast vs. slow algorithm',bool)
%disp(['TWODHIST fast algorithm was ' num2str(dt2./dt1) ' times faster than direct algorithm.'])

if exist('parpool')==2
    mat3=twodhist(xdata,ydata,xbin,ybin,'parallel');
    bool=aresame(mat1,mat3,1e-10);
    reporttest('TWODHIST fast vs. parallel algorithm',bool)
end

xdata=-3*abs(rand(L,1));
ydata=-3*abs(rand(L,1));
xbin=(-2:.1:0);
ybin=(-2:.2:0);
mat1=twodhist(xdata,ydata,xbin,ybin);
mat2=twodhist(xdata,ydata,xbin,ybin,'slow');
bool=aresame(mat1,mat2,1e-10);
reporttest('TWODHIST fast vs. slow algorithm, negative bins',bool)

xdata=randn(L,1);
ydata=randn(L,1);
xbin=(-2:.1:2);
ybin=(-2:.2:2);
mat1=twodhist(xdata,ydata,xbin,ybin);
mat2=twodhist(xdata,ydata,xbin,ybin,'slow');
bool=aresame(mat1,mat2,1e-10);
reporttest('TWODHIST fast vs. slow algorithm, crossing zero',bool)


L=10000;
xdata=3*abs(rand(L,1));
ydata=3*abs(rand(L,1));
xbin=(0:.1:2);
ybin=(0:.2:2);
tic;
mat1=twodhist(xdata,ydata,xbin,ybin);
dt1=toc;
tic
mat2=twodhist(xdata,ydata,xbin,ybin,'slow');
dt2=toc;
bool=aresame(mat1,mat2,1e-10);
reporttest('TWODHIST fast vs. slow algorithm',bool)

if exist('histcounts2')==2
    xdata=randn(L,1);
    ydata=randn(L,1);
    xbin=(-2:.1:2);
    ybin=(-2:.2:2);
    tic;mat1=twodhist(xdata,ydata,xbin,ybin,'histcounts2');toc;
    tic;mat2=twodhist(xdata,ydata,xbin,ybin,'jlab');toc;
    bool=aresame(mat1,mat2,1e-10);
    reporttest('TWODHIST JLAB vs. HISTCOUNTS algorithm, crossing zero',bool)
end

