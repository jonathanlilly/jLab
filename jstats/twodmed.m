function[mat,xmid,ymid]=twodmed(varargin)
%TWODMED  Median value of a function of two variables.
%   __________________________________________________________________
%
%   *|* twodmed.png --- Median speed from the global drifter dataset.  
%   Type 'jhelp twodmed' to view this image. *|*
%   __________________________________________________________________
%
%   MED=TWODMED(X,Y,Z,XBIN,YBIN) where X, Y and Z are arrays of the same
%   length, forms the median of Z over the XY plane.
%
%   If XBIN and YBIN are length N and M, respectively, then MED is of 
%   size M-1 x N-1.  Bins with no data are assigned a value of NAN.
%
%   XBIN and YBIN must be monotonically increasing. 
%
%   MED=TWODMED(X,Y,Z,N) uses N bins in the X and Y directions, linearly
%   spaced between the minimum and maximum values.  MED is N-1 x N-1.
%
%   MED=TWODMED(X,Y,Z,[XMIN XMAX],[YMIN YMAX],N) uses N bins, linearly
%   spaced between the designated X and Y values.  MED is N-1 x N-1. 
%
%   [MED,XMID,YMID]=TWODMED(...) optionally returns the midpoints XMID
%   and YMID of the bins.
%
%   X, Y, and Z can also be cell arrays of numerical arrays, in which case 
%   all data values are concatented prior to finding the histogram.
%
%   TWODMED, TWODSTATS, and TWODHIST are three related functions for 
%   computing statistics as a function two variables using very fast
%   algorithms that avoid any loops through efficient use of indexing. 
%   __________________________________________________________________
%
%   Algorithms
%
%   TWODMED uses a fast (exact) algorithm which is particularly efficient 
%   for large arrays.
%
%   The values of Z are sorted into bins according to the associated (X,Y)
%   value, with bin edges specified by XBIN and YBIN, and the median of 
%   finite values of Z in each bin is computed looplessly using indexing.
%
%   TWODMED(...,'slow') uses a slow algorithm which uses less memory.  
%   By default, TWODMED uses a fast but memory-intensive algorithm.  
%   Use this flag if you get an out-of-memory error.
%   __________________________________________________________________
%
%   See also TWODHIST, TWODSTATS.
%   
%   'twodmed --f' generates the sample figure shown above.
%   'twodmed --t' runs a test.
%
%   Usage: med=twodmed(x,y,z,N);
%          med=twodmed(x,y,z,[xmin xmax],[xmin xmax],N);
%          med=twodmed(x,y,z,xbin,ybin);
%          [med,xmid,ymid]=twodmed(x,y,z,xbin,ybin);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2015 J.M. Lilly --- type 'help jlab_license' for details    

 
if strcmpi(varargin,'--t')
   twodmed_test;return
elseif strcmpi(varargin,'--f')
   type makefigs_twodmed
   makefigs_twodmed;
   return
end

str='fast';
xdata=varargin{1};
ydata=varargin{2};
zdata=varargin{3};

if iscell(xdata)
    [xdata,ydata,zdata]=cell2col(xdata,ydata,zdata);
end
if ~isreal(xdata)||~isreal(ydata)||~isreal(zdata)
    error('X, Y, and Z must be real-valued.');
end
vcolon(xdata,ydata,zdata);
if ~aresame(size(xdata),size(ydata))||~aresame(size(xdata),size(zdata))
     error('X, Y, and Z should have the same number of points.')
end
bool=isfinite(xdata)&isfinite(ydata)&isfinite(zdata);
xdata=xdata(bool);
ydata=ydata(bool);
zdata=zdata(bool);

if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
end

if length(varargin)==5
    xbin=varargin{4};
    ybin=varargin{5};
elseif length(varargin)==4
    N=varargin{4};
    xbin=linspace(minmin(xdata),maxmax(xdata),N);
    ybin=linspace(minmin(ydata),maxmax(ydata),N);
elseif length(varargin)==6
    N=varargin{6};
    xbin=linspace(varargin{4}(1),varargin{4}(2),N);
    ybin=linspace(varargin{5}(1),varargin{5}(2),N);
end
     
if nargin==5
  str='fast';
end

vcolon(xbin,ybin);
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
zdata=zdata(bool);

if ~isempty(zdata)
    if ~isempty(strfind(str,'fast'))
      mat=twodmed_fast(xdata,ydata,zdata,xbin,ybin);
    else
      mat=twodmed_slow(xdata,ydata,zdata,xbin,ybin);
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

%vsize(xdata,ydata,xbin,ybin)
hist=twodhist(xdata,ydata,xbin,ybin);
index=find(hist==0);
if ~isempty(index)
    mat(index)=nan;
end

function[mat]=twodmed_fast(xdata,ydata,zdata,xbin,ybin)

[xnum,xi,xmid]=bindata(xbin,xdata);
[ynum,yi,ymid]=bindata(ybin,ydata);

[mat,matz]=vzeros([length(ybin)-1,length(xbin)-1],'nan');

nani=find(~isnan(xnum)&~isnan(ynum));

if ~isempty(nani)
    index=sub2ind([length(ybin)-1,length(xbin)-1],ynum(nani),xnum(nani));
    [index,sorter]=sort(index);
    
    index=nonnan(index);
    if ~isempty(index)
        [L,ia,ib,numblock]=blocklen(index);
        L=L(ia);

        matz=zdata(nani);
        matz=matz(sorter);

        colbreaks(numblock,matz);
        col2mat(numblock,matz);
        matz=sort(matz,1);

        medz=zeros(size(ia));

        colindex=find(isodd(L));
        if ~isempty(colindex)
            ijindex=sub2ind(size(matz),(L(colindex)+1)/2,colindex);
            medz(colindex)=matz(ijindex);
        end

        colindex=find(iseven(L));
        if ~isempty(colindex)
            ijindex1=sub2ind(size(matz),L(colindex)/2,colindex);
            ijindex2=sub2ind(size(matz),L(colindex)/2+1,colindex);
            medz(colindex)=matz(ijindex1)/2+matz(ijindex2)/2;
        end
        mat(index(ia))=medz;
    end
end


function[mat]=twodmed_slow(xdata,ydata,zdata,xbin,ybin)
vcolon(xdata,ydata,zdata,xbin,ybin);
index=find(isfinite(xdata)&isfinite(ydata)&isfinite(zdata));
vindex(xdata,ydata,zdata,index,1);

mat=0*oprod(ybin,xbin); 
[xbinb,ybinb]=vshift(xbin,ybin,1,1);
for i=1:length(xbin)
   for j=1:length(ybin)
         index=find(xdata>xbin(i)&xdata<=xbinb(i)&ydata>ybin(j)&ydata<=ybinb(j));
         if ~isempty(index)
             mat(j,i)=median(zdata(index));
         end
   end
end
mat=mat(1:end-1,:);
mat=mat(:,1:end-1); 

function[]=twodmed_test
L=10000;
xdata=3*abs(rand(L,1));
ydata=3*abs(rand(L,1));
zdata=randn(L,1);
xbin=(0:.01:2);
ybin=(0:.02:2);
tic;
mat1=twodmed(xdata,ydata,zdata,xbin,ybin,'fast');
dt1=toc;
tic;
mat2=twodmed(xdata,ydata,zdata,xbin,ybin,'slow');
dt2=toc;
bool=aresame(mat1,mat2,1e-10);
reporttest('TWODMED fast vs. slow algorithm',bool)
disp(['TWODMED fast algorithm was ' num2str(dt2./dt1) ' times faster than direct algorithm.'])

xdata=-3*abs(rand(L,1));
ydata=-3*abs(rand(L,1));
xbin=(-2:.01:0);
ybin=(-2:.02:0);
mat1=twodmed(xdata,ydata,zdata,xbin,ybin,'fast');
mat2=twodmed(xdata,ydata,zdata,xbin,ybin,'slow');
bool=aresame(mat1,mat2,1e-10);
reporttest('TWODMED fast vs. slow algorithm, negative bins',bool)


