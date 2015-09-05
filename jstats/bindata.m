function[xnum,xi,xmid]=bindata(xbin,xdata,str)
%BINDATA  Rapidly sort data into adjacent bins.
%
%   BINDATA is a low-level function called by the TWOD* functions.
%
%   As of Matlab 2015b, note that this function is no longer recommended.
%   Matlab's HISTCOUNTS2 accomplishes the same thing but vastly faster. 
%   BINDATA is included for compatibility with older versions of Matlab. 
% 
%   [NUM,XI]=BINDATA(BIN,X) sorts the array X by its values into adjactent
%   bins, with bin edges given by BIN. 
%
%   NUM and XI are arrays of the same size as X.  NUM gives the bin 
%   number to which the corresponding elements of X are assigned, and 
%   XI is the midpoint value of that bin.
%
%   [NUM,XI,MID]=BINDATA(BIN,X) also returns an array MID of the bin
%   midpoints, which has length LENGTH(BIN)-1.
%  
%   BINDATA uses a clever algorithm to run amazingly fast. 
% 
%   BINDATA is useful for making histograms of large datasets.
%
%   'bindata --t' runs a test.
%
%   See also TWODHIST, TWODSTATS, TWODMED.
%
%   Usage: [num,xi]=bindata(bin,x);
%          [num,xi,xmid]=bindata(bin,x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(xbin, '--t')
    bindata_test,return
end

xbin=sort(xbin(:));
xdata=xdata(:);

if nargin<3
    if length(xdata)>1e6
        str='fast';
    else
        str='slow';
    end
end

xmid=(xbin+vshift(xbin,1,1))./2;
xmid=xmid(1:end-1);

%Turns out the slow algorithm is actually faster for small # data points
if ~isempty(strfind(str,'fast'))
    xnum=bindata_fast(xbin,xdata,xmid);
elseif ~isempty(strfind(str,'slow'))
    xnum=bindata_slow(xbin,xdata,xmid);
end


if nargout>1
    xi=nan*xdata;
    index=find(~isnan(xnum));
    if ~isempty(index)
        xi(index)=xmid(xnum(index));
    end
end

function [xnum]=bindata_fast(xbin,xdata,xmid)

if (length(xbin)==2)&&(minmin(xdata)>xbin(1))&&(maxmax(xdata)<xbin(2))
    %Special case for only one bin
    xnum=ones(size(xdata));
    xnum(isnan(xdata))=nan;
else    
    [xsort,index]=sort(xdata);
    xnum=0*xsort;

    ilast=1;
    for i=1:length(xbin)
        %   Find first
        i1=find(xsort(ilast:end)>=xbin(i),1)+(ilast-1);
        if ~isempty(i1)
            xnum(i1)=xnum(i1)+1;
            ilast=i1;
        end
    end
    xnum=cumsum(xnum);
    vswap(xnum,0,nan);
    vswap(xnum,length(xbin),nan);
    xnum(index)=xnum;  %This puts it back where it came from
end


function [xnum]=bindata_slow(xbin,xdata,xmid)

if  (length(xbin)==2)&&(minmin(xdata)>xbin(1))&&(maxmax(xdata)<xbin(2))
    %Special case for only one bin
    xnum=ones(size(xdata));
    xnum(isnan(xdata))=nan;
else
    xnum=nan*xdata;
    for i=1:length(xdata)
        for j=2:length(xbin)
            if xdata(i)>=xbin(j-1)&&xdata(i)<xbin(j)
                xnum(i)=j-1;
            end
        end
    end
end


function[]=bindata_test


xbin=[0 1 2 2.5 3]';
xmid1=[0.5 1.5 2.25 2.75]';
xdata=[-1   .1 .1 1.5 2.75 4]';
xnum1=  [nan  1  1  2   4   nan]';
xi1=  [nan  0.5 0.5 1.5 2.75 nan]';

%xnum1=  [nan   1  2   4   nan]';

[xnum2,xi2,xmid2]=bindata(xbin,xdata);

reporttest('bindata simple',aresame(xnum1,xnum2)&&aresame(xmid1,xmid2)&&aresame(xi1,xi2))


xdata=randn(100000,1);
xbin=(-3:.25:3);
tic;[xnum1,xi1,xmid1]=bindata(xbin,xdata,'fast');t1=toc;
tic;[xnum2,xi2,xmid2]=bindata(xbin,xdata,'slow');t2=toc;


reporttest('bindata slow vs. fast algorithm',aresame(xnum1,xnum2)&&aresame(xmid1,xmid2)&&aresame(xi1,xi2))
%disp(['Fast algorithm was ' num2str(t2./t1) ' times faster.'])
