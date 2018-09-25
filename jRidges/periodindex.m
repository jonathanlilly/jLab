function[index]=periodindex(varargin)
%PERIODINDEX  Returns time index in increments of instantaneous period.
%
%   INDEX=PERIODINDEX(DT,OMEGA,N) returns a time index INDEX that skips 
%   every N times the instantaneous period 2*pi/OMEGA.
%
%   OMEGA is a column vector of instantaneous frequency in radians per unit 
%   time as computed by INSTMOM. DT is the sample time, a scalar.  The
%   units of 2*pi/OMEGA should be the same as the units of DT.
%
%   PERIODINDEX constructs an index into a vector of the same length as 
%   OMEGA that skips every N times the instantaneous period 2*pi/OMEGA.
%   Thus N=1 returns an index with one sample per period, etc.  The first N
%   periods, and last N periods, are both omitted from the index. 
%
%   If NaNs are found within OMEGA, this is interpreted as being separate
%   ridges output by RIDGEWALK.  Then PERIODINDEX applies itself to each 
%   cell separately and returns the result in one long array, with no NaNs.
%
%   INDEX=PERIODINDEX(OMEGA,N) also works, with DT defaulting to unity.  In
%   this case OMEGA must have units of radians per sample interval.
%
%   PERIODINDEX is useful with ELLIPSEPLOT for plotting ellipses a
%   specified number of periods apart.
%
%   See also ELLIPSEPLOT.
%
%   Usage: index=periodindex(dt,omega,N);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2018 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--f')
   periodindex_fig,return
end

if nargin==3
    dt=varargin{1};
    varargin=varargin(2:3);
else
    dt=1;
end

om=varargin{1};
N=varargin{2};

if iscell(om)
    index=cell(length(om),1);
    for i=1:length(om)
        index{i,1}=periodindex(dt,om{i},N);
    end
else
    if ~any(isnan(om))
        index=periodindex_loop(dt,abs(om),N);
    else
        [num,a,b]=blocknum(cumsum(isnan(om))); 
        a=a+1;a(1)=1;  %a is now heads of chunks
        om=col2cell(om);
        index=cell(length(om),1);%length(om)
        for i=1:length(om)
            index{i}=a(i)-1+periodindex_loop(dt,abs(om{i}),N);
        end
        index=cell2col(index);
        index=index(isfinite(index));
    end
end
    

function[index]=periodindex_loop(dt,om,N)

skip=ceil(N*frac(2*pi,om.*dt));
%minmin(skip)

index=[];
if skip(1)<length(om)
    index=skip(1);
    %index(end)+skip(index(end)),length(om)
    while index(end)+skip(index(end))<length(om)   
     %   index(end),skip(index(end))
        index(end+1,:)=index(end)+skip(index(end));
    end
end

function[]=periodindex_fig

load npg2006
use npg2006

%Decide on frequencies
fs=2*pi./(logspace(log10(10),log10(100),50)');

%Compute wavelet transforms using generalized Morse wavelets
[wx,wy]=wavetrans(real(cx),imag(cx),{1,2,4,fs,'bandpass'},'mirror');
[wp,wn]=vectmult(tmat,wx,wy);

%Form ridges of component time series
[ir,jr,wr,fr]=ridgewalk(dt,wn,fs,1.5,'phase'); 


index=periodindex(dt,fr,1);
figure
plot(2*pi./fr),hold on
plot([1;index(1:end-1)],diff([1;index])*dt,'o')
title('Difference between index points from PERIODINDEX for NPG dataset')

