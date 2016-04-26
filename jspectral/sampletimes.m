function[meandt,sigdt,meddt,maxdt]=sampletimes(varargin)
%SAMPLETIMES  Computes mean sampling intervals and their statistics.
%
%   DT=SAMPLETIMES(NUM) where NUM is a column vector of times, returns the 
%   mean interval between sample times DT.
%
%   [DT,SIGDT,MEDDT,MAXDT]=SAMPLETIMES(NUM) also returns the standard 
%   deviation of sample intervals SIGDT, the median sample interval MEDDT,
%   and the maximum sample interval MAXDT.
%
%   IF NUM is a cell array of numerical arrays, the output of SAMPLETIMES
%   will be arrays of size LENGTH(NUM) x 1.
%   __________________________________________________________________
%
%   Parallelization
%
%   SAMPLETIMES(NUM,'parallel') parallelizes the computation using a PARFOR 
%   loop.  This requires that Matlab's Parallel Computing Toolbox be 
%   installed, and is useful for very large datasets.
%   __________________________________________________________________
%
%   Usage: dt=sampletimes(num); 
%          [dt,sigdt,meddt,maxdt]=sampletimes(num);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2015 J.M. Lilly --- type 'help jlab_license' for details
 

num=varargin{1};
cores='serial';
if ischar(varargin{end})
    if strcmpi(varargin{end}(1:3),'par')||strcmpi(varargin{end}(1:3),'ser')
        cores=varargin{end};
    end
    varargin=varargin(1:end-1);
end
if strcmpi(cores(1:3),'par')
    if exist('parpool')~=2
        disp('Sorry, parallel algorithm requires the Parallel Computing Toolbox.')
        disp('Defaulting to the standard algorithm.')
        cores='serial';
    end
end

nout=nargout;
if ~iscell(num)
    [meandt,sigdt,meddt,maxdt]=sampletimes_one(num);
else
    [meandt,sigdt,meddt,maxdt]=vzeros(length(num),1);
    if strcmpi(cores(1:3),'par')
        parfor i=1:length(num)
            if nout>3
                [meandt(i),sigdt(i),meddt(i),maxdt(i)]=sampletimes_one(num{i});
            elseif nout>2
                [meandt(i),sigdt(i),meddt(i)]=sampletimes_one(num{i});
            elseif nout>1
                [meandt(i),sigdt(i)]=sampletimes_one(num{i});
            else
                meandt(i)=sampletimes_one(num{i});
            end
        end
    else
        for i=1:length(num)
            if nout>3
                [meandt(i),sigdt(i),meddt(i),maxdt(i)]=sampletimes_one(num{i});
            elseif nout>2
                [meandt(i),sigdt(i),meddt(i)]=sampletimes_one(num{i});
            elseif nout>1
                [meandt(i),sigdt(i)]=sampletimes_one(num{i});
            else
                meandt(i)=sampletimes_one(num{i});
            end
        end
    end
end

function[meandt,sigdt,meddt,maxdt]=sampletimes_one(num)

if isempty(num)||length(num)==1
   meandt=nan;
   sigdt=nan;
   meddt=nan;
   maxdt=nan;
else
    dt=num(2:end,:)-num(1:end-1,:);
    
    meandt=vmean(dt,1);
    if nargout>1
        sigdt=vstd(dt,1);
    end
    if nargout >2
        meddt=vmedian(dt,1);
    end
    if nargout >3
        maxdt=maxmax(dt);
    end
end
%function[]=sampletime_test
 
%reporttest('SAMPLETIME',aresame())
