function[mult,ztot]=ridgemult(varargin)
%RIDGEMULT  Ridge multiplicity, the number of simultaneous ridges present.
%
%   MULT=RIDGEMULT(IR) where IR is a column array index into wavelet ridge
%   temporal locations, as output by RIDGEWALK, returns the number of 
%   ridges present at each time, MULT, which is the same size as IR.
%
%   Following the convention of RIDGEWALK, IR may contain multiple ridges,
%   which are separated by a NaN at the end of each ridge.
%
%   [MULT,ZTOT]=RIDGEMULT(IR,ZHAT), where ZHAT is the estimated signal
%   associated with the ridges, returns the total signal ZTOT. All input
%   and output arguments are the same size. 
%   _______________________________________________________________________
%
%   Cell array input
%
%   Often one has an array of K different time series, which are grouped 
%   into a cell array Z with one times series per cell.  In this case,
%   calling EDDYRIDGES leads to output fields with one ridge per cell.  
%
%   In the output of EDDYRIDGES, IR is now a cell array of indices into
%   time locations, with one ridge per cell, KR denotes  a cell array of
%   indices into the K different time series, and ZHAT is a cell array of
%   the corresponding estimated signals.
%
%   Note that KR will take on one value along the entirety of any ridge.
%
%   MULT=RIDGERECON(IR,KR) and [MULT,ZHAT]=RIDGERECON(IR,KR,ZHAT) also 
%   work in this case, where not MULT and ZHAT are also cell arrays.
%
%   Usage: mult=ridgemult(ir);
%          [mult,ztot]=ridgemult(ir,zhat);
%          mult=ridgemult(ir,kr);
%          [mult,ztot]=ridgemult(ir,kr,zhat);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2019--2021 J.M. Lilly --- type 'help jlab_license' for details
 
ir=varargin{1};
bcell=iscell(ir);

if bcell
    ir=varargin{1};
    kr=varargin{2};
    zhat=varargin{3};
else
    kr=1+0*ir;
    zhat=varargin{2};
end

if ~bcell
    tempir=ir;
    tempkr=kr;
    tempzhat=zhat;
    clear ir kr zhat
    ir{1}=tempir;
    kr{1}=tempkr;
    zhat{1}=tempzhat;
end    

ira=cellmin(ir);
irb=cellmax(ir);
ks=cellfirst(kr);
mult=ir;

cell2col(ir,kr,zhat);
if nargout==1
    for i=1:length(ks)
        bool=(kr==ks(i));
        mult{i}=ridgemult_one(ira(i),irb(i),ir(bool));
    end
else
    ztot=mult;
    for i=1:length(ks)
        bool=(kr==ks(i));
        %vsize(ir,kr,zhat,find(bool))
        [mult{i},ztot{i}]=ridgemult_one(ira(i),irb(i),ir(bool),zhat(bool));
    end
end


function[mult,ztot]=ridgemult_one(ira,irb,ir,zhat)

ybin=ira-1/2:1:irb+1/2;
xbin=[-1/2 1/2]';
if nargin==1
    mult=twodhist(0*ir,ir,xbin,ybin);
else
    [zbar,~,~,mult]=twodstats(0*ir,ir,zhat,xbin,ybin);
    ztot=zbar.*mult;
end

