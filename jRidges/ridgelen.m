function[lr]=ridgelen(varargin)
%RIDGELEN  Wavelet ridge length expressed as number of full cycles.
%
%   LEN=RIDGELEN(IR,FR) determines the length of the ridges given by the
%   ridge parameters IR and FR, with LEN expressed as number of cycles
%   completed along the ridge.  FR has units of radians per sample.
%
%   IR is an array of the ridge time indices, while FR is the frequency of 
%   the wavelet transform value along the ridge, as output by RIDGEWALK. 
%
%   LEN is a column vector of the same size as IR and FR.  Each element of
%   LEN gives the number of cycles contained in the current ridge. 
%
%   IR and FR can also be cell arrays of ridges, in which case the output
%   LEN is also a cell array of the same size as the input variables.
%   _____________________________________________________________________
%   
%   Sample interval
%
%   LEN=RIDGELEN(DT,IR,FR) uses the sample interval DT in calculating the 
%   ridge length.  DT has the default value DT=1.
%
%   In this case, FR is expected to have units of radians per DT.
%
%   If IR and FR are cell arrays, DT may either be a scalar or an array of
%   length LENGTH(IR).
%   _______________________________________________________________________
%
%   RIDGELEN internally computes an ID that takes on a single value for all
%   points in the same ridge.  If such an ID is already available, calling
%   LEN=RIDGELEN(DT,ID,IR,FR) also works.  This form is used by RIDGEWALK.  
%   _______________________________________________________________________
%
%   'ridgelen --t' runs some tests.
%
%   Usage: lr=ridgelen(ir,fr);
%          lr=ridgelen(dt,ir,fr);
%          lr=ridgelen(dt,id,ir,fr);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2014 J.M. Lilly --- type 'help jlab_license' for details

%   RIDGELEN is a low-level function called by RIDGEWALK.

if strcmpi(varargin{1},'--t')
    ridgelen_test;return
end

args=varargin;
lr=[];
%In call from within RIDGEWALK, ID is already known and the call is 
%lr=ridgelen(1,id,ir,fr);
dt=1;
id=[];
if nargin==4
    dt=args{1};
    id=args{2};
    args=args(3:end);
elseif nargin==3
    dt=args{1};
    args=args(2:end);
end

ir=args{1};
fr=args{2};
    
if ~isempty(ir)
    if iscell(ir)
        for k=1:length(ir)
            if isempty(id)
                idk=[];
            else
                idk=id{k};
            end
            if length(dt)==1
                dtk=dt;
            else
                dtk=dt(k);
            end
            lr{k,1}=ridgelenloop(dtk,idk,ir{k},fr{k});
        end
    else
        lr=ridgelenloop(dt,id,ir,fr);
    end
else
    lr=ir;
end

function[lr]=ridgelenloop(dt,id,ir,fr) 
if isempty(id)
    id=cumsum(isnan(ir),1);
end
index=~isnan(ir);
lr=nan*ir;
if ~isempty(index)
    lr(index)=ridgelen1(dt,id(index),ir(index),fr(index));
end

    
function[lr]=ridgelen1(dt,id,ir,fr) 

%vsize(dt,id,ir,fr)
fr=fr./(2*pi).*dt;  %Convert to cyclic frequency per sample interval
 
[num,a,b]=blocknum(id);
%deal with start and end nans
%npoints=~isnan(fr);
vswap(fr,nan,0);
%angr=cumsum(fr.*dt,1);
ar=cumsum(fr,1);  %Ridge age
%figure,plot(fr)
%figure,plot(ar)
lena=abs(ar(b)-ar(a));

len1=zeros(size(id));
len1(a)=lena;
len1=cumsum(len1);

len2=zeros(size(id));
len2(a)=[0;lena(1:end-1)];
len2=cumsum(len2);

lr=len1-len2;  %maximum age of each ridge

function[]=ridgelen_test

ir=[1:101]';
fr=frac(2*pi,50)+0*ir;
lr=ridgelen(ir,fr);
reporttest('RIDGELEN', aresame(lr,2+0*ir,1e-8))

ir=[1:101]';
fr=frac(2*pi,50)/10+0*ir;
lr2=ridgelen(10,ir,fr);
reporttest('RIDGELEN sample interval invariance', aresame(lr,lr2,1e-8))




