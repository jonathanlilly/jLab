function[lr]=ridgelen(varargin)
%RIDGELEN  Wavelet ridge length expressed as number of full cycles.
%
%   LEN=RIDGELEN(FR) determines the length of the ridges given by the
%   ridge frequency FR, with LEN expressed as number of cycles completed
%   along the ridge.  FR has units of radians per sample.
%
%   LEN is a column vector of the same size as FR.  Each element of LEN 
%   gives the number of cycles contained in the current ridge. 
%
%   LEN=RIDGELEN(DT,FR) uses the sample interval DT in calculating the 
%   ridge length.  DT has the default value DT=1.  In this case, FR is 
%   expected to have units of radians per DT.
%
%   'ridgelen --t' runs some tests.
%
%   Usage: lr=ridgelen(fr);
%          lr=ridgelen(dt,fr);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2018 J.M. Lilly --- type 'help jlab_license' for details

%   RIDGELEN is a low-level function called by RIDGEWALK.

%   RIDGELEN internally computes an ID that takes on a single value for all
%   points in the same ridge.  If such an ID is already available, calling
%   LEN=RIDGELEN(DT,ID,IR,FR) also works.  This form is used by RIDGEWALK.  

%   This is no longer supported, but if you want it back uncomment the
%   commented-out code below.
%  
%   FR can also be cell arrays of ridges, in which case the output LEN is 
%   also a cell array of the same size as the input variables.
%
%   If FR is cell array, DT may either be a scalar or an array of length 
%   LENGTH(IR).

if strcmpi(varargin{1},'--t')
    ridgelen_test;return
end

lr=[];
%From within RIDGEWALK, ID is already known and the call is lr=ridgelen(1,id,fr);
dt=1;
id=[];
if nargin==3
    dt=varargin{1};
    id=varargin{2};
    fr=varargin{3};
elseif nargin==2
    dt=varargin{1};
    fr=varargin{2};
elseif nargin==1
    fr=varargin{1};
end

if ~isempty(fr)
    lr=ridgelenloop(dt,id,fr);
end

% if ~isempty(fr)
%     if iscell(fr)
%         for k=1:length(fr)
%             if isempty(id)
%                 idk=[];
%             else
%                 idk=id{k};
%             end
%             if length(dt)==1
%                 dtk=dt;
%             else
%                 dtk=dt(k);
%             end
%             lr{k,1}=ridgelenloop(dtk,idk,fr{k});
%         end
%     else
%         lr=ridgelenloop(dt,id,fr);
%     end
% end

function[lr]=ridgelenloop(dt,id,fr) 
if isempty(id)
    id=cumsum(isnan(vshift(fr,-1,1)),1);
end

index=~isnan(fr);
lr=nan*fr;
if ~isempty(index)
    lr(index)=ridgelen1(dt,id(index),fr(index));
end

    
function[lr]=ridgelen1(dt,id,fr) 

%vsize(dt,id,ir,fr)
fr=fr./(2*pi).*dt;  %Convert to cyclic frequency per sample interval
 
[num,a,b]=blocknum(id);
%deal with start and end nans
%npoints=~isnan(fr);
vswap(fr,nan,0);
%angr=cumsum(fr.*dt,1);
ar=cumsum(fr,1); %Ridge age
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
lr=ridgelen(fr);
reporttest('RIDGELEN', aresame(lr,2+0*ir,1e-8))

ir=[1:101]';
fr=frac(2*pi,50)/10+0*ir;
lr2=ridgelen(10,fr);
reporttest('RIDGELEN sample interval invariance', aresame(lr,lr2,1e-8))




