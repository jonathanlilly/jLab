function[varargout]=ridgemap(varargin)
%RIDGEMAP  Maps ridge quantities back onto the time series.
%
%   X=RIDGEMAP(M,XR,IR) where IR is a ridge index and XR is a quantity 
%   along the ridge, maps the values of XR to their correct row locations 
%   in a time series of length M, and returns the result in the array X.
%
%   If IR and XR contain L different ridges separated by NaNs, as output by
%   RIDGEWALK, then X is M x L with the values of XR from each ridge in a 
%   separate column. Values not specified by the IR are left as NaNs. 
%
%   [X1,X2,...,XM]=RIDGEMAP(M,X1R,X2R,...,XPR,IR) also works for any P
%   different ridge quantities X1R--XPR.
%
%   When using RIDGEWALK's joint ridges algorithm, in which some quantities
%   have more than one column, they should be passed to RIDGEMAP 
%   individually, for example [X1,X2]=RIDGEMAP(M,X(:,1),X(:,2),IR).
%   __________________________________________________________________
%
%   Collapsing 
%
%   X=RIDGEMAP(...'collapse') combines values from all the ridges using
%   a power-weighted mean.  Then X is a column vector of size M x 1.
%   __________________________________________________________________
% 
%   Ridge multiplicity
%
%   [...,MULT]=RIDGEMAP returns the ridge multiplicity MULT after all the
%   expected output quantities.  MULT is a column vector with size M x 1.
% 
%   The ridge multiplicity is the number of ridges present at each time.     
%   __________________________________________________________________
%
%   See also RIDGEWALK.
%
%   'ridgemap --t' runs some a test.
%
%   Usage:   x=ridgemap(M,xr,ir);
%            [x,f]=ridgemap(M,xr,fr,ir);
%            [x,f]=ridgemap(M,xr,fr,ir,'collapse');
%            [x,mult]=ridgemap(M,xr,ir);
%            [x,f,mult]=ridgemap(M,xr,fr,ir);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2019 J.M. Lilly --- type 'help jlab_license' for details

if strcmpi(varargin{1}, '--t')
    ridgemap_test,return  
end

if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
else
    str='all';
end

M=varargin{1}(1);
N=1;
bNinput=false;
if length(M)==2
    N=M(2);
    M=M(1);
    bNinput=true;
end

varargin=varargin(2:end);
if bNinput
    kr=varargin{end};
    ir=varargin{end-1};
    varargin=varargin(1:end-2);
else
    ir=varargin{end};
    kr=1+0*ir;
    varargin=varargin(1:end-1);
end

if ~isempty(ir)
    if iscell(ir)
        ir=cell2col(ir);
        for i=1:length(varargin)
            varargin{i}=cell2col(varargin{i});
        end
    end
    for i=1:length(varargin)
        [~,varargin{i}]=col2mat(ir,varargin{i});
    end
    ir=col2mat(ir);
    for i=1:length(varargin)
        varargout{i}=nan*zeros(M,size(ir,2));
%        size(varargout{i})
 %       size(ir)
        for k=1:size(ir,2)
           %ir(isfinite(ir(:,k)),k)
            varargout{i}(ir(isfinite(ir(:,k)),k),k)=varargin{i}(isfinite(ir(:,k)),k);
        end
    end
    mult=vsum(0+isfinite(varargout{1}),2);
else
    for i=1:length(varargin)
        varargout{i}=nan*zeros(M,1);
    end
    mult=zeros(M,1);
end

varargout{end+1}=mult;

% for i=1:length(varargin)
%     varargout{i}=vswap(varargout{i},inf+sqrt(-1)*inf,nan+sqrt(-1)*nan);
%     varargout{i}=vswap(varargout{i},inf,nan);
% end

%Calculate multiplicity
%mult=~isnan(vswap(varargout{1}(:,1,:),0,nan));
%The "+0" is to convert the logical into a numerical value

if strfind(str,'col')        
    for i=length(varargin):-1:1
        if i==1
            varargout{i}=vsum(varargout{i},2);
        else
            varargout{i}=vmean(varargout{i},2,squared(varargout{1}));
        end
    end
end

for i=1:length(varargin)        
    varargout{i}=squeeze(varargout{i});
    %anyany(isinf(varargout{i}))
end


function[]=ridgemap_test

load npg2006
use npg2006

%Decide on frequencies
fs=2*pi./(logspace(log10(10),log10(100),50)');

%Compute wavelet transforms using generalized Morse wavelets
[wx,wy]=wavetrans(real(cx),imag(cx),{1,2,4,fs,'bandpass'},'mirror');
[wp,wn]=vectmult(tmat,wx,wy);

%Form ridges of component time series
[wr,ir,jr,fr]=ridgewalk(dt,wn,fs); 
[wa,fa,mult]=ridgemap(length(wn),wr,fr,ir);
reporttest('RIDGEMAP has one column per ridge, non-joint ridges',size(fa,2)==size(wa,2)&&size(fa,2)==length(find(~isfinite(ir))))

%[ir,jr,wpr,wnr,fpr,fnr]=ridgewalk(dt,wp,wn,fs,{0,0});   
%[wpa,wna,fpa,fna,mult]=ridgemap(length(wn),wpr,wnr,fpr,fnr,ir);


