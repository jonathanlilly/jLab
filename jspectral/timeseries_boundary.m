function[y]=timeseries_boundary(varargin)
%TIMESERIES_BOUNDARY  Apply boundary conditions to data before transform.
%
%   TIMESERIES_BOUNDARY is a low-level function called by WAVETRANS and
%   ANATRANS.
%
%   TIMESERIES_BOUNDARY applies periodic, zero-padded, or mirror boundary
%   conditions to a time series.  See ANATRANS or WAVETRANS.
%
%   This is a low-level function, not meant to be user called by users.
%
%   TIMESERIES_BOUNDARY(X,STR,DETRENDSTR) applies boundary conditions 
%   specified by STR to time series X.  If DETRENDSTR='detrend', X is first
%   detrended.  STR may be 'periodic', 'zeros', 'mirror', or 'reverse':
%
%      'periodic' wraps the beginning around to the end, and vice versa
%      'zeros'    extends the time series with zeros at both edges
%      'mirror'   reflects the time series about the beginning and the end
%      'reverse'  is like mirror, but with a change in sign, and a shift
%                   to make the time series be continuous at both edges
%
%   TIMESERIES_BOUNDARY(X,DIM,STR) applies the boundary conditions
%   specified by STR along dimension DIM.  By default, DIM=1.  Note that
%   no detrending is possible for DIM not equal to one.  Also, 'reverse'
%   only works for the default case of DIM=1.  
%
%   Usage: y=timeseries_boundary(x,str,detrendstr);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2015 J.M. Lilly --- type 'help jlab_license' for details
 

if strcmpi(varargin{1}, '--t')
    timeseries_boundary_test,return
end


%x,str,bdetrend
x=varargin{1};

detrendstr='nodetrend';
str='periodic';

for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'nod')||strcmpi(varargin{end}(1:3),'det')
            detrendstr=varargin{end};
        else
            str=varargin{end};
        end
        varargin=varargin(1:end-1);
    end
end

if length(varargin)==2
    dim=varargin{2};
else
    dim=1;
end

%Prepare data by applying boundary condition

if dim==1
    for i=1:size(x,2)
        ai=find(~isnan(x(:,i)),1,'first');
        bi=find(~isnan(x(:,i)),1,'last');
        if isempty(ai) && isempty(bi)
            disp(['Warning: Data column ', int2str(i), ' contains no finite values.'])
            a(i)=1;
            b(i)=size(x,1);
        elseif any(~isfinite(x(ai:bi,i)))
            error(['Data contains interior NANs or INFs in column ', int2str(i), '.'])
        else
            a(i)=ai;
            b(i)=bi;
        end
    end
end
sizex=size(x);
M=sizex(dim);
if ~strcmpi(str,'periodic')
   sizex(dim)=sizex(dim)*3;
end
y=zeros(sizex);

if dim==1
    for i=1:size(x,2)
        index{i}=a(i):b(i);
        indexy{i}=(M+a(i)-length(index{i}):M+2*length(index{i})+a(i)-1);
        xi=x(index{i},i);
        
        if strcmpi(detrendstr(1:3),'det')
            xi=detrend(xi);
        end
        
        if strcmpi(str,'zeros')
            y(indexy{i},i)=[0*xi;xi;0*xi];
        elseif strcmpi(str,'mirror')
            y(indexy{i},i)=[flipud(xi);xi;flipud(xi)];
        elseif strcmpi(str,'reverse')
            xia=-flipud(xi);
            xia=xia-xia(end)+xi(1)-(xi(2)-xi(1));
            xib=-flipud(xi);
            xib=xib-xib(1)+xi(end)-(xi(end-1)-xi(end));
            y(indexy{i},i)=[xia;xi;xib];
        elseif strcmpi(str,'periodic')
            y(index{i},i)=xi;
        else
            error(['Transform option STR = ''',str,''' is not supported.']);
        end
    end
else
       if strcmpi(str,'zeros')
            y=vindexinto(y,x,M+1:2*M,dim);
        elseif strcmpi(str,'mirror')
            y=vindexinto(y,flip(x,dim),1:M,dim);
            y=vindexinto(y,x,M+1:2*M,dim);
            y=vindexinto(y,flip(x,dim),2*M+1:3*M,dim);
        elseif strcmpi(str,'periodic')
            y=x;
        else
            error(['Transform option STR = ''',str,''' is not supported.']);
       end
end
y=vswap(y,nan,0);


function[]=timeseries_boundary_test

load solomon 
use solomon

x=[x y z];
zp=timeseries_boundary(x,'periodic');
zm=timeseries_boundary(x,'mirror');
zz=timeseries_boundary(x,'zeros');

zp2=timeseries_boundary(permute(x,[3 2 1]),3,'periodic');
zm2=timeseries_boundary(permute(x,[3 2 1]),3,'mirror');
zz2=timeseries_boundary(permute(x,[3 2 1]),3,'zeros');

bool(1)=aresame(zp,permute(zp2,[3 2 1]));
bool(2)=aresame(zm,permute(zm2,[3 2 1]));
bool(3)=aresame(zz,permute(zz2,[3 2 1]));

reporttest('TIMESERIES_BOUNDARY application to alternate dimension',allall(bool))



