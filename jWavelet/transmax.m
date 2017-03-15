function[varargout]=transmax(varargin)
% TRANSMAX  Locates modulus maximum points of wavelet transform.
%
%   This function is part of 'element analysis' described in Lilly (2017), 
%   "Element analysis: a wavelet-based method for analyzing time-localized
%   events in noisy time series", submitted.  Available at www.jmlilly.net.
%  
%   [INDEX,WW,FF]=TRANSMAX(FA,W) gives the indicies INDEX, transform values
%   WW, and scale frequencies FF of all the maxima points in time and scale
%   of the wavelet transform W, which has been performed at frequencies FS.
%
%   W is oriented with time in rows and different scales or frequencies in
%   its columns, as output by WAVETRANS.  FS is a length SIZE(W,2) array
%   giving the scale frequencies at which the wavelet transform was taken.
% 
%   The output variables are all column arrays.
%
%   TRANSMAX applies a central difference to both rows and columns to 
%   locate maxima, then applies a quadratic interpolation in columns to 
%   obtain a refined estimate of the scale locations.  The output fields FF 
%   and WW represent interpolated values at these refined scale locations.
%   
%   Note that W can have up to four dimensions.  This is useful for
%   datasets containing many replicate time series. Use 
% 
%           [II,JJ,KK,RR]=IND2SUB(SIZE(W),INDEX)
%
%   to convert INDEX into row numbers II, column numbers JJ, and so forth.
%
%   TRANSMAX(...,CHI) returns only those maxima points having magnitudes
%   exceeding CHI, i.e. ABS(W)>CHI.  CHI may be a scalar, or an array of
%   length SIZE(W,2).  In the latter case, different cutoffs are applied 
%   to maximum points located in the different columns of W.
%   _______________________________________________________________________
%
%   Missing data fraction
%
%   TRANSMAX can also return a measure of how much each detected transform
%   maximum overlaps periods of interpolated or missing data. 
%
%   It is sometimes the case that one is interested in transform maxima for
%   time series that have missing data points or gaps.  In this case, one
%   can just interpolate over all the missing data points, pass the
%   interpolated data to TRANSMAX, and then find the missing data fraction.
%
%   [INDEX,WW,FF,RR]=TRANSMAX(FS,W,{GAMMA,BETA,BOOL}), where the third 
%   input argument is a cell array, returns the missing data fraction RR.  
%   Here BOOL is a boolean variable of the same length as the analyzed time
%   series, which is true if a datapoint was originally missing.  
% 
%   If GAMMA and BETA are the parameters of the Morse wavelet used in the
%   transform, the wavelet time-domain width at scale frequency FS is
%   approximately L=2*SQRT(2)*SQRT(BETA*GAMMA)./FS, see Lilly (2017).
%
%   The missing data fraction RR is defined as the ratio of the number of
%   missing points in a window of half-width ROUND(L/2) centered on each
%   maximum, to the total number of points in the window.  Points in the 
%   window but outside the edges of the time series are considered missing.
%
%   The missing data fraction can then be used as a constraint to exclude,
%   for example, transform maxima with RR > 0.1 or 10% missing data.
%   _______________________________________________________________________
%
%   See also MAXPROPS, TRANSMAXDIST, ISOMAX, MAX2EDDY.
%
%   Usage: [index,ww,ff]=transmax(fs,w);
%          [index,ww,ff]=transmax(fs,w,chi);
%          [index,ww,ff,rr]=transmax(fs,w,{ga,be,bool});
%          [index,ww,ff,rr]=transmax(fs,w,{ga,be,bool},chi);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2017 J.M. Lilly --- type 'help jlab_license' for details        

if strcmpi(varargin{1},'--f')
    transmax_fig;return
end

fs=varargin{1};
w=varargin{2};
varargin=varargin(3:end);
Lcell=[];
chi=0;

for i=1:2
    if length(varargin)~=0
        if iscell(varargin{1})
            Lcell=varargin{1};
        else
            chi=varargin{1};
        end
        varargin=varargin(2:end);
    end
end

[index,ww]=transmax1(w,chi);
[ww,ff,jj]=transmax_interp1(fs,w,index);
    


varargout{1}=index;
varargout{2}=ww;
varargout{3}=ff;

rr=[];
if ~isempty(Lcell)&&~isempty(index)
    [ii,jj,nn,kk]=ind2sub(size(w),index);
    ga=Lcell{1};
    be=Lcell{2};
    bool=Lcell{3};
    %cutoff=Lcell{5};
    L=2*sqrt(2)*sqrt(ga*be)./ff;
    rr=zeros(size(ii));
    for i=1:length(ii)
        a=max(1,ii(i)-round(L(i)/2));
        b=min(size(bool,1),ii(i)+round(L(i)/2));
        %This is the missing fraction, not quite what I wanted
        %rr(i)=frac(length(find(bool(a:b,nn(i),kk(i)))),b-a+1);
        %This is the fraction that is either missing or overlapping edges
        rr(i)=1-frac(length(find(~bool(a:b,nn(i),kk(i)))),2*round(L(i)/2)+1);
    end
end
varargout{4}=rr;


function[ww,ff,jj]=transmax_interp1(fs,w,index)

if ~isempty(index)
    
    [ii,jj,mm,kk]=ind2sub(size(w),index);
    indexo=sub2ind(size(w),ii,jj,mm,kk);
    indexp=sub2ind(size(w),ii,jj-1,mm,kk);
    indexn=sub2ind(size(w),ii,jj+1,mm,kk);
    
    [temp,jjhat]=quadinterp(jj-1,jj,jj+1,abs(w(indexp)),abs(w(indexo)),abs(w(indexn)));
    ww=quadinterp(jj-1,jj,jj+1,w(indexp),w(indexo),w(indexn),jjhat);
    ff=interp1(1:length(fs),fs,jjhat);  %Scale frequency inside transfrom
    jj=jjhat;
else
    ww=[];
    ff=[];
    jj=[];
end

function[index,xx]=transmax1(x,chi)
boolo=~(x==0);
x=x+randn(size(x))*eps;%Add noise to prevent duplicate max for very smooth transforms
%Check magnitude of X versus its neighbors  

xshift=vshift(x,1,1);  bool=(abs(x)>abs(xshift));
xshift=vshift(x,-1,1); bool=bool&(abs(x)>abs(xshift));
xshift=vshift(x,1,2);  bool=bool&(abs(x)>abs(xshift));
xshift=vshift(x,-1,2); bool=bool&(abs(x)>abs(xshift));

%xshift=vshift(x,1,1);  bool=bool.*(abs(x)>=abs(xshift));
%xshift=vshift(x,-1,1); bool=bool.*(abs(x)>=abs(xshift));
%xshift=vshift(x,1,2);  bool=bool.*(abs(x)>=abs(xshift));
%xshift=vshift(x,-1,2); bool=bool.*(abs(x)>=abs(xshift));

% xshift=vshift(x,1,1); bool=bool.*(x>=xshift);%1
% xshift=vshift(x,1,2); bool=bool.*(x>=xshift);%2
% xshift=vshift(x,-1,1);bool=bool.*(x>=xshift);%3
% xshift=vshift(x,-1,1);bool=bool.*(x>=xshift);%4
% xshift=vshift(x,-1,2);bool=bool.*(x>=xshift);%5
% xshift=vshift(x,-1,2);bool=bool.*(x>=xshift);%6
% xshift=vshift(x,1,1); bool=bool.*(x>=xshift);%7
% xshift=vshift(x,1,1); bool=bool.*(x>=xshift);%8
% 
%6 5 4
%7 x 3
%8 1 2

%Exclude points on edges
bool([1 end],:,:,:)=false;
bool(:,[1 end],:,:)=false;
bool=bool&boolo;

%[ii,jj]=find(bool);
%figure,plot(ii,jj,'.')

%Find value of X at isolated maxima
index=find(bool);
if ~isempty(index)
    xx=x(index);
    
    %Allow for CHI being of the same size as the # of columns 
    if length(chi(:))==size(x,2)
        [ii,jj]=ind2sub(size(x),index);
        chi=chi(jj);
    end
    
    %Remove points below cutoff
    vindex(xx,index,abs(xx)>chi,1);
    
    %%Sort remainer
    %if ~isempty(xx)
    %    [temp,sorter]=sort(-xx);
    %    vindex(xx,index,sorter,1);     
    %end
else
    xx=[];
end

%if isempty(ii)
%    disp('TRANSMAX found no maxima points.')
%aend
