function[varargout]=griddrifters(varargin)
%GRIDDRIFTERS  Average drifter velocities onto a space/time 3D grid.
%
%   STRUCT=GRIDDRIFTERS(NUM,LAT,LON,U,V,FILLED,SOURCE,LATO,LONO,YO,MO,N)
%   returns a structure containing velocities averaged onto a 3D grid.  
%
%   The grid has bin edges given by LATO and LONO, and has month-long bins
%   in time with a semimonthly spacing beginning with year YO and month MO,
%   and having N total time steps.  
%
%   The grid size is thus LENGTH(LATO)-1 x LENGTH(LONO)-1 x N.  The first
%   time bin corresponds to year YO and month MO, i.e. it is centered on
%   the midpoint of this month.  Odd slices correspond to month bins, e.g.
%   Jan 1992, while even slices run between midpoints of adjacent months.
%
%   NUM,LAT,LON,U,V,FILLED,SOURCE are all cell arrays of Lagrangian drifter
%   or float data, with one instrument per cell, such as those used by the
%   FLOATS.MAT and DRIFTERS.MAT dataset.  These are all the same size. 
%
%   NUM is the date in Matlab's DATENUM format, LAT and LON are latitude
%   and longitude, U and V are eastward and northward velocity components, 
%   FILLED is a flag that is true if the data point is filled, and SOURCE
%   contains integers indicating a source dataset.
%
%   FILLED and SOURCE are both optional.  To omit one or both of these,
%   replace with the empty array, [], in the input list. 
%
%   If FILLED is included, only non-filled data is used for the averages.
%
%   The output structure STRUCT has the following fields:
%
%       STRUCT.NUM       Date of bin midpoint in DATENUM format
%       STRUCT.LAT       Latitudes of bin centers
%       STRUCT.LON       Longitudes of bin centers
%       STRUCT.U         Bin-averaged eastward velocity
%       STRUCT.V         Bin-averaged northward velocity
%       STRUCT.EPSUU     Instantaneous eastward local variance in each bin 
%       STRUCT.EPSVV     Instantaneous northward local variance in each bin
%       STRUCT.EPSUV     Instantaneous local covariance in each bin
%       STRUCT.COUNT     Number of data points averaged in each bin
%
%   These are all 3D arrays of the same size.  Array entries in which there
%   is no data are filled with NaNs, apart from COUNT which will be zero.
%
%   The "instantaneous local variance" is the variance relative to the 
%   local mean velocity computed withing each 3D bin.  See Lilly and Perez-
%   Brunius (2021), "A gridded surface current product for the Gulf of
%   Mexico from consolidated drifter measurements", for details.
%
%   If SOURCE is input, then STRUCT also includes another field:
%
%       STRUCT.COUNTS    Number of data points from each source in each bin
%
%   which is a 4D array of LENGTH(LATO)-1 x LENGTH(LONO)-1 x N x M, where M
%   is the maximum value occurring in SOURCE.
%
%   Usage: struct=griddrifters(num,lat,lon,u,v,[],[],lato,lono,yo,mo,N);
%          struct=griddrifters(num,lat,lon,u,v,filled,source,...
%                                                     lato,lono,yo,mo,N);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2020 J.M. Lilly --- type 'help jlab_license' for details
 
%if strcmp(varargin{1}, '--t')
%    griddrifters_test,return
%end
 
num=varargin{1};
lat=varargin{2};
lon=varargin{3};
u=varargin{4};
v=varargin{5};
filled=varargin{6};
source1=varargin{7};
lato=varargin{8};
lono=varargin{9};
yo=varargin{10};
mo=varargin{11};
N=varargin{12};

cv=cellpair(u,v);

if ~isempty(source1)
    sourceinput=true;
    if ~iscell(source1)
        %num=col2cell(num);
        clear source
        for i=1:length(num)
            source{i,1}=source1(i)+zeros(size(num{i}));
        end
        %cell2col(num,source);
        %cell2col(num,lon,lat,cv,filled,depth);
        %vsize(num,lat,lon,filled,depth,cv,source)
    else
        source=source1;
    end
else
    source=num;
    sourceinput=false;
end
% %acceleration
% cvcell=col2cell(cv);
% cva=cvcell;
% for i=1:length(cvcell)
%     cva{i}=vdiff(3600,cvcell{i},1);
% end
% cell2col(cva);
% clear cvcell
if ~isempty(filled)
    %vsize(num,lon,lat,cv,filled,source)
    %iscell(source),iscell(num),iscell(lat)
    cell2col(num,lon,lat,cv,filled,source);
    vindex(num,lon,lat,cv,source,filled==0&isfinite(cv),1);
else
    cell2col(num,lon,lat,cv,source);
end 
    
    %vindex(num,lon,lat,cv,cva,source,filled==0&isfinite(cv),1);

%number of months
[yf,mf]=yearfrac(num);
monum=(floor(yf)-yo)*12+mf-mo+1;%month number . fraction from start date
[matmat,u,v,epsuu,epsvv,epsuv]=vzeros(length(lato)-1,length(lono)-1,N,nan);
[nums,mnums]=vzeros(N,1,nan);

count=vzeros(length(lato)-1,length(lono)-1,N);
counts=vzeros(length(lato)-1,length(lono)-1,N,maxmax(source));

for i=1:N
    %nums(i)=yf2num(yo+mo+i/24+(mo-1)/12); %not quite right
    if isodd(i)
        bool=(floor(monum)==(i-1)/2+1);%[1,3,5,7] => [1,2,3,4]
        nums(i)=(1/2)*(datenum(yo,mo+(i-1)/2,1)+datenum(yo,mo+1+(i-1)/2,1)-1);
    elseif iseven(i)  %i/2  %[2,4,6] => [1,2,3]
        bool=(monum>=i/2+1/2)&(monum<i/2+3/2);  %half-months [>1.5<2.5, etc.]
        nums(i)=datenum(yo,mo+i/2,1);
    end
    %datestr(nums)
          
      % [ length(find(bool)),length(find(isfinite(abs(cv(bool)))))]
   
    if ~isempty(find(bool))
        [matmat(:,:,i),xmid,ymid,count(:,:,i)]=twodstats(lon(bool),lat(bool),abs(cv(bool)),lono,lato);
        [mz,~,~,numz,cov]=twodstats(lon(bool),lat(bool),[real(cv(bool)),imag(cv(bool))],lono,lato);
      %  figure,jpcolor(bool)
        u(:,:,i)=mz(:,:,1);v(:,:,i)=mz(:,:,2);
        epsuu(:,:,i)=cov(:,:,1,1);
        epsvv(:,:,i)=cov(:,:,2,2);
        epsuv(:,:,i)=cov(:,:,1,2);
        %mz=twodstats(lon(bool),lat(bool),[real(cva(bool)),imag(cva(bool))],lono,lato);
        %ua(:,:,i)=mz(:,:,1);va(:,:,i)=mz(:,:,2);
        mnums(i)=vmean(num(bool),1);
        count(:,:,i)=numz;
        if sourceinput
            for m=1:maxmax(source)
                boolm=bool&source==m;
                counts(:,:,i,m)=twodhist(lon(boolm),lat(boolm),lono,lato);
            end     
        end
    end
   % figure,jpcolor(u(:,:,i))
%  end
end

%datestr(nums(1))

% [mat,mu,mv]=vmean(matmat,u,v,3);
% ubar=interp2(xmid,ymid,mu,lon,lat,'nearest');
% vbar=interp2(xmid,ymid,mv,lon,lat,'nearest');
% 
% for i=1:N
%     if isodd(i)
%         bool=(floor(monum)==(i-1)/2+1);%[1,3,5,7] => [1,2,3,4]
%     elseif iseven(i)  %i/2  %[2,4,6] => [1,2,3]
%         bool=(monum>=i/2+1/2)&(monum<i/2+3/2);  %half-months
%     end
%     if ~isempty(find(bool))
%         siguu(:,:,i)=twodstats(lon(bool),lat(bool),...
%             squared(real(cv(bool))-ubar(bool)),lono,lato);
%         sigvv(:,:,i)=twodstats(lon(bool),lat(bool),...
%             squared(imag(cv(bool))-vbar(bool)),lono,lato);
%         siguv(:,:,i)=twodstats(lon(bool),lat(bool),...
%             (real(cv(bool))-ubar(bool)).*(imag(cv(bool))-vbar(bool)),lono,lato);
%     end
% end

%[msiguu,msigvv,msiguv]=vmean(siguu,sigvv,siguv,3);
%[mepsuu,mepsvv,mepsuv]=vmean(epsuu,epsvv,epsuv,3);
    
struct.num=nums;
struct.lat=ymid;
struct.lon=xmid;
struct.u=u;
struct.v=v;
%struct.siguu=siguu;
%struct.sigvv=sigvv;
%struct.siguv=siguv;
struct.epsuu=epsuu;
struct.epsvv=epsvv;
struct.epsuv=epsuv;
struct.count=count;
if sourceinput 
    struct.counts=counts;
end

varargout{1}=struct;

    
%function[]=griddrifters_test
 
%reporttest('GRIDDRIFTERS',aresame())
