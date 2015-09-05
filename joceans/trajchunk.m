function[varargout]=trajchunk(varargin)
%TRAJCHUNK  Converts Lagrangian trajectories into chunks based on the Coriolis period.
%
%   TRAJCHUNK is used to split float or drifter data into chuncks such that
%   the length of each chunk is a fixed multiple of the average Coriolis
%   frequency. This is useful in spectral analysis. 
%
%   [NUMO,LATO]=TRAJCHUNK(NUM,LAT,P), where NUM and LAT are arrays of date
%   number and latitude for float or drifter data, re-organizes these into 
%   chunks such that the average Coriolis frequency f_C in each chunk is at 
%   least P times the Rayleigh frequency f_R for that chunk.
%
%   Recall that the Rayleigh frequency is f_R=2*pi/(DT*N), in units of 
%   radians per unit time, where DT is the sample interval and N is the
%   number of samples.
%
%   The output variables NUM and LAT are now cell arrays of numerical
%   arrays, with each cell is truncated such that its length is just long
%   enough such that f_C > P * f_R.  
%
%   Trajectories that are not long enough to satisfy this criterion are 
%   discarded, as are short residual segments at the end of trajectories.   
%
%   Each input trajectory is thus split into zero, one, or more than one
%   cells in the output variables.   
%   __________________________________________________________________
%
%   Additional options
%
%   [NUMO,LATO,Y1,Y2,...,YM]=TRAJCHUNK(NUM,LAT,X1,X2,...XM,N) chunks the M
%   input arrays X1, X2,... XM in the same manner, and returns these as Y1,
%   Y2,... YM.  The YM are cell arrays of the same size as NUMO.
%
%   TRAJCHUNK(...,P,LMIN) also specifies a mininum number of points LMIN
%   for each chunk.      
%
%   TRAJCHUNK with no output arguments overwrites the original named output
%   variables. 
%
%   [..,II]=TRAJCHUNK(...), with an extra final output argument, 
%   outputs a cell array II of indices to the original time series. 
%
%   As an example, LAT(II{1}) gives the latitudes of the data in the first 
%   cell of the output, LATO{1}.
%   __________________________________________________________________
%   
%   Keeping short data segments
%
%   By default, any data is cells shorter than the specified length are 
%   discarded, as are data segments at the end of the trajectories.  
%
%   TRAJCHUNK(...,'keep') keeps these instead.  Short cells are returned in
%   their own chunks, and leftover segments are appended to the end of the 
%   preceding chunk.  
%
%   This preserves the number of data points, while favoring a requested 
%   length, if possible, from each trajectory. 
%   __________________________________________________________________
%
%   Cell array input 
%
%   The input variables NUM and LAT are cell arrays of numerical arrays,
%   with one trajectory per cell, as with FLOATS.MAT and DRIFTERS.MAT.  For
%   details on these datasets, see ABOUT_FLOATS and ABOUT_DRIFTERS.
%
%   In this case, the output variables are still 1D cell arrays of numeric
%   arrays.  Chunks drawn from each successive trajectory are appended to 
%   the end of the output cell arrays.   
%
%   [...,II,KK]=TRAJCHUNK(...) in this case also outputs the indices of
%   the data locations within the input cells.  KK is not a cell array
%   like the other output arguments, but rather a row array of LENGTH(II).
%
%   As an example, LAT(KK(1))(II{1}) gives the latitudes of the data in the 
%   first cell of the output, LATO{1}.
%   __________________________________________________________________
%   
%   Overlap
%
%   TRAJCHUNK(...,'overlap',PCT) outputs chunks with a percentage PCT 
%   overlap.  For example, TRAJCHUNK(...,'overlap',50) outputs chunks that 
%   overlap by 50%.  The default behavior gives chunks with no overlap.
%   __________________________________________________________________   
%
%   See also CELLCHUNK.
%
%   'trajchunk --t' runs a test.
%
%   Usage: [num,lat]=trajchunk(num,lat,P);
%          [num,lat,lon,cv]=trajchunk(num,lat,lon,cv,P);
%          [num,lat,lon,cv]=trajchunk(num,lat,lon,cv,P,lmin);
%          [num,lat,lon,cv]=trajchunk(num,lat,lon,cv,P,'overlap',50);
%          trajchunk(num,lat,lon,cv,P);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2015 J.M. Lilly --- type 'help jlab_license' for details
 

%   Equivalently, this means that the inertial period 2*pi/f_C is just less
%   than 1/M times the chunk duration DT*N, or 2*pi/f_C < (1/M) * (DT*N). 
%   Thus M inertial oscillations fit into each chunk.
%   Not completely sure about the wording of this due to the definition of 
%   'mean Coriolis frequency' etc.

if strcmpi(varargin{1}, '--t')
    trajchunk_test,return
end

factor=1;
opt='nokeep';     %What to do with leftover points

for i=1:2
    if ischar(varargin{end})
        if strcmpi(varargin{end}(1:3),'kee')||strcmpi(varargin{end}(1:3),'nok')
            opt=varargin{end};
            varargin=varargin(1:end-1);
        end
    end
    if ischar(varargin{end-1})
        if strcmpi(varargin{end-1}(1:3),'ove')
            factor=1-varargin{end}/100;
            varargin=varargin(1:end-2);
        end
    end
end

if ~iscell(varargin{end-1})&&length(varargin{end-1})==1
    N=varargin{end-1};
    M=varargin{end};
    varargin=varargin(1:end-2);
else
    N=varargin{end};
    M=0;
    varargin=varargin(1:end-1);
end
        
num=varargin{1};
lat=varargin{2};
%Note: leave num and lat as first two entries also
na=length(varargin);

%na,N,M,size(num),size(lat),str,opt
if ~iscell(lat)
    %Create ii index 
    varargin{na+1}=[1:length(lat)]';
    index=trajchunk_index(num,lat,N,M,opt,factor);
    n=0;
    if ~isempty(index)
        for k=1:length(index)
            n=n+1;
            for j=1:na+1
                varargout{j}{n,1}=varargin{j}(index{k});
            end
        end
    end
else        
    %Create ii and kk indices
    for i=1:length(lat)
        varargin{na+1}{i}=[1:length(lat{i})]';
        varargin{na+2}{i}=i+0*lat{i};
    end
    index=[];
    for i=1:length(lat)
        if length(lat)>1000
            if res((i-1)/1000)==0
                disp(['TRAJCHUNK working on cells ' int2str(i) ' to ' int2str(min(i+1000,length(lat))) ' of ' int2str(length(lat)) '.'])
            end
        end
        index{i,1}=trajchunk_index(num{i},lat{i},N,M,opt,factor);
    end     
    n=0;
    for i=1:length(lat)
        if ~isempty(index{i})
            for k=1:length(index{i})
                n=n+1;
                for j=1:na+2
                    varargout{j}{n,1}=varargin{j}{i}(index{i}{k});
                end
            end
        end
    end
    
    %i,j,k,n
    %Convert cell array kk into numeric array
    temp=varargout{na+2};
    varargout{na+2}=zeros(1,length(varargout{na+1}));
    for i=1:length(varargout{na+1})
        varargout{na+2}(i)=temp{i}(1);
    end
end



eval(to_overwrite(na));




function[index]=trajchunk_index(num,lat,N,M,opt,factor)

dt=num(2)-num(1);

fo=abs(corfreq(lat))*24;
meanfo=frac(cumsum(fo),[1:length(fo)]');
fr=frac(2*pi,dt*[1:length(fo)]');
                %aresame(meanfo(end),vmean(fo,1))

bdone=false;
n=0;

x=[1:length(fo)]';
index=[];
while ~bdone
    ii=max(find(N*fr<meanfo,1,'first'),M);
    n=n+1;
    if isempty(ii)||ii>=length(x)
        bdone=1;
    else
        index{n}=x(1:ii);
        %N*fr(ii)<meanfo(ii)
        if ii<length(x)
            x=x(floor(ii*factor)+1:end);
            fo=fo(floor(ii*factor)+1:end);
            %fr=fr(ii+1:end)-fr(ii+1)+frac(2*pi,dt);
            meanfo=frac(cumsum(fo),[1:length(fo)]');
            fr=frac(2*pi,dt*[1:length(fo)]');
            %aresame(meanfo(end),vmean(fo,1))
        end
    end
end

% if ~isempty(index)
%     if strcmpi(opt(1:3),'kee')
%         index{end}=index{end}(1):x(end);
%     end
% end

if strcmpi(opt(1:3),'kee')
    if ~isempty(index)
        index{end}=index{end}(1):x(end);  %Append leftovers
    else
        index{1}=x;  %Keep short segments
    end
end

function[]=trajchunk_test
 

load ebasnfloats

use ebasnfloats
dt=1;
trajchunk(num,lat,lon,32);
meanfo=vmean(abs(corfreq(col2mat(cell2col(lat)))),1)'*24*dt;
fr=frac(2*pi,dt*cellength(lat));
reporttest('TRAJCHUNK no overlap',allall(meanfo>32*fr))
%length(num),length(cell2col(num))

use ebasnfloats
dt=1;
trajchunk(num,lat,lon,32,'overlap',50);
meanfo=vmean(abs(corfreq(col2mat(cell2col(lat)))),1)'*24*dt;
fr=frac(2*pi,dt*cellength(lat));
reporttest('TRAJCHUNK overlap',allall(meanfo>32*fr))
%length(num),length(cell2col(num))

% with minimum length 
%
%     trajchunk(num,lat,lon,32,35);
%
% load drifters
% use drifters 
% %trajchunk(num,lat,lon,32,'overlap');
% trajchunk(num,lat,lon,128);
% 
% %meanfo=zeros(length(num),1);
% %fr=zeros(length(num),1);
% 
% %dt=1/4;
% 
% %Same as below but faster
% tic
% meanfo=vmean(abs(corfreq(col2mat(cell2col(lat)))),1)'*24*dt;
% fr=frac(2*pi,dt*cellength(lat));
% toc
% 
% % tic
% % for i=1:length(num)
% %     meanfo(i)=vmean(abs(corfreq(lat{i})),1)*24*dt;
% %     fr(i)=frac(2*pi,dt*length(lat{i}));
% % end
% % toc
% 
% figure,plot(meanfo,'b.'),hold on,plot(32*fr,'ro')
% figure,plot(2*pi./(meanfo/32),'b.'),hold on,plot(2*pi./fr,'ro')
% figure,plot(2*pi./(meanfo/32),2*pi./fr,'ro'),axis equal

