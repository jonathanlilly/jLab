function[varargout]=make_jason(varargin)
%MAKE_JASON  Create a reformatted version of the Beckley altimeter dataset.
%
%   MAKE_JASON creates JasonAlongTrack, a reformatted version of the 
%   Integrated Multi-Mission Ocean Altimeter Data for Climate Research
%   Version 5.1, by Brian Beckley and Richard Ray, itself available from 
%
%     https://podaac.jpl.nasa.gov/dataset/MERGED_TP_J1_OSTM_OST_ALL_V51
% 
%   In JasonAlongTrack, altimeter passes have been rearranged and sorted
%   geographically within each cycle for convenience. Data from a given 
%   cycle can then be meaningfully plotted as a function of along-track 
%   location and track number. 
%
%   'make_jason --create' generates the dataset from source files.
%   'make_jason --f' makes a sample figure using JasonAlongTrack.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2022--2023 J.M. Lilly --- type 'help jlab_license' for details


if strcmp(varargin{1}, '--t')
    make_jason_test,return
end

if strcmp(varargin{1}, '--create')
    make_jason1, return
elseif strcmp(varargin{1}, '--f')
    make_jason_figure,return
end
 

%https://podaac.jpl.nasa.gov/dataset/MERGED_TP_J1_OSTM_OST_ALL_V51

%Download with https://github.com/podaac/data-subscriber
%podaac-data-downloader -c MERGED_TP_J1_OSTM_OST_ALL_V51 -d . --start-date 1992-09-25T00:00:00Z --end-date 2992-10-02T00:00:00Z -e ""     

%Data hosted by the PO.DAAC is openly shared, without restriction, in 
%accordance with NASA's Earth Science program Data and Information Policy.


function[]=make_jason_figure


%/*************************************************************************
%reading in jason data
jason.filename='JasonAlongTrack.nc';
jason.lat=ncread(jason.filename,'lat');
jason.lon=ncread(jason.filename,'lon');
%jason.time_offset=ncread(jason.filename,'time_offset')+datenum(1950,1,1);
jason.cycle_time=ncread(jason.filename,'cycle_time')+datenum(1950,1,1);
%jason.mss=ncread(jason.filename,'mss');
jason.time=ncread(jason.filename,'time')+datenum(1950,1,1);
jason.sla=ncread(jason.filename,'sla')*100;
jason.atd=ncread(jason.filename,'atd');
%\*************************************************************************


%/*************************************************************************
use jason
%these next two lines take a few minutes to compute the statistics, sorry
sla_std=vstd(sla,3);
sla_count=sum(isfinite(sla),3);

figure
%--------------------------------------------------------------------------
ax(1)=subplot(2,1,1);
mat=log10(100*(1-sla_count./size(sla,3)));mat(mat==2)=nan;
jpcolor(1:size(sla,2),vmean(atd,2),sla_std)
colorquant(0.5),vlines(127.5,'0.3k:'),axis tight,colormap(ax(1),flipud(crameri('davos')))
title('The JasonAlongTrack Dataset')
%--------------------------------------------------------------------------
ax(2)=subplot(2,1,2);
mat=log10(100*(1-sla_count./size(sla,3)));mat(mat==2)=nan;
jpcolor(1:size(sla,2),vmean(atd,2),mat)
caxis([log10(2) log10(90)]),vlines(127.5,'0.3k:'),axis tight,colormap(ax(2),flipud(crameri('davos')))
xlabel(['(Descending Tracks 1--127)\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,' ...
    'Track Number\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,(Ascending Tracks 128--254)'])
%--------------------------------------------------------------------------
ax=packfig(2,1,'rows');
for i=1:2
   axes(ax(i))
   %text(245,18,['(' char( real('a')+i-1) ')'])
   hc=colorbar('EastOutside');
   ylabel('Along-Track Distance (Mm)')
   if i==1
        hc.Label.String='SSH Standard Deviation (cm)';
       hc.Ticks=[0:5:50];
   elseif i==2
       hc.Label.String='Percent of Missing Data';
       hc.Ticks=log10([1 2 5 10 20 50 90]);
       hc.TickLabels={'1','2','5','10','20','50','90'};
   end
   %plot([92+0*1i 92+2*1i],'w','linewidth',2)
   %plot([92+0*1i 92+2*1i],'k','linewidth',1.5)
   %plot([92+maxmax(atd)*1i 92+(maxmax(atd)-2)*1i],'w','linewidth',2)
   %plot([92+maxmax(atd)*1i 92+(maxmax(atd)-2)*1i],'k','linewidth',1.5)
   set(gcf,'color','w');set(gca,'color',0.6*[1 1 1]);set(gcf, 'InvertHardcopy', 'off')
   vlines(127.5,'0.3k')
   xticks([25:25:300]),yticks([2:2:18])%,fixlabels([0,-1])
   %yticks([0:.5:3.5]),
   posax=get(gca,'Position');
   pos=hc.Position;
   hc.Position=[pos(1) pos(2)+0.5*pos(3) pos(3)/2 pos(4)-pos(3)];
   set(gca,'Position',posax)%thin x colorbar
   set(gca,'TickLen',get(gca,'TickLen')/2)
end
%--------------------------------------------------------------------------
fontsize 8 8 8 8
%set(gcf,'paperposition',[1 1 8.5 9]),
set(gcf,'paperposition',[1 1 7 5]),
jprint(pwd,'jason_distribution_twopanels','png','-r300')
%\*************************************************************************


function[]=make_jason1
 
disp('Creating JasonAlongTrack.  This may take while...')

tic;% Track how long the whole thing takes
addpath('/Volumes/PlanetaryScience/beckley')
%addpath('/Users/lilly/Downloads')
filename='Merged_TOPEX_Jason_OSTM_Jason-3_Version_V5.1.nc';

ncdisp(filename)

lat=double(ncread(filename,'lat')');
%lon=deg180(double(ncread(filename,'lon')'));
lon=double(ncread(filename,'lon')');
mss=double(ncread(filename,'mssh')')/1000;  %Convert mm to m
topo=double(ncread(filename,'Bathymetry')')/1000;%Convert m to km  
dfc=double(ncread(filename,'Distance_to_coast')')/1000; %Convert m to km
stf=double(ncread(filename,'Surface_Type')');

vinfo=ncinfo(filename,'time');
N=vinfo.Size(1);
K=vinfo.Size(2);  

rev=vrep([1:127],size(lat,1),1);  %revolution number


% indexo=[1:length(lat(:))]';
% vcolon(lat,lon,mss)
% ismaxlat=lat>vshift(lat,1,1)&lat>vshift(lat,-1,1);
% vindex(lat,lon,indexo,ismaxlat,[1686:length(lat) 1:1685],1);

%/*************************************************************************
vcolon(lat,lon,mss,topo,dfc,stf,rev);
newindex=(1:length(lat))';
%--------------------------------------------------------------------------
%remove four bad points in each cycle 

vindex(lat,lon,mss,topo,dfc,stf,rev,newindex,[1687:length(lat) 1:1686],1);

length(find(isnan(lat(3375+1685:6750:end))))%31
length(find(isnan(lat(3375+1686:6750:end))))%127
length(find(isnan(lat(3375+1687:6750:end))))%127
length(find(isnan(lat(3375+1688:6750:end))))%127
length(find(isnan(lat(3375+1689:6750:end))))%127
length(find(isnan(lat(3375+1690:6750:end))))%0

%so there is a set of say, five (or maybe four) bad points 

mean(lat(3375+1690:6750:end)-lat(3375+1684:6750:end),'omitnan')%0.086265038960159
mean(lat(3375+1690:6750:end)-lat(3375+1685:6750:end),'omitnan')%0.049236635334789
mean(lat(3375+1691:6750:end)-lat(3375+1690:6750:end),'omitnan')%0.049047062924227
%ok, this means I should keep the 31 bad points to have the right delta-lat

bool=true(size(lat));
bool(3375+1686:6750:end)=false;
bool(3375+1687:6750:end)=false;
bool(3375+1688:6750:end)=false;
bool(3375+1689:6750:end)=false;

vindex(lat,lon,mss,topo,dfc,stf,rev,newindex,bool,1);
%--------------------------------------------------------------------------
%check locations of max and min
ismaxlat=lat(:)>=vshift(lat(:),1,1)&lat(:)>=vshift(lat(:),-1,1);
%how many data points to the nearest maximum?

N=3373;
ismaxm1=length(find(ismaxlat(2*N-1:2*N:end))) %0
ismaxm1=length(find(ismaxlat(2*N:2*N:end))) %53
ismax0=length(find(ismaxlat(1:2*N:end))) %83  %<-- expected
ismaxp1=length(find(ismaxlat(2:2*N:end))) %0
%so, all maxima are within one data point

isminlat=lat(:)<=vshift(lat(:),1,1)&lat(:)<=vshift(lat(:),-1,1);
%how many data points to the nearest minimum?

ismaxm2=length(find(isminlat(N-1:2*N:end))) %0
ismaxm1=length(find(isminlat(N:2*N:end))) %44 
ismax0=length(find(isminlat(N+1:2*N:end))) %93<-- expected
ismaxp1=length(find(isminlat(N+2:2*N:end))) %0

%so, all minima are also within one data points
%--------------------------------------------------------------------------
%reshape
lat=reshape(lat,N,254);
lon=reshape(lon,N,254);
mss=reshape(mss,N,254);
topo=reshape(topo,N,254);
dfc=reshape(dfc,N,254);
stf=reshape(stf,N,254);
rev=reshape(rev,N,254);
newindex=reshape(newindex,N,254);
%--------------------------------------------------------------------------
%sort by initial longitude
a=lon(1,:);
[~,sorter]=sort(a);
vindex(lat,lon,mss,topo,dfc,stf,rev,newindex,sorter,2);

%sort into descending and ascending
vindex(lon,lat,mss,topo,dfc,stf,rev,newindex,[1:2:254 2:2:254],2);

%flip first half
lat(:,1:127)=flipud(lat(:,1:127));
lon(:,1:127)=flipud(lon(:,1:127));
mss(:,1:127)=flipud(mss(:,1:127));
topo(:,1:127)=flipud(topo(:,1:127));
dfc(:,1:127)=flipud(dfc(:,1:127));
stf(:,1:127)=flipud(stf(:,1:127));
rev(:,1:127)=flipud(rev(:,1:127));
newindex(:,1:127)=flipud(newindex(:,1:127));

% %Sort by initial longitude
% [~,sorter]=sort(lon(1,1:127));
% vindex(lat,lon,mss,topo,dfc,stf,rev,newindex,[sorter 128:254],2);
% [~,sorter]=sort(lon(1,128:254));
% vindex(lat,lon,mss,topo,dfc,stf,rev,newindex,[1:127 sorter+127],2);

%fill bad points 
lat=fillbad(lat);
lon=fillbad(lon);
mss=fillbad(mss);
topo=fillbad(topo);
dfc=fillbad(dfc);
stf=fillbad(stf,'nearest');
%\*************************************************************************

%figure,plot(lon(1,:),'.')  %sorted by descending / ascending and by lon

%/*************************************************************************
%Test of index
lato=double(ncread(filename,'lat')');
lono=double(ncread(filename,'lon')');

lat2=nan*zeros(size(lat));
lon2=nan*zeros(size(lat));
lat2=lato(newindex);
%lon2(nani)=deg180(lono(newindex(nani)));
lon2=lono(newindex);

nani=~isnan(lon2);

reporttest('Reshaping index matches expectation for latitude and longitude',aresame(lon(nani),lon2(nani))&&aresame(lat(nani),lat2(nani)))
%\*************************************************************************

%/*************************************************************************
%Read in time, ssh, and flag, with reshaping

%From ncdisp(filename)
%Seconds from epoch 01/01/1992 00:00:00 (MJD = 48621.0)
%Conversion to Modifed Julian Day + fraction of day (accurate to a second): MJD = time/86400.0 + 48621.0
%
%  :time_coverage_start = "1992-09-25T04:49:31.0";
%  :time_coverage_end = "2021-06-23T20:06:23.0";

numo=double(ncread(filename,'time'))./86400+datenum(1992,1,1)-datenum(1950,1,1);
ssho=double(ncread(filename,'ssha'))/1000; %Convert mm to m
flago=double(ncread(filename,'flag'));

% accurate to within a second
%[datestr(minmin(numo+datenum(1992,1,1))),datestr(maxmax(numo+datenum(1992,1,1)))]
[datestr(minmin(numo+datenum(1950,1,1))),datestr(maxmax(numo+datenum(1950,1,1)))]
%'25-Sep-1992 04:49:31' '23-Jun-2021 20:06:23'

numall=nan*zeros(size(lat,1),size(lat,2),K);
ssh=nan*zeros(size(lat,1),size(lat,2),K);
flag=zeros(size(lat,1),size(lat,2),K);

ii=double(ncread(filename,'index'));
jj=double(ncread(filename,'reference_orbit'));
index=sub2ind(size(ncread(filename,'lat')'),ii,jj);

for k=1:K
    k
    
    %Map alongtrack data into Beckley format matrix
    numk=nan*zeros(6750,127);
    numk(index)=numo(:,k);

    sshk=nan*zeros(6750,127);
    sshk(index)=ssho(:,k);
    
    flagk=zeros(6750,127);
    flagk(index)=flago(:,k);

    %Map Beckley format into Lilly format
    newnumk=nan*zeros(size(lat));
    newnumk(nani)=numk(newindex(nani));
    numall(:,:,k)=newnumk;
    
    newsshk=nan*zeros(size(lat));
    newsshk(nani)=sshk(newindex(nani));
    ssh(:,:,k)=newsshk; 

    newflagk=zeros(size(lat));
    newflagk(nani)=flagk(newindex(nani));    
    flag(:,:,k)=newflagk;
end
vsize(flag,numall,ssh)
datestr(minmin(numall)+datenum(1950,1,1))
datestr(maxmax(numall)+datenum(1950,1,1))
clear numo ssho flago
%\*************************************************************************

%/*************************************************************************
% Fix a problem I found with time in cycle 269.  This was reported and 
% confirmed by Brian Beckley and should be taken care of in the next release
if 0
    time1=nan*zeros(6750,127);
    timek=nan*zeros(6750,127);
    timen=nan*zeros(6750,127);
    time1(index)=double(ncread(filename,'time',[1 268],[N 1]));
    timek(index)=double(ncread(filename,'time',[1 269],[N 1]));
    timen(index)=double(ncread(filename,'time',[1 270],[N 1]));

    figure
    plot(timek(:,1:end-1)-time1(:,1:end-1)),hold on
    plot(timek(:,end)-time1(:,end))

    %no such offset in 270-268
    figure
    plot(timen(:,1:end-1)-time1(:,1:end-1)),hold on
    plot(timen(:,end)-time1(:,end))
    %plot(timen-time1),hold on

    figure,
    plot([timek(:,end) timen(:,1)]/60-dt*24*60-minmin(time1(:,end))/60,'o'),hold on
    plot([time1(:,end) timek(:,1)]/60-minmin(time1(:,end))/60,'linewidth',2),hold on
    xlabel('Along-pass location')
    ylabel('Time (minutes) relative to start of Pass 127')
    hlines(minmin(timek(:,1))/60-minmin(time1(:,end))/60,'k:')
    legend('Cycle 269 Pass 127','Cycle 270 Pass 1','Cycle 268 Pass 127','Cycle 269 Pass 1','Start of Pass 1')
    title('Time problem with Cycle 269 Pass 127')
    axis tight
    fontsize 14 14 14 14
    jprint(pwd,'jason-time-problem')

    max(time1(:))-min(timek(:))%about 53 seconds
    max(timek(:))-min(timen(:))%about 385 seconds = 6 minute
end
%time1=double(ncread(filename,'time',[1 268],[N 1]));  
%timek=double(ncread(filename,'time',[1 269],[N 1]));  
%timen=double(ncread(filename,'time',[1 270],[N 1]));  

%shifted time in cycle 269 is associated with very last track
%up to 350 seconds off from normal  

%figure,plot([time1 timek-dt*24*3600 timen-2*dt*24*3600])
%figure,plot([timek;timen]/24/3600-dt0), hold on,
%plot([time1;timek]/24/3600)
%vlines(5.98128*1e5), vlines(length(time1))
%There appears to be a problem with the original time file...
%it does not make any sense that time should slow down, only 
%to catch back up again later. Thus, this appears to be an error

%dt0=median(vcolon(numall(:,:,269)-numall(:,:,268)),'omitnan');
%plot((numall(:,:,269)-numall(:,:,268)-dt0)*24*3600)
% indexk=find((numall(:,:,269)-numall(:,:,268)-dt0)*24*3600<-1);
% dtime=(numall(:,:,269)-numall(:,:,268)-dt0)*24*3600;
% this is happening at tracks 12, 197, and 207
% figure,plot(dtime),hold on,plot(dtime(:,[12 197 207]),'ko')
%
% plot((numall(:,12,269)-numall(:,12,268)-dt0)*24*3600)
% plot(rev(:,12))
% plot((numall(:,197,269)-numall(:,197,268)-dt0)*24*3600)
% plot(rev(:,197))
% plot((numall(:,207,269)-numall(:,207,268)-dt0)*24*3600)
% plot(rev(:,207))
% %so the problem is entirely associated with revolution #127

%length(find(~isnan(numall(:,12,268))))
%length(find(~isnan(numall(:,12,269))))
%length(find(~isnan(numall(:,197,268))))
%length(find(~isnan(numall(:,197,269))))
%length(find(~isnan(numall(:,207,268))))
%length(find(~isnan(numall(:,207,269))))
%these lengths are all the same, so I can just sub dt
dt0=median(vcolon(numall(:,:,269)-numall(:,:,268)),'omitnan');
%numall(:,[12 197 207],269)=numall(:,[12 197 207],268)+dt0;
%numall(:,[70 80 139],269)=numall(:,[70 80 139],268)+dt0;%180 version 
%numall(:,[7 17 202],269)=numall(:,[7 17 202],268)+dt0;%360 version
numall(:,[75 133 143],269)=numall(:,[75 133 143],268)+dt0;%360 version
%jpcolor(numall(:,:,269)-numall(:,:,268))
%\*************************************************************************

%/*************************************************************************
%Additional processing for time variable

%first, find the slices with the least data dropout
ngood=sum(isfinite(reshape(numall,size(numall,1)*size(numall,2),K)),1);
max(ngood)%602834
goodindex=find(ngood==max(ngood));
length(goodindex)%990

mnum=zeros(size(goodindex));
for i=1:length(goodindex)
   mnum(i)=mean(vcolon(numall(:,:,goodindex(i))),'omitnan');
end
dmnum=diff(mnum);
dt=mean(dmnum(dmnum<11))%9.915644056430926
length(find(dmnum<11)) %911
std(dmnum(dmnum<11))./dt%4.711066134136858e-06
std(dmnum(dmnum<11))./dt*3600*24% about a 1/2 second standard deviation

%my reference time should go right through the middle of the 883 cycles 
%that have the maximum number of good points and are adjacent to others 

num=(0:K-1)'*dt;
dnum=numall-vrep(vrep(permute(num,[3 2 1]),size(numall,1),1),size(numall,2),2);

num0=mean(vcolon(dnum(:,:,goodindex(dmnum<11))),'omitnan');
num=num+num0;
dnum=dnum-num0;

%find series of mean values of dnum for goodindex
% mdnum=zeros(size(goodindex));
% for i=1:length(goodindex)
%    mdnum(i)=mean(vcolon(dnum(:,:,goodindex(i))),'omitnan');
% end
%plot(24*3600*mdnum)
%the average over the whole cycle agrees with num0 to within 
%about 80 seconds whenever I have the max number of observations

minnum=min(reshape(numall,size(numall,1)*size(numall,2),K),[],1)';
maxnum=max(reshape(numall,size(numall,1)*size(numall,2),K),[],1)';

%figure,plot(maxnum-num),hold on,plot(minnum-num),
%ok, we are in the middle since we're less than 5 days on either side
%max(maxnum-num)-min(minnum-num)

num_uniform=num;

%Now I refine num as the mean value at each slice
dnum=mean(dnum,3,'omitnan');%average difference from the mean time at each k
for k=1:K
    num(k)=mean(vcolon(numall(:,:,k)-dnum),'omitnan');
end

%Note that cycle 521 has no good data and thus cycle_time is undefined.
%Thanks to André Palóczy for pointing this out.
%If desired could write num(521)=0.5*num(520)+0.5*num(522);

ddnum=nan*zeros(size(numall));
for k=1:K
   ddnum(:,:,k)=numall(:,:,k)-num(k)-dnum;
   %ddnum(:,:,k)=numall(:,:,k)-num_uniform(k)-dnum;
end

%what is the standard deviation
mddnum=nan*num;
maxddnum=nan*num;
stdnum=nan*num;

% for k=1:K
%    mddnum(k)=mean(vcolon(ddnum(:,:,k)),'omitnan');
%    stdnum(k)=std(vcolon(ddnum(:,:,k)),'omitnan');
%    maxddnum(k)=max(abs(vcolon(ddnum(:,:,k))));
% end
%figure,plot([maxddnum stdnum]*3600*24),ylim([0 1]*2)
[sqrt(mean(squared(ddnum(:)),'omitnan')) maxmax(abs(ddnum(:)))]*3600*24
%rms error from num+dnum is 0.30 seconds and max is 3.8 seconds
clear ddnum
%\*************************************************************************


%%/*************************************************************************
%%convert flagword to boolean array
%%maximum value of flag is 32767. This is equal to sum(2.^[0:14])
%%in other words I have 15 different booleans
%%see about bitwise flags, https://mplnet.gsfc.nasa.gov/about-flags

% using bitget
% boolflag=false(size(flag,1),size(flag,2),size(flag,3),15);
% for i=1:15
%     boolflag(:,:,:,i)=bitget(flag,i);
% end

% %the hard way
% boolflag=false(size(flag,1)*size(flag,2)*size(flag,3),15);
% flagn=flag;
% for i=15:-1:1
%     i
%     index=find(floor(flagn./2.^(i-1))==1);
%     boolflag(index,i)=true;
%     flagn(index)=flagn(index)-2.^(i-1);
% end
% %b15 = bitget(flag,15);
% %aresame(b15,boolflag(:,:,:,15))
% %maxmax(flagn)==0 %true
% %boolflag=reshape(boolflag,size(flag,1),size(flag,2),size(flag,3),15);
%clear flagn flag
%jpcolor(mean(boolflag(:,:,:,1),5,'omitnan'))
%%\*************************************************************************


if 0
%/*************************************************************************
%compute my own version of depth from Smith and Sandwell
[topo2,lattopo,lontopo]=readtopo([-180 180 -80 80]);
dlon=lontopo(2)-lontopo(1);

lontopo=[lontopo(1)-2*dlon lontopo(1)-dlon lontopo lontopo(end)+dlon lontopo(end)+2*dlon];
topo2=[topo2(:,end-1) topo2(:,end) topo2 topo2(:,1) topo2(:,2)];

depth=-interp2(lontopo,lattopo,topo2,lon,lat);
boolland=depth<0; %land according to Smith and Sandwell

%similar to, but not identical with, the surface type land flag from Beckley
%jpcolor(boolland-(stf==1|stf==3))

%Need to put this into a surface type flag... maybe another boolean

%could also put high std, and low % coverage, as additional flags
%\*************************************************************************
end


atd=spheredist(lat,lon)/1000;
%plot(atd(end,:)-atd(end,1),'.')
%this is not completely right!!! ... not sure what this means
%/*************************************************************************
if exist('tpjaos_rejects.mat')
    load tpjaos_rejects
    use tpjaos_rejects
    ssh(rejected_index)=nan;
else
    %despiking with running mean absolute deviation
    ssh0=ssh;

    L=21;
    rmax=6;%cutoff ratio
    nmin=10; %cutoff number of data points

    tic
    parfor j=1:size(ssh,2)
        j
        sshj=squeeze(ssh(:,j,:));
        for i=1:3
            medsshj=vfilt(sshj,L,'median','mirror');
            devj=abs(sshj-medsshj);
            meddevj=vfilt(devj,L,'median','mirror');
            devratj=devj./meddevj;
            ngoodj=vfilt(isfinite(sshj),ones(L,1),'mirror');
            %length(find((devratj>=rmax|ngoodj<=nmin)&isfinite(sshj)))
            sshj(devratj>=rmax|ngoodj<=nmin)=nan;
        end
        ssh(:,j,:)=sshj;
    end
    toc %24 minutes
    rejected_index=find(isnan(ssh)&~isnan(ssh0));
    rejected_values=ssh0(rejected_index);
    clear ssh0

    matsave tpjaos_rejects rejected_index rejected_values
    %ssh(rejected_index)=nan;
    % dfcref=vrep(dfc,size(ssh,3),3);
    % [n,x]=hist(dfcref(bool),[0:10:maxmax(dfc)]);
    % plot(x,cumsum(n)./sum(n)),vlines(100)
    % %so one-third are within 100~km of the coasts
end
%\*************************************************************************

if 1 %skip the wavelet noise estimate for now
%/*************************************************************************
if 0%exist('tpjaos_noiselevel.mat')
    load jason_noiselevel
    use jason_noiselevel
else
    allarebad=sum(isfinite(ssh),3)<1;
    %Estimating noise standard deviation
    wasfilled=~isfinite(ssh);

    be=2;ga=3;
    fs=morsespace(ga,be,{0.1,pi},2*pi/200,4);
    %--------------------------------------------------------------------------
    %highest-frequency transform
    jj=1;
    psi=morsewave(100,ga,be,fs(jj));
    %[psi,psif]=morsewave(100,ga,be,fs);
    %figure,uvplot(psi./max(psi)),hold on,plot(abs(psi./max(psi)))

    wtrans=nan*zeros(size(ssh));
    wfilled=nan*zeros(size(ssh));%wavelet transform of filled values
    tic
    for k=1:size(ssh,3)
        k
        sshk=vswap(fillbad(ssh(:,:,k),nan,inf),nan,0);
        wtrans(:,:,k)=wavetrans(sshk,{ga,be,fs(jj),'bandpass'});
        wfilled(:,:,k)=abs(wavetrans(double(wasfilled(:,:,k)),{ga,be,fs(jj),'bandpass'}));
    end
    toc
    %3.5 minutes

    %remove portions that experienced filled data
    100*length(find(wfilled./max(abs(psi))>0.1))./N %  13.970629925294961
    100*length(find(wfilled./max(abs(psi))>0.2))./N %  10.877284988749047

    % %Noise level estimated by after removing only initially bad data points
    % wtrans(wasfilled)=nan;
    % %The standard deviation in the trans wavelet band
    % sigma=sqrt(vmean(squared(wtrans(:)),1))  %6.422695395567049
    % %The corresponding estimated signal standard standard deviation
    % sigmaeps=sigma*sqrt(frac(morsefreq(ga,be),fs(1)*ffun))%7.340155396404459

%     %Noise level estimated by after removing contaminated points at 0.2 level
%     wtrans(wasfilled|wfilled./max(abs(psi))>0.2)=nan;
%     %The standard deviation in the trans wavelet band
%     sigma0=sqrt(vmean(squared(wtrans(:)),1))%2.349409003386438
%     %The corresponding estimated signal standard standard deviation
%     ffun=morseffun(ga,be,0);%See (4.10) of Lilly (2017)
%     sigma=sigma0*sqrt(frac(morsefreq(ga,be),fs(jj)*ffun))%2.68501401864250


    %exploring median.  But, I don't think I want to do this because the
    %transform already bandpasses it
    %sigma_trans=1.4826*vmedian(abs(wtrans),3);
    %how to convert the median of a chi-squared variable 
    %https://en.wikipedia.org/wiki/Median_absolute_deviation
    %https://dsp.stackexchange.com/questions/25321/estimating-the-variance-of-noise-using-median-function-applies-on-the-gradient-o
    %[mat,xmid,ymid]=twodhist(sigma_trans_mean*100,sigma_trans*100,[0:.01:5],[0:.01:5]);
    %jpcolor(xmid,ymid,log10(mat))

    %Noise level estimated by after removing contaminated points at 0.2 level
    wtrans(wasfilled|wfilled./max(abs(psi))>0.2)=nan;

    sigma_trans=sqrt(vmean(squared(wtrans),3));
    ffun=morseffun(ga,be,0);%See (4.10) of Lilly (2017)
    sigma=sigma_trans*sqrt(frac(morsefreq(ga,be),fs(jj)*ffun));

%     %is there a tendency for more rarely-sampled locations to have higher
%     %(possibly spurious) noise levels? No
%     pctgood=vsum(isfinite(ssh),3)./size(ssh,3);
%     [mat,xmid,ymid]=twodhist(pctgood,sigma*100,[0:.01:1],[0:.05:50]);
%     jpcolor(xmid,ymid,log10(mat))
% 
%     %The higher noise levels are near the coasts. Very interesting!
%     [mat,xmid,ymid]=twodhist(dfc,sigma*100,[0:.01:1]*1000,[0:.05:50]);
%     jpcolor(xmid,ymid,log10(mat))

    sigma(allarebad)=nan;

    matsave jason_noiselevel be ga fs sigma
end
%\*************************************************************************
end

%/*************************************************************************
%Interpolate to get mdt from dtu.  Thanks to André Palóczy for this idea.

if 0
dtulat=double(ncread('DTU15MSS_1min.nc','lat'));
dtulon=double(ncread('DTU15MSS_1min.nc','lon'));
dtumss=double(ncread('DTU15MSS_1min.nc','mss'))';
mss2=interplatlon(dtulat,dtulon,[],dtumss,lat,lon,'cubic');

%How close is mss2 to mss?  mss is the dtu15 mean sea surface provided 
%by Beckley, and mss2 is the one I compute via interpolation
mssdiff=abs(mss2-mss)*1000;
mssdiff(mssdiff>1000)=nan; 
%jpcolor(mssdiff),caxis([0 1])
median(mssdiff(isfinite(mssdiff))) %median absolute deviation in mm 
%0.5 mm, less than one mm.  Greater if I use linear. Looks interpolation noise. 
end

dtulat=double(ncread('DTU15MDT_1min.mdt.nc','lat'));
dtulon=double(ncread('DTU15MDT_1min.mdt.nc','lon'));
dtumdt=double(ncread('DTU15MDT_1min.mdt.nc','mdt'))';
mdt=interplatlon(dtulat,dtulon,[],dtumdt,lat,lon,'cubic');
%\*************************************************************************

%/*************************************************************************
%start saving the NetCDF file

%Some useful links
%https://foundations.projectpythia.org/core/data-formats/netcdf-cf.html
%http://cfconventions.org/cf-conventions/cf-conventions.html
%https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3
%https://www.unidata.ucar.edu/software/udunits/udunits-2.2.28/udunits2.html#Database
%https://github.com/Unidata/netcdf4-python/issues/442  %no reference year 0

%coverage_content_type	An ISO 19115-1 code to indicate the source of the data 
%(image, thematicClassification, physicalMeasurement, auxiliaryInformation, 
%qualityInformation, referenceInformation, modelResult, or coordinate).
%https://ngdc.noaa.gov/wiki/index.php/ISO_19115_and_19115-2_CodeList_Dictionaries#MD_CoverageContentTypeCode

writedir=['/Users/lilly/Desktop/Dropbox/NetCDF/'];
name_nc='JasonAlongTrack.nc';
%--------------------------------------------------------------------------
% Define the global attributes

ncid = netcdf.create([writedir name_nc],'NETCDF4');
varid = netcdf.getConstant('GLOBAL');

%GeosatAlongTrack, ERSAlongTrack, etc. later 

%Highly recommended
netcdf.putAtt(ncid,varid,'title','JasonAlongTrack: A reformatted version of the Integrated Multi-Mission Ocean Altimeter Data for Climate Research Version 5.1');
netcdf.putAtt(ncid,varid,'summary','Geo-registered along-track sea surface height anomalies with respect to the DTU15 mean sea surface at 1-second intervals from Jason-class altimeters, reformatted for convenience into a 3D array with dimensions of along-track direction by geographically sorted track number by cycle.');
netcdf.putAtt(ncid,varid,'id','JasonAlongTrack');
netcdf.putAtt(ncid,varid,'keywords','Oceans, Sea Surface Topography, Sea Surface Height');
netcdf.putAtt(ncid,varid,'Conventions','CF-1.9, ACDD-1.3');

%Recommended and suggested
netcdf.putAtt(ncid,varid,'product_version','1.0.2');  %suggested
netcdf.putAtt(ncid,varid,'naming_authority','edu.psi');
netcdf.putAtt(ncid,varid,'history','05-Nov-2023 version 1.0.2 created adding noise estimate and modified documentation | 02-Nov-2023 version 1.0.1 created adding uniform_time variable | 30-Oct-2023 reformatted version 1.0 created');
netcdf.putAtt(ncid,varid,'source','TOPEX/Poseidon MGDR_B: Benada, J.R. 1997. PO.DAAC Merged GDR (TOPEX/POSEIDON) Generation B Users Handbook, Version 2.0 JPL D-11007;Jason-1 GDR_E: AVISO and PODAAC User Handbook. IGDR and GDR Jason Products SMM-MU-M5-OP-13184-CN (AVISO), JPL D-21352 (PODAAC) Edition 4;OSTM GDR_D: OSTM/Jason-2 Products Handbook, CNES: SALP-MU-M-OP-15815-CN, EUMETSAT: EUM/OPS-JAS/MAN/08/0041, JPL: OSTM-2 9-1237, NOAA/NESDIS: Polar Series/OSTM J400, Issue 1 rev 8, December 1, 2011;Jason-3 GDR_D: Jason-3 Products Handbook, CNES : SALP-MU-M-OP-16118-CN, Issue: 2 rev 1, June 21, 2021');
netcdf.putAtt(ncid,varid,'processing_level','2');
netcdf.putAtt(ncid,varid,'featureType','trajectory');
netcdf.putAtt(ncid,varid,'comment','This is a reformatted version of the Brian Beckley and Richard Ray''s ''Integrated Multi-Mission Ocean Altimeter Data for Climate Research complete time series Version 5.1'', available from https://podaac.jpl.nasa.gov/dataset/MERGED_TP_J1_OSTM_OST_ALL_V51. The changes are as follows.  Altimeter passes are sorted according to their initial longitude, then split into descending and ascending potions with all descending tracks preceding all ascending tracks.  Descending tracks are then flipped so that latitude increases in the alongtrack direction for all tracks.  This leads to a 3373 x 254 matrix of observational locations, with the first dimension being the along-track location and the second dimension being the track index.  Sea surface height anomaly, time, and flag values are then placed into their correct locations within this matrix, such that these three variables are all of size 3373 x 254 x K where K is the number of cycles, currently 1087.  A very good approximation to the time at each of the 3373 x 254 x K observation points is constructed with a length K array of cycles times together with a 3373 x 254 array of time offsets.  A median-based editing criterion in introduced to identify a small number of suspect data points.  These are set to a value of NaN in sla, but their positions and values are recorded in rejected_index and rejected_values, respectively. The DTU15 mean dynamic topography (mdt) is included, in addition to the mean sea surface field already provided, interpolated onto the track locations using bicubic interpolation.  Finally, an estimate of the small-scale noise level, sigma, is produced using a wavelet transform filter.');
%Finally, an estimate of the small-scale noise level, sigma, is produced using a wavelet transform filter.']);
% Finally, a despiking flag is included.
netcdf.putAtt(ncid,varid,'acknowledgment','The creation of this NetCDF data product was carried out by J. M. Lilly with the support of NASA grant 80NSSC21K1823.');
netcdf.putAtt(ncid,varid,'license','Creative Commons Attribution 4.0 International (CC BY 4.0), https://creativecommons.org/licenses/by/4.0/');
netcdf.putAtt(ncid,varid,'standard_name_vocabulary','CF Standard Name Table v79');
netcdf.putAtt(ncid,varid,'date_created','2023-10-30T');
%netcdf.putAtt(ncid,varid,'creator_name','J. M. Lilly');
%netcdf.putAtt(ncid,varid,'creator_url','https://www.jmlilly.net');
%netcdf.putAtt(ncid,varid,'creator_email','jmlilly@psi.edu');
%netcdf.putAtt(ncid,varid,'creator_institution','Planetary Science Institute'); %suggested
netcdf.putAtt(ncid,varid,'creator_name','Beckley, B.; Zelensky, N.P.; Holmes, S.A.; Lemoine, F.G.; Ray, R.D.; Mitchum, G.T.; Desai, S.; Brown, S.T.');
netcdf.putAtt(ncid,varid,'creator_url','https://podaac.jpl.nasa.gov/dataset/MERGED_TP_J1_OSTM_OST_CYCLES_V51');
%netcdf.putAtt(ncid,varid,'creator_type','person');
netcdf.putAtt(ncid,varid,'creator_type','institution');%as listed in the original
netcdf.putAtt(ncid,varid,'creator_email','Brian.D.Beckley@nasa.gov; Richard.D.Ray@nasa.gov');
netcdf.putAtt(ncid,varid,'creator_institution','KBR, Inc.; NASA/Goddard Space Flight Center'); %suggested
netcdf.putAtt(ncid,varid,'institution','NASA/Goddard Space Flight Center');
netcdf.putAtt(ncid,varid,'project','Making Earth System Data Records for Use in Research Environments (MEaSUREs, NASA grant NNH06ZDA001N); Eddy Dynamics from Along-Track Altimetry (NASA grant 380NSSC21K182)');
netcdf.putAtt(ncid,varid,'publisher_name','J. M. Lilly');
netcdf.putAtt(ncid,varid,'publisher_email','jmlilly@psi.edu');
netcdf.putAtt(ncid,varid,'publisher_url','http://www.jmlilly.net');
netcdf.putAtt(ncid,varid,'publisher_institution','Planetary Science Institute'); %suggested
netcdf.putAtt(ncid,varid,'contributor_name','J. M. Lilly');
netcdf.putAtt(ncid,varid,'contributor_role','reformatted dataset and added median-based editing criterion, DTU15 mean dynamic topography, and small-scale noise estimate');
%netcdf.putAtt(ncid,varid,'contributor_name','Beckley, B.; Zelensky, N.P.; Holmes, S.A.; Lemoine, F.G.; Ray, R.D.; Mitchum, G.T.; Desai, S.; Brown, S.T.');
%netcdf.putAtt(ncid,varid,'contributor_role','created the multi-mission along-track dataset upon which this reformatted dataset is based');
netcdf.putAtt(ncid,varid,'geospatial_lat_min',minmin(lat));
netcdf.putAtt(ncid,varid,'geospatial_lat_max',maxmax(lat));
%netcdf.putAtt(ncid,varid,'geospatial_lat_resolution','0.25 degree'); %suggested
netcdf.putAtt(ncid,varid,'geospatial_lon_min',minmin(lon));
netcdf.putAtt(ncid,varid,'geospatial_lon_max',maxmax(lon));
%netcdf.putAtt(ncid,varid,'geospatial_lon_resolution','0.25 degree');  %suggested
netcdf.putAtt(ncid,varid,'time_coverage_start','1992-09-25T04:49:31');
netcdf.putAtt(ncid,varid,'time_coverage_end','2022-02-28T11:25:06');%make sure to update this
%netcdf.putAtt(ncid,varid,'time_coverage_duration','P0001-00-00T00:00:00');%one year
netcdf.putAtt(ncid,varid,'time_coverage_resolution','P0000-00-00T00:00:01');%one second
netcdf.putAtt(ncid,varid,'references','Lemoine, et. al. 2010, Advances in Space Research, Beckley et. al. 2010, Marine Geodesy');%suggested

%Additional attributes included from original dataset
netcdf.putAtt(ncid,varid,'reference_document','Integrated Multi-Mission Ocean Altimeter Data for Climate Research TOPEX/Poseidon, Jason-1, OSTM/Jason-2 and Jason-3 Users Handbook, Version 5.1, https://doi.org/10.5067/ALTUG-TJ151.');
netcdf.putAtt(ncid,varid,'platform','Integrated T/P, Jason-1, OSTM/Jason-2 and Jason-3 Altimetry');
netcdf.putAtt(ncid,varid,'instrument','TOPEX-A, TOPEX-B, Poseidon-1, Poseidon-2, Poseidon-3, Poseidon-3B');
netcdf.putAtt(ncid,varid,'radiometer_sensor_name', 'TMR (TOPEX/Poseidon Microwave Radiometer), JMR (Jason-1 Microwave Radiometer), AMR (OSTM Advanced Microwave Radiometer), AMR (Jason-3 Advanced Microwave Radiometer)');

%Notes to self: these are not needed because they are assumed by default
%netcdf.putAtt(ncid,varid,'geospatial_lat_units','degree_north');
%netcdf.putAtt(ncid,varid,'geospatial_lon_units','degree_east');

% Define the dimensions
dimitime = netcdf.defDim(ncid,'cycle_time',length(num));
dimitrack = netcdf.defDim(ncid,'track_number',size(lat,2));
dimialong = netcdf.defDim(ncid,'along_track',size(lat,1));
dimirej = netcdf.defDim(ncid,'rejected_points',size(rejected_index,1));

time_ID = netcdf.defVar(ncid,'time','double',[dimialong dimitrack dimitime]);
netcdf.putAtt(ncid,time_ID,'standard_name','time');
netcdf.putAtt(ncid,time_ID,'long_name','time of observations');
netcdf.putAtt(ncid,time_ID,'coordinates','lat lon cycle_time');
netcdf.putAtt(ncid,time_ID,'units','days since 1950-01-01 00:00:00');
netcdf.putAtt(ncid,time_ID,'time_zone','UTC') ;
netcdf.putAtt(ncid,time_ID,'calendar','standard');
netcdf.putAtt(ncid,time_ID,'coverage_content_type','referenceInformation');
netcdf.defVarFill(ncid,time_ID,false,nan);

uniform_time_ID = netcdf.defVar(ncid,'uniform_time','double',dimitime);
netcdf.putAtt(ncid,uniform_time_ID,'standard_name','time');
netcdf.putAtt(ncid,uniform_time_ID,'long_name','a uniformly progressing time variable for each cycle that increases at the average inter-cycle interval of 9.915644056430926 days');
netcdf.putAtt(ncid,uniform_time_ID,'units','days since 1950-01-01 00:00:00');
netcdf.putAtt(ncid,uniform_time_ID,'axis','T')
netcdf.putAtt(ncid,uniform_time_ID,'time_zone','UTC') ;
netcdf.putAtt(ncid,uniform_time_ID,'calendar','standard');
netcdf.putAtt(ncid,uniform_time_ID,'cell_methods','along_track: track_number: mean');%see 7.3.1 CF 
netcdf.putAtt(ncid,uniform_time_ID,'coverage_content_type','coordinate');

cycle_time_ID = netcdf.defVar(ncid,'cycle_time','double',dimitime);
netcdf.putAtt(ncid,cycle_time_ID,'standard_name','time');
netcdf.putAtt(ncid,cycle_time_ID,'long_name','a nonuniformly progressing time variable for each cycle defined as the mean difference, within each cycle, between the exact time at each valid observation point and the cycle-averaged time_offset');
netcdf.putAtt(ncid,cycle_time_ID,'units','days since 1950-01-01 00:00:00');
netcdf.putAtt(ncid,cycle_time_ID,'axis','T')
netcdf.putAtt(ncid,cycle_time_ID,'time_zone','UTC') ;
netcdf.putAtt(ncid,cycle_time_ID,'calendar','standard');
netcdf.putAtt(ncid,cycle_time_ID,'cell_methods','along_track: track_number: mean');%see 7.3.1 CF 
netcdf.putAtt(ncid,cycle_time_ID,'coverage_content_type','coordinate');

time_offset_ID = netcdf.defVar(ncid,'time_offset','double',[dimialong dimitrack]);
netcdf.putAtt(ncid,time_offset_ID,'standard_name','time');
netcdf.putAtt(ncid,time_offset_ID,'long_name','the mean difference, averaged over all cycles, between the time at each observation point and uniform_time');
netcdf.putAtt(ncid,time_offset_ID,'comment','For the ith along-track location, jth track, and kth cycle, time_offset(i,j)+cycle_time(k) approximates time(i,j,k) very closely, with an RMS error of 0.30 seconds and a maximum error of 3.8 seconds.');
netcdf.putAtt(ncid,time_offset_ID,'coordinates','lat lon');
netcdf.putAtt(ncid,time_offset_ID,'units','days');
netcdf.putAtt(ncid,time_offset_ID,'time_zone','UTC') ;
netcdf.putAtt(ncid,time_offset_ID,'calendar','standard');
netcdf.putAtt(ncid,time_offset_ID,'cell_methods','cycle_time: mean');%see 7.3.1 CF 
netcdf.putAtt(ncid,time_offset_ID,'coverage_content_type','referenceInformation');
netcdf.defVarFill(ncid,time_offset_ID,false,nan);

% %Having coordinate variables for non-spatial coordinates is not required, and feels cluttered
% along_track_ID = netcdf.defVar(ncid,'along_track','int',dimialong);
% netcdf.putAtt(ncid,along_track_ID,'long_name','along_track dimension');
% netcdf.putAtt(ncid,along_track_ID,'comment','data point location within a track, proceeding from south to north');
% netcdf.putAtt(ncid,along_track_ID,'coverage_content_type','coordinate');
% netcdf.putAtt(ncid,along_track_ID,'axis','Y')
% netcdf.putAtt(ncid,along_track_ID,'valid_min',1);
% netcdf.putAtt(ncid,along_track_ID,'valid_max',3373);
% 
% track_number_ID = netcdf.defVar(ncid,'track_number','int',dimitrack);
% netcdf.putAtt(ncid,track_number_ID,'long_name','track_number dimension');
% netcdf.putAtt(ncid,track_number_ID,'comment','Tracks are sorted in order of increasing longitude of their initial points, that is, of the southernmost point for ascending tracks and of the northernmost point for descending tracks.');
% netcdf.putAtt(ncid,track_number_ID,'coverage_content_type','coordinate');
% netcdf.putAtt(ncid,track_number_ID,'axis','X')
% netcdf.putAtt(ncid,track_number_ID,'valid_min',1);
% netcdf.putAtt(ncid,track_number_ID,'valid_max',254);
% 
% rej_dim_ID = netcdf.defVar(ncid,'rejected_points','int',dimirej);
% netcdf.putAtt(ncid,rej_dim_ID,'long_name','rejected_points dimension');
% netcdf.putAtt(ncid,rej_dim_ID,'comment','Suspect data points are identified with a median-based criterion, as follows.  Points whose absolute deviation from the running 21-point median is 6 or more times the running 21-point median absolute deviation, or having 10 or fewer valid data points within the 21-point window, are rejected by setting their values to NaNs within sla.  This condition is iterated three times.  The positions of the points removed in this manner is recorded in rejected_index, and the corresponding sla values in rejected_values.  Thus, sla(rejected_index)=rejected_values undoes the editing process.');
% netcdf.putAtt(ncid,rej_dim_ID,'coverage_content_type','referenceInformation');
% netcdf.putAtt(ncid,rej_dim_ID,'axis','X')
% netcdf.putAtt(ncid,rej_dim_ID,'valid_min',1);
% netcdf.putAtt(ncid,rej_dim_ID,'valid_max',102004072);

lat_ID = netcdf.defVar(ncid,'lat','double',[dimialong dimitrack]);
netcdf.putAtt(ncid,lat_ID,'standard_name','latitude');
netcdf.putAtt(ncid,lat_ID,'long_name','latitude');
netcdf.putAtt(ncid,lat_ID,'units','degrees_north');
netcdf.putAtt(ncid,lat_ID,'coverage_content_type','referenceInformation');

lon_ID = netcdf.defVar(ncid,'lon','double',[dimialong dimitrack]);
netcdf.putAtt(ncid,lon_ID,'standard_name','longitude');
netcdf.putAtt(ncid,lon_ID,'long_name','longitude');
netcdf.putAtt(ncid,lon_ID,'units','degrees_east');
netcdf.putAtt(ncid,lon_ID,'coverage_content_type','referenceInformation');

% lat_ID = netcdf.defVar(ncid,'lat','double',[dimialong dimitrack dimitime]);
% netcdf.putAtt(ncid,lat_ID,'standard_name','latitude');
% netcdf.putAtt(ncid,lat_ID,'long_name','latitude');
% netcdf.putAtt(ncid,lat_ID,'units','degrees_north');
% netcdf.putAtt(ncid,lat_ID,'coverage_content_type','referenceInformation');
% 
% lon_ID = netcdf.defVar(ncid,'lon','double',[dimialong dimitrack dimitime]);
% netcdf.putAtt(ncid,lon_ID,'standard_name','longitude');
% netcdf.putAtt(ncid,lon_ID,'long_name','longitude');
% netcdf.putAtt(ncid,lon_ID,'units','degrees_east');
% netcdf.putAtt(ncid,lon_ID,'coverage_content_type','referenceInformation');

rev_ID = netcdf.defVar(ncid,'rev','int',[dimialong dimitrack]);
netcdf.putAtt(ncid,rev_ID,'long_name','reference orbit revolution number');
netcdf.putAtt(ncid,rev_ID,'coordinates','lat lon');
netcdf.putAtt(ncid,rev_ID,'valid_min',1);
netcdf.putAtt(ncid,rev_ID,'valid_max',127);
netcdf.putAtt(ncid,rev_ID,'comment','specifies the number of the revolution within a near 10-day repeat orbit')
netcdf.putAtt(ncid,rev_ID,'coverage_content_type','referenceInformation');

atd_ID = netcdf.defVar(ncid,'atd','double',[dimialong dimitrack]);
netcdf.putAtt(ncid,atd_ID,'long_name','great circle along-track distance from southernmost point');
netcdf.putAtt(ncid,atd_ID,'coordinates','lat lon');
netcdf.putAtt(ncid,atd_ID,'units','Mm')
netcdf.putAtt(ncid,atd_ID,'coverage_content_type','referenceInformation');
netcdf.putAtt(ncid,atd_ID,'comment','1 Mm = 1000 km')

dfc_ID = netcdf.defVar(ncid,'dfc','double',[dimialong dimitrack]);
netcdf.putAtt(ncid,dfc_ID,'long_name','along-track distance from coast');
netcdf.putAtt(ncid,dfc_ID,'coordinates','lat lon');
netcdf.putAtt(ncid,dfc_ID,'units','km')
netcdf.putAtt(ncid,dfc_ID,'coverage_content_type','auxillaryInformation');

stf_ID = netcdf.defVar(ncid,'stf','int',[dimialong dimitrack]);
netcdf.putAtt(ncid,stf_ID,'long_name','surface type flag');
netcdf.putAtt(ncid,stf_ID,'comment','0 - ocean; 1 and 3 -land; 2 - inland sea or lake');
netcdf.putAtt(ncid,stf_ID,'coordinates','lat lon');
netcdf.putAtt(ncid,stf_ID,'coverage_content_type','auxiliaryInformation');

depth_ID = netcdf.defVar(ncid,'depth','double',[dimialong dimitrack]);
netcdf.putAtt(ncid,depth_ID,'standard_name','sea_floor_depth_below_geoid');
netcdf.putAtt(ncid,depth_ID,'long_name','sea floor depth');
netcdf.putAtt(ncid,depth_ID,'coordinates','lat lon');
netcdf.putAtt(ncid,depth_ID,'source','1 arc-minute bathymetry/topography grid ETOPO1 (Amante and Eakins, 2009');
netcdf.putAtt(ncid,depth_ID,'units','km')
netcdf.putAtt(ncid,depth_ID,'positive','down')
netcdf.putAtt(ncid,depth_ID,'coverage_content_type','auxillaryInformation');

flag_ID = netcdf.defVar(ncid,'flag','short',[dimialong dimitrack dimitime]);
netcdf.putAtt(ncid,flag_ID,'long_name','flag value for representing different errors in the sea surface height anomaly estimates');
netcdf.putAtt(ncid,flag_ID,'comment','0 for good data - other states are described by individual bits as in flag_meanings; states are defined only when sla is valid.');
netcdf.putAtt(ncid,flag_ID,'coordinates','cycle_time lat lon');
netcdf.putAtt(ncid,flag_ID,'flag_masks','1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384');
netcdf.putAtt(ncid,flag_ID,'flag_meanings','abs(GOT4.8-FES04)_ocean_tide>2.0cm Radiometer_Observation_is_Suspect Attitude_Out_of_Range Sigma0_Ku_Band_Out_of_Range Possible_Rain_Contamination Sea_Ice_Detected Significant_Wave_Height>8m Cross_Track_slope>10cm/km Cross_Track_Distance>1km Any_Applied_SSH_Correction_Out_of_Limits Contiguous_1Hz_Data Sigma_H_of_fit>15cm Distance_to_Land<50km Water_Depth<200m Single_Frequency_Altimeter');
netcdf.putAtt(ncid,flag_ID,'coverage_content_type','qualityInformation');

noise_ID = netcdf.defVar(ncid,'sigma','double',[dimialong dimitrack]);
netcdf.putAtt(ncid,noise_ID,'standard_name','sea_surface_height_above_reference_ellipsoid standard_error');
netcdf.putAtt(ncid,noise_ID,'long_name','small-scale noise standard deviation in sla estimated using a wavelet filter');
netcdf.putAtt(ncid,noise_ID,'units','m');
netcdf.putAtt(ncid,noise_ID,'coordinates','lat lon');
netcdf.putAtt(ncid,noise_ID,'coverage_content_type','qualityInformation');
netcdf.defVarFill(ncid,noise_ID,false,nan);

%for the mss and sla need to propagate information
struct=ncinfo(filename);
%--------------------------------------------------------------------------
mss_ID = netcdf.defVar(ncid,'mss','float',[dimialong dimitrack]);
netcdf.putAtt(ncid,mss_ID,'units','m');
netcdf.putAtt(ncid,mss_ID,'coordinates','lat lon');
netcdf.defVarFill(ncid,mss_ID,false,nan);
%found a better standard name in the table: https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
netcdf.putAtt(ncid,mss_ID,'standard_name','sea_surface_height_above_reference_ellipsoid');

for i=[2 5 8:10] %all but _FillValue, shortname, and units and coordinates
    struct.Variables(5).Name;
    att=struct.Variables(5).Attributes(i).Name;
    val=struct.Variables(5).Attributes(i).Value;
    if i==2,val=lower(val);end
    netcdf.putAtt(ncid,mss_ID,att,val);
end
%--------------------------------------------------------------------------
mdt_ID = netcdf.defVar(ncid,'mdt','float',[dimialong dimitrack]);
netcdf.putAtt(ncid,mdt_ID,'units','m');
netcdf.putAtt(ncid,mdt_ID,'coordinates','lat lon');
netcdf.defVarFill(ncid,mdt_ID,false,nan);
netcdf.putAtt(ncid,mdt_ID,'long_name','mean ocean dynamic topography');

%propagating these from mss
netcdf.putAtt(ncid,mdt_ID,'source','DTU15')
netcdf.putAtt(ncid,mdt_ID,'coverage_content_type','physicalMeasurement')
netcdf.putAtt(ncid,mdt_ID,'Reference1','Ole B. Andersen, Gaia Piccioni, L. Stenseng, P. Knudsen: The DTU15 Mean Sea Surface and Mean Dynamic Topography - focusing on Arctic issues and development. OSTST Meeting, Reston VA USA, October 19-23, 2015');
netcdf.putAtt(ncid,mdt_ID,'Reference2','Non-Ocean Component: Farr, T. G., et al. 2007, The Shuttle Radar Topography Mission, Rev. Geophys., 45, RG2004, doi:10.1029/2005RG000183.');
%--------------------------------------------------------------------------
ssh_ID = netcdf.defVar(ncid,'sla','float',[dimialong dimitrack dimitime]);
netcdf.putAtt(ncid,ssh_ID,'units','m');
netcdf.putAtt(ncid,ssh_ID,'coordinates','cycle_time lat lon');
netcdf.defVarFill(ncid,ssh_ID,false,nan);
%found a better standard name in the table: https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
netcdf.putAtt(ncid,ssh_ID,'standard_name','sea_surface_height_above_mean_sea_level');

for i=[2 4:34] %all but _FillValue and units, valid_min, and valid_max
    struct.Variables(11).Name;
    att=struct.Variables(11).Attributes(i).Name;
    val=struct.Variables(11).Attributes(i).Value;
    netcdf.putAtt(ncid,ssh_ID,att,val);
end
netcdf.putAtt(ncid,ssh_ID,'long_name','sea surface height anomaly relative to the mean sea surface given in mss');

rej_index_ID = netcdf.defVar(ncid,'rejected_index','int',dimirej);
netcdf.putAtt(ncid,rej_index_ID,'long_name','index locations for rejected data points');
netcdf.putAtt(ncid,rej_index_ID,'comment','Suspect data points are identified with a median-based criterion, as follows.  Points whose absolute deviation from a running 21-point median is 6 or more times the running 21-point median absolute deviation, or having 10 or fewer valid data points within the 21-point window, are rejected by setting their values to NaNs within sla.  This condition is iterated three times.  The positions of all points removed in this manner is recorded in rejected_index, and the corresponding sla values in rejected_values.  Thus, sla(rejected_index)=rejected_values undoes the editing process.');
netcdf.putAtt(ncid,rej_index_ID,'coverage_content_type','qualityInformation');

rej_values_ID = netcdf.defVar(ncid,'rejected_values','double',dimirej);
netcdf.putAtt(ncid,rej_values_ID,'standard_name','sea_surface_height_above_mean_sea_level');
netcdf.putAtt(ncid,rej_values_ID,'units','m');
netcdf.putAtt(ncid,rej_values_ID,'long_name','sea surface height anomaly values of rejected data points');
netcdf.putAtt(ncid,rej_values_ID,'comment','see comment under rejected_index');
netcdf.putAtt(ncid,rej_values_ID,'coverage_content_type','qualityInformation');

%found a better short name in the table: https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
%--------------------------------------------------------------------------
netcdf.endDef(ncid);
%--------------------------------------------------------------------------
netcdf.putVar(ncid,uniform_time_ID,num_uniform);
netcdf.putVar(ncid,cycle_time_ID,num);
%netcdf.putVar(ncid,along_track_ID,[1:3373]');
%netcdf.putVar(ncid,track_number_ID,[1:254]');
netcdf.putVar(ncid,time_offset_ID,dnum);
netcdf.putVar(ncid,time_ID,numall);
%netcdf.putVar(ncid,lat_ID,vrep(lat,length(num),3));
%netcdf.putVar(ncid,lon_ID,vrep(lon,length(num),3));
netcdf.putVar(ncid,lat_ID,lat);
netcdf.putVar(ncid,lon_ID,lon);
netcdf.putVar(ncid,rev_ID,rev);
netcdf.putVar(ncid,atd_ID,atd);
netcdf.putVar(ncid,dfc_ID,dfc);
netcdf.putVar(ncid,stf_ID,stf);
netcdf.putVar(ncid,depth_ID,-topo);
netcdf.putVar(ncid,flag_ID,flag);
netcdf.putVar(ncid,noise_ID,sigma);
netcdf.putVar(ncid,rej_index_ID,rejected_index);
netcdf.putVar(ncid,rej_values_ID,rejected_values);
netcdf.putVar(ncid,mss_ID,mss);
netcdf.putVar(ncid,mdt_ID,mdt);
netcdf.putVar(ncid,ssh_ID,ssh);
netcdf.close(ncid)
%\*************************************************************************




%try loading ssh

% tic
% sla=ncread([writedir name_nc],'sla');
% toc
%netcdf.putAtt(ncid,depth_ID,'long_name','sea floor depth');

%useful discussion from the CF conventions
%https://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#terminology
%  
% auxiliary coordinate variable
% Any netCDF variable that contains coordinate data, but is not a coordinate 
% variable (in the sense of that term defined by the NUG and used by this 
% standard - see below). Unlike coordinate variables, there is no 
% relationship between the name of an auxiliary coordinate variable and the 
% name(s) of its dimension(s).
% 
% coordinate variable
% We use this term precisely as it is defined in the NUG section on 
% coordinate variables. It is a one- dimensional variable with the same name
% as its dimension [e.g., time(time) ], and it is defined as a numeric data
% type with values that are ordered monotonically. Missing values are not
% allowed in coordinate variables.
%
% There are two methods used to identify variables that contain coordinate 
% data. The first is to use the NUG-defined "coordinate variables." The use
% of coordinate variables is required for all dimensions that correspond to
% one dimensional space or time coordinates . In cases where coordinate 
% variables are not applicable, the variables containing coordinate data are 
% identified by the coordinates attribute.
%
% What I gauge from this reference 
%   https://ngdc.noaa.gov/wiki/index.php/ISO_19115_and_19115-2_CodeList_Dictionaries#MD_CoverageContentTypeCode
% is that coverage_content_type is supposed to apply to coordiante
% variables, not auxiliary coordinate variables.  Further, the example for 
% referenceInformation says:
% 
% reference information used to support the calculation or use of the 
% physicalMeasurement coverages in the dataset (e.g. grids of 
% latitude/longitude used to geolocate the physical measurements).
%
% based on this, I think that lat/lon were misclassified in the PO.DAAC
% version as "coordinates" but that they are actually referenceInformation
%
% also, coordinates are not allowed to have missing variables

%deflate seems to slow things down a lot
%netcdf.defVarDeflate(ncid,ssh_ID,true,true,1);%shuffle and deflate level 1
%netcdf.defVarDeflate(ncid,ssh_ID,true,true,2);%shuffle and deflate level 2

%recommended deflate level is 5 according to 
%http://climate-cms.wikis.unsw.edu.au/NetCDF_Compression_Tools
%https://www.unidata.ucar.edu/blogs/developer/entry/netcdf_compression
%--------------------------------------------------------------------------

 
% %/********************************************************
% %dx=1000*5.75;
% %dx=1000*abs(vdiff(atd,1));   %Into meters
% %c=frac(1,100)*9.81./(abs(corfreq(lat))/3600)./dx;
% %vswap(c,nan,0);
% %veke=c.*vdiff(ssh,1);
% %veke=c.*diff([0*ssh(1,:);ssh]+[ssh;0*ssh(1,:)],1,1)/2;
% 
% %description='Integrated alongtrack altimetry dataset from Brian Beckley';
% %link='http://podaac.jpl.nasa.gov/dataset/MERGED_TP_J1_OSTM_OST_CYCLES_V2';
% %creator=[jlab_version ', function L_BECKLEY'];
% %timestamp=datestr(now);
% %\**********************************************************
% 
% %How does my definition of depth relate to that from Beckley? 
% [topo2,lattopo,lontopo]=readtopo([-180 180 -80 80]);
% dlon=lontopo(2)-lontopo(1);
% 
% lontopo=[lontopo(1)-2*dlon lontopo(1)-dlon lontopo lontopo(end)+dlon lontopo(end)+2*dlon];
% topo2=[topo2(:,end-1) topo2(:,end) topo2 topo2(:,1) topo2(:,2)];
% 
% depth2=-interp2(lontopo,lattopo,topo2,lon,lat);
% 
% ddepth=depth2+topo;
% 
% %how many "ocean" data points do I find with negative depth?
% length(find(stf==0&depth2<0))  %742, not many
% 
% %how many "land" data points do I find with positive depth?
% length(find((stf==1|stf==3)&depth2>0))  %2652, not many
% 
% [jj,ii]=find(stf==0&depth2<0);
% jpcolor(depth2+topo),caxis([-1 1]),hold on,plot(ii,jj,'r*')
% %only along the coasts, as one might expect
% 
% %do these points tend to look particularly bad?
% sshstd=vstd(ssh,3);
% 
% %\**********************************************************
% 
% %plot(lon(:,1:127),lat(:,1:127),'b.'),hold on,plot(lon(:,128:end),lat(:,128:end),'r.')
% if 0
% %Set to NaNs over land
% l1=length(find(isnan(ssh)));
% for i=1:size(ssh,3)
%     sshi=ssh(:,:,i);
%     sshi(topo==0)=nan;
%     ssh(:,:,i)=sshi;
% end
% l2=length(find(isnan(ssh)));
% disp([int2str(l2-l1) ' data points changed to NaNs over land.'])
% 
% ngood=vsum(~isnan(ssh)+0,3);
% %figure,jpcolor(ngood)
% 
% % %Set to NaNs where there is data missing 3/4 of the time or more (mostly ice)
% % index=find(ngood./size(ssh,3)<1/4);
% % for i=1:size(ssh,3)
% %     sshi=ssh(:,:,i);
% %     sshi(index)=nan;
% %     ssh(:,:,i)=sshi;
% % end
% l3=length(find(isnan(ssh)));
% %disp([int2str(l3-l2) ' data points changed to NaNs where more than 3/4 of the data is missing.'])
% %disp([num2str(100*(l3-l2)./length(find(~isnan(ssh)))) ' percent of the data affected.'])
% %disp([int2str(length(find(~isnan(ssh)))) ' good data points remain.'])
% 
% %Despiking with 3-point difference statistics
% %This needs to happen on heimdall
% %/*************************************************************************
% end
