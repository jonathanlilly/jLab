function[id,ii,jj,xr,fr]=ridgechains(fs,N,bool,x,f,alpha,mask,rq)
%RIDGECHAINS  Forms ridge curves by connecting transform ridge points.
%
%   RIDGECHAINS is a low-level function called by RIDGEWALK.
%
%   [ID,IR,JR,WR]=RIDGECHAINS(N,BOOL,W) forms chains of ridge points of
%   wavelet transform W.  
%
%   Ridge points are all points of W for which BOOL, a matrix of the same
%   size as W, is true.  Only ridges with at least N periods are returned.
%   
%   ID is a unique ID number assigned to each ridge.  IR and JR are the 
%   time- and scale-indices along the ridges.  WR is the wavelet transform 
%   along the ridge. 
% 
%   All output variables are the same size.
% 
%   See also RIDGEWALK.
%
%   Usage:  [id,ii,jj,xr]=ridgechains(N,bool,x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2018 J.M. Lilly --- type 'help jlab_license' for details

if isempty(find(bool,1))
    [id,ii,jj,xr,fr]=vempty;return
end

indexridge=find(bool);
[ii,jj]=ind2sub(size(bool),indexridge);
[ii,sorter]=sort(ii);
vindex(jj,indexridge,sorter,1);

%New version as of April 2018, using quadratic interpolation
%/*************************************************************************
dfdt=vdiff(f,1);
%Interpolate frequency within ridge for better estimation
[fr,dfdt,xr]=ridgeinterp(rq,ii,jj,f,dfdt,x);

%Old version
%dfdt=vdiff(f,1);
%dfdt=dfdt(indexridge);
%xr=x(indexridge);             %Transform value along ridge
%fr=f(indexridge);             %Frequency along ridge

%New version actually makes a significant difference. 
%Probably will help to minimize ridge breaking
%fro=f(indexridge);figure,plot(fr,'.'),hold on,plot(fro,'o')
%Try uncommenting this on the example figure

fsr=fs(jj);       %Scale frequency along ridge
fr_next=fr+dfdt;  %Predicted frequency at next point
fr_prev=fr-dfdt;  %Predicted frequency at previous point
%\*************************************************************************


cumbool=cumsum(bool,2);
J=maxmax(cumbool(:,end));   %Maximum number of ridges at any one time

[indexmat,nextindexmat,iimat,jjmat,fsmat,frmat,fr_nextmat,fr_prevmat]=vzeros(size(f,1),J,'nan');

%Indices for this point
indexmat(sub2ind(size(indexmat),ii,cumbool(indexridge)))=(1:length(ii));
nonnanindex=find(~isnan(indexmat));

%Don't overwrite the original variables
[ii1,jj1,fsr1,fr1,fr_next1,fr_prev1]=vindex(ii,jj,fsr,fr,fr_next,fr_prev,indexmat(nonnanindex),1);
vindexinto(iimat,jjmat,fsmat,frmat,fr_nextmat,fr_prevmat,ii1,jj1,fsr1,fr1,fr_next1,fr_prev1,nonnanindex,0);
clear ii1 jj1 fsr1 fr1 fr_next1 fr_prev1

clear iimat

%Scale frequency difference from points here to next points
fsmat3=vrep(fsmat,J,3);
frmat3=vrep(frmat,J,3);

%Predicted minus actual frequency at this point
fr_nextmat3=vrep(fr_nextmat,J,3);
df1=permute(circshift(frmat3,-1,1),[1 3 2])-fr_nextmat3;
df1=frac(df1,frmat3);

clear fr_nextmat3

%Expected minus actual frequency at next point
fr_prevmat3=vrep(fr_prevmat,J,3);
df2=permute(circshift(fr_prevmat3,-1,1),[1 3 2])-frmat3;
df2=frac(df2,frmat3);

df=frac(abs(df1)+abs(df2),2);

%Note: I had stopped using alpha, but it actually does help 
df(df>alpha)=nan;

clear fr_prevmat3 fr_mat3 df1 df2 

%Keep when they are the same 
%Set df to nan except when min along one direction
[mindf,jjmin]=min(df,[],2);
iimin=vrep((1:size(df,1))',size(df,2),2);
kkmin=vrep(1:size(df,2),size(df,1),1);
df=nan*df;
df(sub2ind(size(df),iimin,squeeze(jjmin),kkmin))=mindf;
clear iimin jjmin kkmin 

[mindf,jjmin]=min(df,[],3);

index=find(~isnan(mindf));

if ~isempty(index)
    mindf=mindf(index);
    [ii2,jj2]=ind2sub(size(indexmat),index);
    [ii2,~,index,~]=vindex(ii2,jj2,index,mindf,find(ii2<size(x,1)),1);
    index2=sub2ind(size(indexmat),ii2+1,jjmin(index));
    nextindexmat(index)=indexmat(index2);
end

id=nan*ii;
%Assign a unique number to all ridge points
for i=1:size(nextindexmat,1)
    id(nonnan(indexmat(i,:)))=nonnan(indexmat(i,:));
end

%Reassign number for linked points
for i=1:size(nextindexmat,1)
    id(nonnan(nextindexmat(i,:)))=id(indexmat(i,~isnan(nextindexmat(i,:))));
end
[id,sorter]=sort(id);
vindex(ii,jj,indexridge,xr,fr,sorter,1);
    
%Remove ridge lines of length shorter than a specified length
lr=ridgelen(1,id,fr);
vindex(id,ii,jj,indexridge,xr,fr,lr,find(lr>N),1);

%Remove isolated points.
bool=~(~isnan(id)&isnan(vshift(id,-1,1))&isnan(vshift(id,1,1)));
vindex(id,ii,jj,indexridge,xr,fr,lr,bool,1);

%This block is for applying the mask
%/**************************************************************
%Keep from chaining ridges through mask
%isempty(mask)
if ~isempty(mask)&&~isempty(ii)
    [~,~,ib]=blocklen(id);  
    jjnext=circshift(jj,-1,1);
    jjnext(ib)=jj(ib);
    
    bool=false(size(ii));
    
    for i=1:length(ii)
        jjindex=min(jj(i),jjnext(i)):max(jj(i),jjnext(i));
        bool(i)=all(mask(ii(i),jjindex));
    end
    %figure,plot(bool)
    vindex(id,ii,jj,indexridge,xr,fr,lr,find(bool),1); 
    
    %Renumber
    breaks=false(size(id));
    breaks(find(diff(ii)>1)+1)=true;
    breaks(find(diff(id)>0)+1)=true;
    id=cumsum(breaks);
  
    %And remove short ridges again
    lr=ridgelen(1,id,fr);
    vindex(id,ii,jj,indexridge,xr,fr,lr,find(lr>N),1);
    
    %Remove isolated points. 
    bool=~(~isnan(id)&isnan(vshift(id,-1,1))&isnan(vshift(id,1,1)));
    vindex(id,ii,jj,indexridge,xr,fr,lr,bool,1);
    %figure,plot(ii,lr,'.'),hlines((2*sqrt(6)/pi)),ylog,N
end
%\**************************************************************

disp(['RIDGEWALK pruning to ' int2str(length(id)) ' ridge points.'])
