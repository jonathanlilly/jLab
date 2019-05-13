function[dsn,xsn,ysn,fullindexn]=polysmooth_cat(sizez,n,K,ds,xs,ys,index)
%POLYSMOOTH_CAT
%
%   POLYSMOOTH_CAT is an auxiliary function for use with POLYSMOOTH.
%
%   [DS,XS,YS,TS,ZS,WS]=POLYSMOOTH_CAT(N,K,DS,XS,YS,INDEX,T,Z,W);
%
%   Usage: [ds,xs,ys,ts,zs,ws]=polysmooth_cat(n,K,ds,xs,ys,index,t,z,w);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2018 J.M. Lilly --- type 'help jlab_license' for details
 
%sort first so we don't have to sort later
[ds,xs,ys,~,index]=polysmooth_presort(ds,xs,ys,[],index);

[ii,jj]=ind2sub(sizez,index);

[fullii,fulljj,fullkk]=vzeros(size(ds,1),size(ds,2),size(ds,3),K);

for k=-(K-1)/2:(K-1)/2
    %k+(K+1)/2
    %size(ii),size(fullii(:,:,k+(K+1)/2))
    fullii(:,:,:,k+(K+1)/2)=ii;
    fulljj(:,:,:,k+(K+1)/2)=jj;
    nplusk=min(max(1,n+k),sizez(3));  %To avoid going out of bounds   
    fullkk(:,:,:,k+(K+1)/2)=nplusk+0*jj;
end
fullindex=sub2ind(sizez,fullii,fulljj,fullkk);
clear fullii fulljj fullkk ii jj

[dsn,xsn,ysn,fullindexn]=vzeros([size(ds,1),size(ds,2),K*size(ds,3)],nan);

%interleave so we don't need to sort
newindex=[1:size(ds,3)]';
for k=1:K
    dsn(:,:,(newindex-1)*K+k)=ds;
    xsn(:,:,(newindex-1)*K+k)=xs;
    ysn(:,:,(newindex-1)*K+k)=ys;
    fullindexn(:,:,(newindex-1)*K+k)=fullindex(:,:,:,k);
end

%dsn=vrep(ds,K,3);
%xsn=vrep(xs,K,3);
%ysn=vrep(ys,K,3);
%fullindexn=reshape(fullindex,[size(ds,1),size(ds,2),size(ds,3)]);

%vsize(ds,xs,ys,fullindex)

%must sort!!! 
%[ds,xs,ys,~,fullindex]=polysmooth_presort(ds,xs,ys,[],fullindex);
%note I can put fullindex into the 'time' argument slot in calling presort
%to correctly sort it (but not into the 'z' slot)
%     
%     if ~isempty(t)
%         tk=t(:,:,n+k);
%     end
%     if ~isempty(z)
%        zk=z(:,:,n+k);
%     end
%     if ~isempty(w)
%        wk=z(:,:,n+k);
%     end
%     
%     [tsk,zsk,wsk]=polysmooth_index(size(ds),index,tk,zk,wk);
%     if ~isempty(t)
%         tk(:,:,k+(K+1)/2)=tsk;
%     end
%     if ~isempty(z)
%         zk(:,:,k+(K+1)/2)=zsk;
%     end
%     if ~isempty(w)
%         wk(:,:,k+(K+1)/2)=wsk;
%     end
% end
% 
% 
% 
% [ts,zs,ws]=vzeros(size(ds,1),size(ds,2),size(ds,3),K);
% 
% [tk,zk,wk]=vempty;
% 
% for k=-(K-1)/2:(K-1)/2
%     if ~isempty(t)
%         tk=t(:,:,n+k);
%     end
%     if ~isempty(z)
%        zk=z(:,:,n+k);
%     end
%     if ~isempty(w)
%        wk=z(:,:,n+k);
%     end
%     
%     [tsk,zsk,wsk]=polysmooth_index(size(ds),index,tk,zk,wk);
%     if ~isempty(t)
%         tk(:,:,k+(K+1)/2)=tsk;
%     end
%     if ~isempty(z)
%         zk(:,:,k+(K+1)/2)=zsk;
%     end
%     if ~isempty(w)
%         wk(:,:,k+(K+1)/2)=wsk;
%     end
% end
% 
% 
% ds=vrep(ds,K,3);
% xs=vrep(xs,K,3);
% ys=vrep(ys,K,3);
% 
% if ~isempty(t)
%     tk=reshape(tk,[size(tk,1),size(tk,2),size(tk,3)*size(tk,4)]);
% end
% if ~isempty(z)
%     zk=reshape(zk,[size(zk,1),size(zk,2),size(zk,3)*size(zk,4)]);
% end
% if ~isempty(w)
%     wk=reshape(wk,[size(wk,1),size(wk,2),size(wk,3)*size(wk,4)]);
% end
% 
