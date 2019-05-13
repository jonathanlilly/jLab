function[B]=polysmooth_bandwidth(P,d,w)
%POLYSMOOTH_BANDWIDTH  Determine bandwidth given population for POLYSMOOTH.
%
%   POLYSMOOTH_BANDWIDTH is a low-level function called by POLYSMOOTH. 
%
%   B=POLYSMOOTH_BANDWIDTH(P,DS,WS) returns the spatially-varying bandwidth 
%   B implied by the population P. DS and WS are fields as output from 
%   TWODSORT and SPHERESORT.  DS is distance and WS is an optional weight.
%
%   P may be a scalar or matrix of size M x N, where M and N are the sizes 
%   of the first two dimensions of DS and WS.  B will always be M x N.
%
%   See also POLYSMOOTH.
%
%   Usage: B=polysmooth_bandwidth(P,ds,ws);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2018 J.M. Lilly --- type 'help jlab_license' for details
 

%Convert number of points input to equivalent bandwidth
P=ceil(P);%What I've actually input is the number of points we want to keep

if minmin(P)<=1
    error('With fixed population algorithm, population must be greater than one.')
end
%If a matrix-valued number of points was entered, need to replicate
%the matrix along the third dimension to make the rest of the code work

%figure,jpcolor(d(:,:,900)),caxis
if isempty(w)  %Set weight to zero if additional weight factor vanishes
    w=ones(size(d));
else
    w=0+(abs(w)>0);
end
%size(d)
if ~ismatrix(d)
    if aresame(size(P),size(d(:,:,1)))
        P=vrep(P,size(d,3),3);
    end
end
w(~isfinite(d))=0; %Important line, makes us handle missing data correctly
%Compute number of data points with potential for weighted points

%vsize(P,d,w)

if  ismatrix(d) %false if ndims>=3
    B=nan;
    cumw=cumsum(double(w),1);
    %vsize(cumw,w,P)
    index=find(cumw>=P,1,'first');
    if ~isempty(index)
        B=squeeze(d(index));
    end
else  %Only do this for the 'fast' algorithm
    cumw=cumsum(double(w),3);
    %figure,plot(squeeze(cumw)')  
    %vsize(d,B,cumw,P)
    %figure,jpcolor(max(cumw./P,[],3))
    %d(cumw<P)=nan;
    %B=min(d,[],3);  %Loopless way to find bandwidth
    d(cumw>P)=nan;
%    size(d)
    B=max(d,[],3);  %Loopless way to find bandwidth
    %This is the same as the above
%    size(P)
%    figure,jpcolor(squeeze(cumw))
%     figure,jpcolor(squeeze(d))
%     for i=1:size(d,1)
%         for j=1:size(d,2)
%             index=find(cumw(i,j,:)>=P,1,'first');
%             
%             if ~isempty(index)
%                 B(i,j)=squeeze(d(i,j,index));
%             else
%                 B(i,j)=squeeze(d(i,j,end));
%             end
%         end
%     end
%    aresame(B,B1)
        
        %figure,plot(squeeze(P))
%       figure,plot(squeeze(B1)),hold on,plot(squeeze(B),'r')
    end


