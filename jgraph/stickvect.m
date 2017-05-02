function[h]=stickvect(arg1,arg2,arg3,arg4,arg5,arg6)
%STICKVECT Plots "stick vectors" for multicomponent velocity time series.
%  
%   STICKVECT(TIME,SCALE,U,V), where U and V are measurements of
%   "eastward" and "northward" velocities respectively, represents the
%   velocity points as vectors strung out along the x-axis.
%
%   The vector TIME is used as the x-axis, and SCALE is the stretching
%   factor used to map U values onto TIME.  SCALE velocity units
%   correspond to one X unit. The aspect ratio of the plot is set so
%   that angles are not distorted: the point (U,V)=(1,1) will make a
%   45 degree angle.
%
%   STICKVECT(TIME,SCALE,U,V,NSTICKS) plots only NSTICKS vectors
%   evenly spaced throughout the length of the input vectors.
%
%   STICKVECT may also be used with multicomponent time series. In
%   this case U and V are matrices the columns of which are distinct
%   time series from (for example) different depths.
% 
%   STICKVECT(TIME,SCALE,U,V,NSTICKS,SHIFT), when U and V are
%   matrices, will offset the stickvector lines in the vertical by
%   amount SHIFT. SHIFT may either be a scalar, denoting the offset
%   between adjacent time series, a vector containing the "baseline"
%   values for each instrument, or a matrix of same size as U and V.
%   TIME may be a vector, or it may be a matrix of the same size as U
%   and V if the time-axes are not all the same.
%	
%   STICKVECT(TIME,CV,...), where CV=U+iV, also works.
%
%   STICKVECT plots a small diagonal line in the upper right-hand
%   corner of the plot, which should be at a 45 degree angle; this
%   provides a check that the aspect ratio is correct.  
%  
%   H=STICKVECT(...) returns a vector containing handles to the
%   stickvector lines of the various time series.
%
%   To adjust the appearance of the plot, first vary SCALE and then
%   vary NSTICKS.
%	
%   As an example,  
%  
%         load bravo94
%         cv=vfilt(bravo94.rcm.cv,100);
%         stickvect(yearfrac(bravo94.rcm.num),180,cv,300,-30);
% 
%   recreates the upper part of Fig. 24 of Lilly et. al JPO 1999. 
%   This may be run by typing 'stickvect --f'.
%		
%   Usage:  stickvect(time,scale,u,v);
%           stickvect(time,scale,u,v,nsticks);
%           stickvect(time,scale,u,v,nsticks,shift);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 1998--2015 J.M. Lilly --- type 'help jlab_license' for details      

if strcmpi(arg1, '--f')
  type makefigs_stickvect
  makefigs_stickvect;
  return
end

na=nargin;

x=arg1(:);
eval(argadvance(nargin))
na=na-1;

if length(arg1)==1
	factor=arg1;
	eval(argadvance(nargin))
	na=na-1;
else
	factor=1;
end


if isreal(arg1)
	jcv=arg1+sqrt(-1)*arg2; 
	eval(argadvance(nargin))
	na=na-1;
else
	jcv=arg1;
end

if isempty(x) || aresame(x,0)
	x=0*jcv;
end
if size(x,2)<size(jcv,2)
	x=x(:,1)*ones(size(jcv(1,:)));
end

if na>2,
	if length(arg3)==size(jcv,2)
		shift=conj(arg3(:)');
	else
		shift=arg3;
	end
else	
	shift=0;
end
if length(shift)==1
	shift=shift*(0:1:size(jcv,2)-1);	
end

if size(shift,1)~=size(jcv,1)
	shift=ones(size(jcv(:,1)))*shift;
end
if na>1
	nstick=arg2;
	dl=floor(length(jcv)/nstick);
	index=(1:dl:nstick*dl)';
	[jcv,shift,x]=vindex(jcv,shift,x,index,1);
end



%/***********************************
%weave together the endpoints of the vectors
weave=zeros(size(jcv,1)*3,size(jcv,2));
index=(1+(0:length(jcv)-1)*3)';
weave(index,:)=x;
weave(index+1,:)=real(jcv)./factor+sqrt(-1)*imag(jcv)+x;
weave(index+2,:)=nan; 

shiftmat=zeros(size(weave));
for i=1:3
	if isreal(shift)
	  shiftmat(index+i-1,:)=shift*sqrt(-1);
	else
	  shiftmat(index+i-1,:)=shift;
	end
	
end
sticks=weave+shiftmat;
%\************************************


h=plot(sticks,'k');
set(gca,'dataaspectratio',[1./factor 1 1]);
set(gca,'dataaspectratiomode','manual');


%/*****************************
%plot a little line in the corner so we know if it's true or not
hold on
jcv=real(jcv)./factor+sqrt(-1)*imag(jcv);
len=max(max(imag(jcv)))/10;
ax=axis;
top=ax(4)*sqrt(-1);
right=ax(2);
plot([top+right-len./factor,top+right-len*sqrt(-1)],'k')
%\*******************************



function[evalme]=argadvance(na,n)

if nargin==1
	n=2;
end
evalme=[];
for j=n:na
	evalme=[evalme 'arg',int2str(j-1),'=arg',int2str(j),';'];
end

