function[h]=packfig(varargin)
%PACKFIG  Squeeze together rows and/or columns of the current figure.
%  
%   PACKFIG(M,N) squeezes together all rows and columns in the current 
%   figure, which has M rows and N columns generated with SUBPLOT. This is
%   used when all subplots in a given column share common x- and y-axes.
%  
%   PACKFIG(M,N,'columns') only packs the columns, and not the rows.
%   PACKFIG(M,N,'rows') only packs the rows, and not the columns.
%
%   When the columns are packed, the Y-axis tick labels for the interior 
%   subplots are removed.  Similarly, when the rows are packed, the X-axis
%   tick labels for the interior subplots are removed.
%
%   H=PACKFIG(M,N) returns a vector of handles H into the subplots.
%
%   PACKFIG(H,M,N) also works, where H is a vector of handles.
%
%   After calling PACKFIG, do not call SUBPLOT again, as this will destroy
%   the subplots; instead, access the subplots through AXES(H(I)) where I 
%   is an integer. 
%   _________________________________________________________________
%
%   Use with AXIS EQUAL
%
%   PACKFIG may not work when AXIS EQUAL and AXIS([A B C D]) are used 
%   together because of the way Matlab constricts the figure with AXIS 
%   EQUAL. Consider altering the aspect ratio with 
%
%       set(gcf,'paperposition',[left bottom width height]) 
%
%   or
%
%       set(gcf,'position',[left bottom width height]) 
%
%   in order to obtain the desired result.  
%   [Thanks to JLAB user Rodrigo Duran for this note.]
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details   

if isstr(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
else
    str='both';
end

if length(varargin)==2
  m=varargin{1};
  n=varargin{2};
  h=subplots(m,n);
elseif length(varargin)==3
  h=varargin{1};
  m=varargin{2};
  n=varargin{3};
end

if ~strcmpi(str(1:3),'col')
    packrows(h,m,n);
end
if ~strcmpi(str(1:3),'row')
    packcols(h,m,n);
end
if nargout==0
  clear h
end


function[h]=packrows(h,m,n)
%PACKROWS   Squeeze together all subplot rows of the current figure.

dd=.01;

if m>1
  for i=1:n
    for j=1:m-1
      index=(j-1)*n+i;
      axes(h(index))
      noxlabels,xlabel('')
    end
  end

  axes(h(1))
  pos1=get(gca,'position');
     
  axes(h(end))
  pos2=get(gca,'position');
      
  dy=pos1(2)+pos1(4)-pos2(2);    %total height to work with 
  
  yheight=frac(dy-dd*(m-1),m);  %This is how tall each axis should be 
  ybottom=zeros(m,1);
  ybottom(1)=pos2(2);
  
  for j=2:m;
      ybottom(j)=ybottom(j-1)+yheight+dd;
  end
  ybottom=flipud(ybottom);
  
  
  for i=1:n
    for j=1:m
      index=(j-1)*n+i;
     
      axes(h(index))
      pos=get(gca,'position');
      isboxoff=strcmpi(get(gca,'box'),'off');
      
      pos(2)=ybottom(j);
      pos(4)=yheight;
     
      set(gca,'position',pos,'box','on')
      if isboxoff
          boxoff
      end
    end
  end
end

if nargout==0
  clear h
end


function[h]=packcols(h,m,n)
%PACKCOLS   Squeeze together all subplot columns of the current figure.
dd=.01;

if n>1  
  for i=2:n
    for j=1:m
      index=(j-1)*n+i;
      axes(h(index))
      noylabels,ylabel('')
    end
  end

  axes(h(1))
  pos1=get(gca,'position');
     
  axes(h(end))
  pos2=get(gca,'position');
      
  dx=pos2(1)+pos2(3)-pos1(1);   %total width to work with 
  
  xwidth=frac(dx-dd*(n-1),n);  %This is how wide each axis should be 
  xleft=zeros(n,1);
  xleft(1)=pos1(1);
  
  for i=2:n;
      xleft(i)=xleft(i-1)+xwidth+dd;
  end
  %xleft=flipud(xleft);
  
  
  for i=1:n
    for j=1:m
      index=(j-1)*n+i;
     
      axes(h(index))
      pos=get(gca,'position');
      isboxoff=strcmpi(get(gca,'box'),'off');
      
      pos(1)=xleft(i);
      pos(3)=xwidth;
     
      set(gca,'position',pos,'box','on')
      if isboxoff
          boxoff
      end
    end
  end
end


if nargout==0
  clear h
end


function[h]=subplots(m,n)
%SUBPLOTS  Return handles to subplots of the current figure.
%
%   H=SUBPLOTS(M,N) returns a vector H containing all handles to the
%   subplots of the current figure, where M is the number of rows and
%   N is the number of columns.
%
%   AXES(H(i)) is then equivalent to SUBPLOT(M,N,i).  
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000, 2004 J.M. Lilly --- type 'help jlab_license' for details      

  
for i=1:n
  for j=1:m
     index=(j-1)*n+i;
     h(index)=subplot(m,n,index);
  end
end  

