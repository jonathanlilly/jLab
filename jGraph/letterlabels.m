function[h]=letterlabels(arg1,arg2,arg3)
%LETTERLABELS	For automatically putting letter labels on subplots.
%
%   LETTERLABELS puts letter labels '(a)','(b)','(c)' etc., on the
%   subplots of the current figure, in the upper left-hand corner of
%   each subplot window.
%
%   LETTERLABELS(I) specifies the labelling position, clockwise from 
%   top left: 
%
%		I=1	Top left
%		I=2	Top right
%		I=3	Bottom right
%		I=4	Bottom left
%
%   LETTERLABELS(H,I), where H is a vector of handles to subplots,
%   puts labels on the subplots indicated by H.
%
%   LETTERLABELS(H,I,'d'), begins the labelling with letter 'd'.
%
%   LETTERLABELS will not put the letters on the right subplots if you
%   click on them before hand (because this changes the order of the
%   subplot handles).  Also if you want to label the subplots in a
%   nonstandard order, do this by reordering H.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2016 J.M. Lilly --- type 'help jlab_license' for details


if strcmpi(arg1,'--t')
    return
end

py=0.08;
ar=get(gca,'plotboxaspectratio');
xstretch=ar(1)./ar(2);
px=py/xstretch*1.8;

firstlet=real('a');
i=1;
axhand=flipud(get(gcf,'children'));

for j=1:nargin
    xx=eval(['arg' int2str(j)]);
    isaxhandle=false;
    if ishandle(xx)
        if strcmp(get(xx,'Type'),'axes')
            isaxhandle=true;
        end
    end
    if ischar(xx)&~isaxhandle
        firstlet=real(xx);
    elseif (length(xx)==1)&&~isaxhandle
        i=xx;
    else
        axhand=xx;
    end
end
nax=length(axhand);
fact=1/100;

for j=1:length(axhand)
    axes(axhand(j))
    ax=axis;
    isrevx=strcmpi(get(gca,'xdir'),'reverse');
    isrevy=strcmpi(get(gca,'ydir'),'reverse');
    
    if ~isrevx
        x1=ax(1);
        x2=ax(2);
    else
        x1=ax(2);
        x2=ax(1);
    end
    if ~isrevy
        y1=ax(3);
        y2=ax(4);
    else
        y1=ax(4);
        y2=ax(3);
    end
    
    if i==1
        t=y2-(y2-y1)*fact;
        b=y2-(y2-y1)*py;
        l=x1+(x2-x1)*fact;
        r=x1+(x2-x1)*px;
    elseif i==2
        t=y2-(y2-y1)*fact;
        b=y2-(y2-y1)*py;
        l=x2-(x2-x1)*px;
        r=x2-(x2-x1)*fact;
    elseif i==3
        t=y1+(y2-y1)*py;
        b=y1+(y2-y1)*fact;
        l=x2-(x2-x1)*px;
        r=x2-(x2-x1)*fact;
    elseif i==4
        t=y1+(y2-y1)*py;
        b=y1+(y2-y1)*fact;
        l=x1+(x2-x1)*fact;
        r=x1+(x2-x1)*px;
    end
    
    %Tried setting log to lin before adding label but didn't help
    if ~strcmpi(get(gca,'yscale'),'log')
        h(j)=patch([l l r r],[b t t b],'w');
        set(h(j),'edgecolor',[1 1 1])
    end
    cx=l/2+r/2;
    cy=t/2+b/2;
    text(l+(r-l)*0.2,cy,['(' char(firstlet-1+j),')']);
    axis(ax)
end

if nargout==0
    clear h
end








