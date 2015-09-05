function[h]=cellplot(varargin)
%CELLPLOT  Rapidly plot all elements of a cell array of numeric arrays.
%
%   CELLPLOT(X) where X is a cell array containing N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   is the same as 
%
%      PLOT(X1), PLOT(X2),..., PLOT(XN).
%
%   However, CELLPLOT is vastly quicker than the obvious loop and also has
%   some handy features, as described below.
%
%   CELLPLOT(X,Y), where X and Y are both cell arrays of N arrays, performs
%
%      PLOT(X1,Y1), PLOT(X1,Y1),..., PLOT(XN,YN).
%
%   H=CELLPLOT(...) returns the handle to the plots.  
%
%   CELLPLOT is compatible with the M_MAP toolbox, as described below.%
%
%   CELLPLOT(X) or CELLPLOT(X,Y) still works if X and Y are not cell 
%   arrays, but are vectors or matrices.  This is useful in accessing the 
%   additional options described below. 
%   __________________________________________________________________
%   
%   Additional options
%
%   CELLPLOT(...,STY) uses the linestyle specified by STY.  STY is a
%   string following the format in LINESTYLE, e.g. STY='2b g r--'. 
%
%   CELLPLOT(...,'M_MAP') will work with Rich Pawlowicz's M_MAP 
%   package by calling M_PLOT.  
%
%   CELLPLOT(...,INDEX) only plots the elements of the cell array 
%   indicated by index.
%
%   The string arguments and INDEX can be combined provided they are 
%   after the one or two input cell array variables.
%   __________________________________________________________________
%
%   CELLPLOT on the sphere
%    
%   When plotting data on the sphere, there is an annoying wrap-around
%   effect when the data crosses the dateline at longitude 180.  This
%   leads to lines extending all the way across the plot.
%
%   CELLPLOT(LONO,LON,LAT) where LONO is the value of a cutoff longitude
%   uses LONO as the right-hand-edge of the plot, and LONO-360 as the 
%   left-hand-edge, with no wraparound effects.  
%
%   This works together also with the 'm_map' option.  In this case LONO 
%   should be between -180 and 180, and one must set the maximum longitude
%   in M_MAP to LONO. For a global plot, this means one will call M_PROJ as
% 
%         M_PROJ(...,'longitudes',[LONO-360 LONO],...).
%
%   If you run CELLPLOT with M_MAP and data doesn't show up, make sure 
%   the longitudes have been set correctly.
%   __________________________________________________________________
%
%   See also LINESTYLE.
%
%   Usage: cellplot(x)
%          cellplot(x,y)
%          h=cellplot(x,y,index,'m_map');
%          cellplot(lono,lat,lon,index);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2014 J.M. Lilly --- type 'help jlab_license' for details
 
linestr=[];
str='matlab';

lono=[];
if ~iscell(varargin{1})&&length(varargin{1})==1
    lono=varargin{1};
    varargin=varargin(2:end);
end

while ischar(varargin{end})
    if ~isempty(strfind(varargin{end},'mat'))||~isempty(strfind(lower(varargin{end}),'m_m'))
        str=lower(varargin{end});
    else
        linestr=varargin{end};
    end
    varargin=varargin(1:end-1);
end

if iscell(varargin{1})
    if  ~iscell(varargin{end})
        ii=varargin{end};
        varargin=varargin(1:end-1);
    else
        ii=inf;
    end
end

x=varargin{1};
if length(varargin)==1
     y=[];
else
     y=varargin{2}; 
end
if length(varargin)<=2
    z=[];
else
    z=varargin{3};
end

% if length(varargin)==0
%     error('X is not a cell array.')
% end
% if ~iscell(varargin{1})
%     error('X is not a cell array.')
% end

if iscell(x)
    if isempty(ii)
        x=x(ii);
    else
        if isfinite(ii)
            x=x(ii);
            if ~isempty(y)
                y=y(ii);
            end
            if ~isempty(z)
                z=z(ii);
            end
        end
    end
end

if iscell(x)
   index=find(cellength(x)>0);
   x=x(index);
   if ~isempty(y)
       y=y(index);
   end
   if ~isempty(z)
       z=z(index);
   end
end

holdstate=ishold;
storestate=get(gcf,'BackingStore');

if ~holdstate
    cla
end



if ~isempty(x)
    if ~iscell(x)
        if size(x,2)~=1&&ndims(x)<3
            x=col2cell(mat2col(x));
            if ~isempty(y);
                 y=col2cell(mat2col(y));
            end
            if ~isempty(z);
                 z=col2cell(mat2col(z));
            end
        end
    end
    
    if iscell(x)
        cell2col(x);
        if ~isempty(y);
            cell2col(y);
        end
        if ~isempty(z);
            cell2col(z);
        end
    end
    
    %/****************************************************************
    %Unwrapping longitude
    if ~isempty(lono)
        infi=isinf(x);
        x=deg360(frac(360,2*pi)*unwrap(angle(rot(frac(2*pi,360)*(x-lono)))))+lono;
        x(infi)=inf;
        
        for i=1:size(x,2)
            while maxmax(x(:,i))>lono
                x(:,i)=x(:,i)-360;
            end
        end
        booljump=abs(x(2:end)-x(1:end-1))>90;
        x(booljump)=inf;
        if ~isempty(y)
            y(booljump)=inf;
        end
        if ~isempty(z);
            z(booljump)=inf;
        end
        xlim([lono-360 lono])
    end
    %/****************************************************************
    
    hold on
    nani=[0;find(isnan(x))];  %Don't forget the first one!
    h=ones(length(nani)-1,1);
    
    %x=datetime(x,'ConvertFrom','datenum');  %Just testing this...
    
    for i=1:length(nani)-1;
        a=nani(i)+1;
        b=nani(i+1)-1;
        if isempty(y)
            h(i)=plot(x(a:b));
        elseif ~isempty(z)
            h(i)=plot3(x(a:b),y(a:b),z(a:b));
        elseif ~isempty(y)&&isempty(z)
            if strcmpi(str(1:3),'mat')
                h(i)=plot(x(a:b),y(a:b));
            elseif  strcmpi(str(1:3),'m_m')
                h(i)=m_plot(x(a:b),y(a:b));
            end
        end 
    end
end
    
if ~isempty(linestr)
    linestyle(h,linestr);
else
    linestyle(h,'default');
end

if nargout==0
    clear h
end

if ~isempty(lono)
    xlim([lono-360 lono])
    ylim([minmin(y(isfinite(y))) maxmax(y(isfinite(y)))])
    boxon
end

set(gcf,'BackingStore',storestate)

if ~holdstate
    hold off
end

