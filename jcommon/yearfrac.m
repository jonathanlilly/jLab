function[yf,mf]=yearfrac(num)
%YEARFRAC  Converts a DATENUM into 'year.fraction' and 'month.fraction'.
%  
%   YF=YEARFRAC(NUM) where NUM is an array of dates in Matlab's 'datenum'
%   format, returns the fraction of the year at each date.
%
%   FLOOR(YF) returns the year.  Note that the actual number of days in 
%   each year is used, including leap years.
%
%   YF=YEARFRAC(NUM) where NUM is a cell array of numeric arrays, also
%   works.  YF is then a cell array of the same size as NUM.
%
%   [YF,MF]=YEARFRAC(NUM) also returns MF, the fraction of the current 
%   month at each date.  FLOOR(MF) is the standard month number. 
%  
%   See also DATENUM, DATEVEC.
%
%   Usage: yf=yearfrac(num);
%          [yf,mf]=yearfrac(num);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 1998--2016 J.M. Lilly --- type 'help jlab_license' for details        
  
if ~iscell(num)
    if strcmpi(num, '--t')
        yf2num('--t'),return
    end  
end

if iscell(num)
    for i=1:length(num)
        if ~isempty(num{i})
            [yf{i,1},mf{i,1}]=yearfrac1(num{i});
        else
            yf{i}=[];
            mf{i}=[];
        end
    end
else
    [yf,mf]=yearfrac1(num);
end

function[yf,mf]=yearfrac1(num)
yf=[];
if ~isempty(num)
    index=find(isnan(num));
    if ~isempty(index)
        num(index)=0;
    end

    [y,mo,d,h,mi,s] = datevec(num);
    
    %Number of days in this year?
    nd=datenum(y,12,31)-datenum(y-1,12,31);
    na=datenum(y,mo,d,h,mi,s)-datenum(y-1,12,31)-1;
    %The minus one is because for Jan 1, I add nothing to year.fraction

    yf=y+na./nd;
    
    %Number of days in this month?
    nd=datenum(y,mo+1,1)-datenum(y,mo,1);  %Yes, this also works if mo = 12, for January
    mf=mo+((d-1)+h/24)./nd;
    
    if ~isempty(index)
      yf(index)=nan;
      mf(index)=nan;
    end
end

function[num]=yf2num(yf)
%YF2NUM  Convert date in 'year.fraction' format to 'datenum' format.
%
%   YF2NUM(YF) where YF is an array of dates in 'year.fraction' format
%   returns the array in Matlab's 'datenum' format.
%
%   See also YEARFRAC, DATENUM, DATEVEC.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2009 J.M. Lilly --- type 'help jlab_license' for details  
    
if strcmpi(yf, '--t')
  yf2num_test,return
end

y=floor(yf);

%Number of days in this year?
d0=datenum(y-1,12,31);
d1=datenum(y,12,31);
nd=d1-d0;
num=d0+nd.*(yf-y)+1;

function[]=yf2num_test

yearf=(1850:(1/360):2004)';
num=yf2num(yearf);
yearf2=yearfrac(num);
bool=maxmax(abs(yearf-yearf2))<1e-10;
reporttest('YF2NUM and YEARFRAC, daily resolution, 1e-10 cutoff',bool)

yearf=(1990:(1/360/24):2004)';
%yearf=(1850:(1/360/24):2004)';
num=yf2num(yearf);
yearf2=yearfrac(num);
bool=maxmax(abs(yearf-yearf2))<1e-10;
reporttest('YF2NUM and YEARFRAC, hourly resolution, 1e-10 cutoff',bool)

  
  
