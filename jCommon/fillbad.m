function[datao]=fillbad(data,i2,i3,i4)
%FILLBAD  Linearly interpolate over bad data points. 
%
%   DATA=FILLBAD(DATA,FLAG), where DATA is a column vector and "bad" points
%   within DATA are equal to FLAG, fills the bad points with a linear 
%   interpolation between the closest adjacent good points. 
%
%   FILLBAD(DATA) uses a value of NAN for FLAG.
%
%   If DATA is a matrix, interpolation is performed along each column. 
%
%   FILLBAD(DATA,FLAG,MAXN) or FILLBAD(DATA,INDEX,MAXN) interpolates over
%   gaps having a maximum length MAXN, which defaults to INF. 
%
%   For complex-valued data, the interpolation is performed on both the 
%   real and imaginary parts separately.  For a FLAG of NAN, this means 
%   the data must contain NAN+SQRT(-1)*NAN.
%
%   DATAC=FILLBAD(DATAC,...) also works when the input DATAC is a cell
%   array of numerical arrays, DATAC{1}=DATA1, DATAC{2}=DATA2, etc. The 
%   output is then a cell array of the identical size.
%
%   'fillbad --t' runs a test.
%
%   Usage: data=fillbad(data,flag);
%          data=fillbad(data,flag,maxnan);   
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2000--2015 J.M. Lilly --- type 'help jlab_license' for details  

%
%   FILLBAD(...,STR) uses interpolation method specified by STR; see
%   INTERP1 for details.

if strcmpi(data, '--t')
  fillbad_test,return
end
    
    
flag=nan;
maxn=inf;
str='linear';

if nargin==2
   if ischar(i2)
       str=i2;
   else
       flag=i2;
   end
elseif nargin==3
   flag=i2;
   if ischar(i3)
       str=i3;
   else
       maxn=i3;
   end
elseif nargin==4
    flag=i2;
    maxn=i3;
    str=i4;
end



%corrected this to ignore missing values at beginning and end
%really should work column by column....
%changed to work with nans
if iscell(data)
    datao=data;
    for i=1:length(data)
        datao{i}=fillbad_loop(data{i},flag,maxn,str);
    end
else
    datao=fillbad_loop(data,flag,maxn,str);
end

function[datao]=fillbad_loop(data,flag,maxn,str)
datao=data;

for i=1:size(data,2)
  if all(isreal(data(:,i)))
     datao(:,i)=fillbad1(data(:,i),flag,maxn,str);
  else
    datar=fillbad1(real(data(:,i)),flag,maxn,str);
    datai=fillbad1(imag(data(:,i)),flag,maxn,str);
    datao(:,i)=datar+sqrt(-1)*datai;
  end
end

function[data]=fillbad1(data,flag,maxn,str)
  
if ~isnan(flag)
   index=find(data==flag);
else 
   index=find(isnan(data));
end

[m,n]=size(data);

%Imply bad blocks outside of time series range
index=[0;index;m+1];

if ~isempty(index)
%	data=data(:);
	ii=find(diff(index)~=1);
	if ~isempty(ii) 
  	  ia=[index(1)-1;index(ii+1)-1];%Last data before bad points
  	  ib=[index(ii)+1;index(end)+1];%First data after bad points
      
      %Exclude beginning and end blocks
      ia=ia(2:end-1);
      ib=ib(2:end-1);
      
  	  lblocks=ib-ia-1;
      if lblocks>0
          if ia(1)==0
             ia=ia(2:end);
             ib=ib(2:end);
             lblocks=lblocks(2:end);
          end	
          if ib(end)==length(data)+1
             ia=ia(1:end-1);
             ib=ib(1:end-1);
             lblocks=lblocks(1:end-1);
          end	
          for i=1:length(lblocks)
             if lblocks(i)<=maxn
                temp=interp1([ia(i) ib(i)],[data(ia(i)) data(ib(i))],ia(i)+1:ib(i)-1);
                data(ia(i)+1:ib(i)-1)=temp;
             end
          end
      end
	end
%	data=reshape(data,m,n);
end


function[]=fillbad_test


tol=1e-10;


x=[1 nan nan nan nan 6]';
y=[1 2    3   4   5  6]';

b=aresame(y,fillbad(x,nan),tol);
reporttest('FILLBAD one column one block',b);


b=aresame([y y y],fillbad([x x x],nan),tol);
reporttest('FILLBAD multiple columns one block',b);

x=[1 nan 3   nan nan 6]';
y=[1 2    3   4   5  6]';

b=aresame(y,fillbad(x,nan));
reporttest('FILLBAD one column multiple blocks',b);


b=aresame([y y y],fillbad([x x x],nan));
reporttest('FILLBAD multiple columns multiple blocks',b);

