function[cax]=colorquant(a,b)
%COLORQUANT  Set the color axis according to quantiles of the data.
%
%   COLORQUANT(Q) sets the limits of the color axis CAXIS to the Qth and 
%   the 100-Qth data quantiles.  In other words, Q percent of the data will
%   be below, and Q percent below, the limits of the color axis.
%
%   This can save you a lot of time versus manually setting color axes.
%
%   COLORQUANT(QA,QB) sets the limits to the QAth and QBth quantiles.
%   
%   Note that the parenthesis are optional when the input arguments are 
%   numbers, so one can type 'colorquant 10.5' or 'colorquant 10 95'.
%
%   COLORQUANT with no input arguments uses the default value of Q=0.1.
%
%   COLORQUANT works by grabbing all surface, contour, and scatter data 
%   values on the current set of axes. 
%
%   CAX=COLORQUANT(...) returns the value of the color axis.
%
%   Usage: colorquant
%          colorquant(Q)
%          colorquant(Q1,Q2)
%          colorquant 10
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2018--2022 J.M. Lilly --- type 'help jlab_license' for details

ao=0.1;
if nargin==1||nargin==2
    if ischar(a)
        a=str2num(a);
    end
end
if nargin==2
    if ischar(b)
        b=str2num(b);
    end
end

if nargin==0
    a=ao;b=100-ao;
elseif nargin==1
    b=100-a;
end

%a,b

h=gca;
hc=h.Children;
z=[];
for i=1:length(hc)
    if strcmpi(hc(i).Type,'Contour')
        zi=hc(i).ZData;
        z=[z;zi(:)];
    elseif strcmpi(hc(i).Type,'Surface')||strcmpi(hc(i).Type,'Scatter')
        zi=hc(i).CData;
        z=[z;zi(:)];
    end
end

%size(z)
z=z(isfinite(z));
z=sort(z(:));
%figure,plot(z)

ii=[1:length(z)]/length(z);
%[a, b]
ai=find(ii>a/100,1,'first');
bi=find(ii<b/100,1,'last');
%[ai bi]
%figure,plot(z),vlines([ai bi])
%z(ai), z(bi)
%ai, bi
cax=[z(ai) z(bi)];
caxis(cax)
if nargout==0
    clear cax
end

