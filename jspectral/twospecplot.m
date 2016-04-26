function[h]=twospecplot(varargin)
%TWOSPECPLOT  Plots a pair of rotary or Cartesian spectra.
%
%   TWOSPECPLOT(F,SPP,SNN) plots the rotary spectra SPP and SNN with log- 
%   log axes.  SNN is on the left, with a reversed x-axis. 
%
%   TWOSPECPLOT(F,SXX,SYY,'cartesian') left the same, but now labels the 
%   and right axes as the zonal and meridional spectra, respectively.
%
%   TWOSPECPLOT(...,'tides') marks the major tidal with dotted lines.  F is
%   assumed to be in units of *radians* per hour.
%
%   TWOSPECPLOT(LAT,F,SPP,SNN) also marks the Coriolis frequency at 
%   latitude LAT, using a red line, on both sides of the spectrum.
%
%   TWOSPECPLOT(...,STR) sets the line style by passing the style string
%   STR to LINESTYLE.
%
%   H=TWOSPECPLOT(...) also returns the handle to the subplots.
%
%   The input quantities F, SPP, and SNN or F, SXX, and SYY may be cell
%   arrays.  See 'Cell array input / output' under MSPEC for details.
%
%   See also MSPEC, SLEPTAP.
%
%   Usage: twospecplot(f,spp,snn)
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2012--2015 J.M. Lilly --- type 'help jlab_license' for details

str='rotary';
xlims=[];
ylims=[];
linestr='notides';
stylestr=[];

for i=1:4
    if ischar(varargin{end-1})&&~ischar(varargin{end})
        if strcmpi(varargin{end-1}(1:3),'xli')
            xlims=varargin{end};
        elseif strcmpi(varargin{end-1}(1:3),'yli')
            ylims=varargin{end};
        end
        varargin=varargin(1:end-2); 
    elseif ischar(varargin{end})
        if length(varargin{end})<3
             stylestr=varargin{end};
        else
            if strcmpi(varargin{end}(1:3),'tid')||strcmpi(varargin{end}(1:3),'not')
                linestr=varargin{end};
            elseif strcmpi(varargin{end}(1:3),'rot')||strcmpi(varargin{end}(1:3),'car')
                str=varargin{end};
            else
                stylestr=varargin{end};
            end
        end
        varargin=varargin(1:end-1);
    end
end

if length(varargin)==3
    lat=[];
else
    lat=varargin{1};
    varargin=varargin(2:end);
end
f=varargin{1};
spp=varargin{2};
snn=varargin{3};


if isempty(xlims)
    if ~iscell(f)
        xlims=[minmin(f) maxmax(f)];
    else
        xlims=[minmin(cell2col(f)) maxmax(cell2col(f))];
    end
end
if isempty(ylims)
    if ~iscell(snn)
        ylims=[min(minmin(snn),minmin(spp))*0.9 max(maxmax(snn),maxmax(spp))*1.1];
    else
        ylims=[min(minmin(cell2col(snn)),minmin(cell2col(spp)))*0.9 max(maxmax(cell2col(snn)),maxmax(cell2col(spp)))*1.1];
    end
end

h(1)=subplot(1,2,1);
if ~iscell(snn)
    hl=plot(f,snn);
else
    hl=cellplot(f,snn);
end
if ~isempty(stylestr)
    linestyle(hl,stylestr);
end
set(gca,'xdir','reverse'),xlog,ylog,axis tight,ylim(ylims),xlim(xlims),boxon,hold on
if strcmpi(linestr(1:3),'tid')
    vlines(tidefreq,':'),
end
if ~isempty(lat)
    vlines(abs(corfreq(lat)),'r')
end
if strcmpi(str(1:3),'rot')
    title('Negative Rotary Spectra')
else
    title('Zonal Spectra')
end

h(2)=subplot(1,2,2);
if ~iscell(spp)
    hl=plot(f,spp);
else
    hl=cellplot(f,spp);
end
if ~isempty(stylestr)
    linestyle(hl,stylestr);
end
set(gca,'YAxisLocation','right'),xlog,ylog,axis tight,ylim(ylims),xlim(xlims),boxon,hold on
if strcmpi(linestr(1:3),'yes')
    vlines(tidefreq,':'),
end
if ~isempty(lat)
    vlines(abs(corfreq(lat)),'r')
end
if strcmpi(str(1:3),'rot')
    title('Positive Rotary Spectra')
else
    title('Meridional Spectra')
end
%h=packfig(1,2,'columns');

if nargout==0
    clear h
end
