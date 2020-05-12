function[varargout]=about_jtopo(varargin)
%ABOUT_JTOPO  One-sixth degree global topography, from Smith and Sandwell + IBCAO.
%   _______________________________________________________________________
%   
%   *|* jtopo.png --- Figure illustrating the JTOPO topography dataset.  
%   Type 'jhelp about_jtopo' to view this image. *|*
%   _______________________________________________________________________
%
%   JTOPO is a matfile containing smoothed one-twelfth degree global 
%   topography based on the Smith and Sandwell database together with the
%   International Bathymetric Chart of the Arctic Ocean (IBCAO).
%
%   LOAD JTOPO loads the structure JTOPO, with the following fields:
%
%       jtopo.about     Pointer to this document  
%       jtopo.lat       Array of latitudes         [2048 x 1] 
%       jtopo.lon       Array of longitudes        [1 x 4320]
%       jtopo.topo      Matrix of topography       [2048 x 4320]
%
%   Typing 'use jtopo' maps these fields into variables in the current
%   workspace, e.g. 'use jtopo, pcolor(lon,lat,topo), shading interp'.
%
%   TOPO is in units of kilometers and is positive for above sea level, 
%   and negative for below sea level.  
%
%   LAT is uniformly spaced from -80.666 to  89.917, and LON is uniformly
%   spaced from -180 to 179.912.  These are *grid-centered* values, that 
%   is, they indicate midpoints of the topography cells.
%
%   JTOPO is distributed with JLAB, available at http://www.jmlilly.net.
%
%   See also READTOPO, which reads in any region of the Smith and Sandwell 
%   data, ABOUT_IBCAO, and TOPOPLOT, which makes simple plots using JTOPO.
%   __________________________________________________________________
%  
%   Processing
%
%   The one-minute Smith and Sandwell data, and 1/2 minute IBCAO data, are
%   smoothed to one-twelfth of a degree by averaging in 1/12 x 1/12 bins. 
%
%   Smith and Sandwell is defined from -80.738 and 80.738, while IBCAO is 
%   defined from 64 N to 90 N.  In the overlap region, from 64 to 84.738 N,
%   the root-mean-square different between the two datasets is about 50 m.
%
%   JTOPO uses a linear blend to resolve the small discrepencies in the 
%   overlap region.  From 64 to 84.738 N, JTOPO transitions from being all
%   Smith and Sandwell, to all IBCAO, with a weighted average in between.
%   __________________________________________________________________
%
%   Data and documentation
%
%   This dataset is based on the Smith and Sandwell Global Topography 
%   Dataset v. 18.1 and IBCAO v. 3.0, which are included with JDATA.
%
%   The source and reference for the Smith and Sandwell dataset are
%
%       http://topex.ucsd.edu/WWW_html/mar_topo.html
%
%      Smith, W. H. F., and D. T. Sandwell, Global seafloor topography 
%         from satellite altimetry and ship TOPO soundings, Science, 
%         v. 277, p. 1957-1962, 26 Sept., 1997.
%
%   The source and reference for the IBCAO dataset are
%
%      http://www.ngdc.noaa.gov/mgg/bathymetry/arctic/grids/version3_0/
%
%      Jakobsson, M., L. A. Mayer, B. Coakley, J. A. Dowdeswell, S. Forbes,
%          B. Fridman, H. Hodnesdal, R. Noormets, R. Pedersen, M. Rebesco,
%          H.-W. Schenke, Y. Zarayskaya A, D. Accettella, A. Armstrong, 
%          R. M. Anderson, P. Bienhoff, A. Camerlenghi, I. Church, 
%          M. Edwards, J. V. Gardner, J. K. Hall, B. Hell, O. B. Hestvik, 
%          Y. Kristoffersen, C. Marcussen, R. Mohammad, D. Mosher, 
%          S. V. Nghiem, M. T. Pedrosa, P. G. Travaglini, and 
%          P. Weatherall, The International Bathymetric Chart of the Arctic
%          Ocean (IBCAO) Version 3.0, Geophysical Research Letters, 
%          doi: 10.1029/2012GL052219
%   __________________________________________________________________
%
%   License and Copyright 
%
%   JTOPO.MAT is distributed with JDATA for RESEARCH AND NON-PROFIT USE 
%   ONLY, in accordance with the copyright statement for the Smith and 
%   Sandwell dataset.  For details, see TOPO_COPYRIGHT.
%
%   No copyright or policy is specified in the IBCAO documentation. 
%   __________________________________________________________________
%
%   Dataset creation
%
%   For completeness, the m-file ABOUT_JTOPO also contains the processing 
%   steps used in the creation of JTOPO.MAT.  
%
%   If you wish to do this yourself, with JLAB on your search path, 
%   'about_jtopo --create' will recreate the JTOPO.MAT dataset by reading  
%   in and averaging the two topographic datasets. This will take a while.
%
%   For this to work you will need to have the JDATA folder containing
%   file 'topo_19.1.img' and 'ibcao.mat' downloaded and on your Matlab 
%   search path. 
%   __________________________________________________________________
%
%   See also READTOPO, TOPOPLOT, JDATA.
%
%   'about_jtopo --f' generates the sample figure shown above.
%
%   Usage: about_jtopo
%          about_jtopo --create
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2020 J.M. Lilly --- type 'help jlab_license' for details
 
if nargin==0
    help about_jtopo
elseif nargin>0
    if strcmpi(varargin{1}, '--create')
        jtopo_create,return
    elseif strcmpi(varargin{1}, '--f')
        type makefigs_jtopo
        makefigs_jtopo;
        return
    end
end

function[]=jtopo_create
      
%/************************************************************************
%Main Smith and Sandwell portion
%delta=1/6;
delta=1/12;

lon=[-180-delta/2:delta:180-delta/2];
lat=[-80.66666-delta/2:delta:80.67+delta/2]';

topo=zeros(length(lat)-1,length(lon)-1);
for i=1:length(lon)-1  
    disp(['Longitude band number ' num2str(i) ' of ' num2str(length(lon)-1) '.' ])
    [topo1,lat1,lon1]=readtopo([lon(i) lon(i)+delta -80.738 80.738]);
    topo1=vmean(topo1,2);  %Average across longitudes
    [mz,latbin,lonbin]=twodstats(lat1,lon(i)+0*lat1+delta/2,topo1,lat,[lon(i) lon(i)+delta]);%Average across latitudes
    topo(:,i)=mz';
end

lat=latbin;
lon=lon(1:end-1)+delta/2;

about='For more information, type ''about_jtopo''.';
matsave jtopo about lat lon topo
%matsave jtopo_onetwelfth about lat lon topo
%\************************************************************************

%/************************************************************************
%Adding Arctic dataset
lon=[-180-delta/2:delta:180-delta/2];
lat=[-80.66666-delta/2:delta:90]';

newtopo=zeros(length(lat)-1,length(lon)-1);
newtopo(1:size(jtopo.topo,1),:)=jtopo.topo;

load ibcao
ibcaotopo=zeros(size(newtopo));
[long,latg]=meshgrid(ibcao.lon,ibcao.lat);
for i=1:length(lon)-1  
    disp(['Longitude band number ' num2str(i) ' of ' num2str(length(lon)-1) '.' ])
    index=find(ibcao.lon>lon(i) & ibcao.lon< lon(i)+delta);
    [mz,latbin,lonbin]=twodstats(latg(:,index),long(:,index),ibcao.topo(:,index),...
        lat(min(find(lat>64)):end),[lon(i) lon(i)+delta]);%Average across latitudes
    ibcaotopo(min(find(lat>64)):end,i)=mz';
end

%topodiff=jtopo.topo(min(find(lat>64)):end,:)-ibcaotopo(1:100,:);
%figure,plot(1000*vmean(abs(topodiff),2))
%sqrt(mean(squared(topodiff(:))))*1000 = 52 meters
lat=lat(1:end-1)+delta/2;
lon=lon(1:end-1)+delta/2;

%Linearly blend between the two datasets
index=min(find(lat>64)):find(lat<max(jtopo.lat),1,'last');
weight=ones(size(lat));
weight(index)=1-[0:length(index)-1]./(length(index)-1);
weight(index(end)+1:end)=0;

%figure,plot(lat,weight),vlines([64 80.6666666])
weight=vrep(weight,size(newtopo,2),2);
topo=newtopo.*weight+ibcaotopo.*(1-weight);
%\************************************************************************

matsave jtopo about lat lon topo
%matsave jtopo_onetwelfth about lat lon topo
disp('JTOPO creation complete.')

