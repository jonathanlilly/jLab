function[lat1,lon1,cv1]=noisedrifters(varargin)
%NOISEDRIFTERS  Create a noise Lagrangian dataset matching mean and variance.
%
%   Given Lagrangian trajectories and their velocity spectra, NOISEDRIFTERS
%   creates a noise dataset that matches, for each trajectory, (i) the
%   starting location, (ii) the mean velocity, and (iii) the approximate
%   variance or eddy kinetic energy, with a velocity spectrum that is an 
%   isotropic version of the spectrum of the input trajectory.
%
%   [LATN,LONN,CVN]=NOISEDRIFTERS(NUM,LAT,LON,CV,SPP,SNN), where LAT and
%   LON are the latitudes and longitudes of Lagrangian trajectories
%   observed at Matlab date number NUM, and with complex velocities CV 
%   rotary velocity spectra SPP and SSN, outputs a noise dataset of 
%   trajectories with latitudes LATN, longitudes LONN, and velocities CVN.
%
%   All input arguments are cell arrays of the same size, with one
%   trajectory per cell, and all output arrays will also be of this size. 
%
%   Note that the units of CV are cm/s, as are those of CVN.  Similarly,
%   the spectral SPP and SNN should be computed with CV having those units.
%
%   NOISEDRIFTERS creates an isotropic spectrum by taking, at each
%   frequency, the minimum of SPP and SNN.  These are as output by MSPEC, 
%   and should be formed with a sufficient degree of smoothing, e.g. not
%   the periodogram estimate.  
%
%   Once the spectral shape is established, NOISEDRIFTERS creates a random
%   Gaussian time series having exactly this spectral shape.  Each random 
%   velocity time series is then set to have the same mean value and 
%   variance as the corresponding original time series. 
% 
%   The random velocity time series are then integrated to give the output
%   trajectories LATN and LONN, using UV2LATLON.  The initial value of LAT 
%   and LON are the same as that of LATN and LONN.  
%
%   These noise trajectories are differenced again using LATLON2UV to
%   produce CVN, as this is how the velocities are produced from 
%   trajectories for the observations.  Because of the differences between 
%   numerically integrating with UV2LATLON and differencing with LATLON2UV,
%   the variancees of CVN and CV are approximately but not exactly equal.
%
%   NOISEDRIFTERS(NUM,LAT,LON,CV) with no spectra input alternately uses
%   white noise velocities to generate the output fields.
%
%   NOISEDRIFTERS(...,'parallel') parallelizes the computation using a 
%   PARFOR loop, which requires the Parallel Computing Toolbox.
%
%   Usage: [latn,lonn,cvn]=noisedrifters(num,lat,lon,cv,spp,snn);
%          [latn,lonn,cvn]=noisedrifters(num,lat,lon,cv);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2019 J.M. Lilly --- type 'help jlab_license' for details

 %   'noisedrifters --t' runs a test.

%if strcmp(varargin{1}, '--t')
%    noisedrifters_test,return
%end
 
%function[]=noisedrifters_test
str='ser';
if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
end
num=varargin{1};
lat=varargin{2};
lon=varargin{3};
cv=varargin{4};
spp=[];
snn=[];
if length(varargin)>4
    spp=varargin{5};
    snn=varargin{6};
end
    
noise=cell(length(lat),1);
lat1=cell(length(lat),1);
lon1=cell(length(lat),1);
cv1=cell(length(lat),1);
smin=cell(length(lat),1);

if strcmp(str(1:3),'ser')
    %serial version
    %----------------------------------------------------------------------
    if isempty(spp)&&isempty(snn)
        for i=1:length(lat)
            noise{i}=randn(size(cv{i}))+1i*randn(size(cv{i}));
        end
    else
        %Shape white noise by square root of minimum rotary spectra
        for i=1:length(lat)
            smin{i}=min([spp{i} snn{i}],[],2);
            if iseven(length(lat{i}))
                spn=[smin{i};flipud(smin{i}(2:end-1,:))];
            else
                spn=[smin{i};flipud(smin{i}(2:end,:))];
            end
            noisei=randn(size(cv{i}))+1i*randn(size(cv{i}));
            noise{i}=ifft(sqrt(spn).*fft(noisei,[],1),[],1);
        end
    end
    
    %Set white noise to have same variance and mean as the real trajectory
    for i=1:length(lat)
        noise{i}=noise{i}.*vstd(cv{i},1)./vstd(noise{i},1);
        phase=rot(angle(mean(cv{i},1))-angle(mean(noise{i},1)));
        noise{i}=noise{i}.*phase; clear phase
        noise{i}=noise{i}-mean(noise{i},1);
        noise{i}=noise{i}+mean(cv{i},1);
    end
    
    %aresame(cellstd(noise),cellstd(cv),1e-6)   %true
    %aresame(cellmean(noise),cellmean(cv),1e-6) %true
    
    for i=1:length(noise)
        [latf,lonf]=...
            uv2latlon(num{i},real(noise{i}),imag(noise{i}),lat{i}(1),lon{i}(1),'forward');
        %[latb,lonb]=...
        %    uv2latlon(num{i},real(noise{i}),imag(noise{i}),lat{i}(end),lon{i}(end),'forward');
        lat1{i,1}=latf;
        lon1{i,1}=lonf;
        %lat1{i}=latb;
        %lon1{i}=lonb;
        %   [xf,yf,zf]=latlon2xyz(latf,lonf);
        %    [xb,yb,zb]=latlon2xyz(latb,lonb);
        %  [lat1{i},lon1{i}]=xyz2latlon((xf+xb)./2,(yf+yb)./2,(zf+zb)./2);
        %[lat1{i},lon1{i}]=xyz2latlon((xf+xb)./2,(yf+yb)./2,(zf+zb)./2);
        %this doesn't seem to help
        %     [lat1{i},lon1{i}]=...
        %         uv2latlon(num{i},real(noise{i}),imag(noise{i}),lat{i}(1),lon{i}(1));
        
        cv1{i,1}=latlon2uv(num{i},lat1{i},lon1{i});
        %cv1{i}=latlon2uv(num{i},lat1{i},lon1{i},'forward');
    end
else
    %parallel version
    %----------------------------------------------------------------------
    parfor i=1:length(lat)
        disp(['NOISEDRIFTERS working on trajectory ' int2str(i) ' of ' int2str(length(lat)) '.'])
        if isempty(spp)&&isempty(snn)
            noise{i}=randn(size(cv{i}))+1i*randn(size(cv{i}));
        else
            %Shape white noise by square root of minimum rotary spectra
            smin{i}=min([spp{i} snn{i}],[],2);
            if iseven(length(lat{i}))
                spn=[smin{i};flipud(smin{i}(2:end-1,:))];
            else
                spn=[smin{i};flipud(smin{i}(2:end,:))];
            end
            noisei=randn(size(cv{i}))+1i*randn(size(cv{i}));
            noise{i}=ifft(sqrt(spn).*fft(noisei,[],1),[],1);
        end
        
        %Set white noise to have same variance and mean as the real trajectory
        
        noise{i}=noise{i}.*vstd(cv{i},1)./vstd(noise{i},1);
        phase=rot(angle(mean(cv{i},1))-angle(mean(noise{i},1)));
        noise{i}=noise{i}.*phase;
        noise{i}=noise{i}-mean(noise{i},1);
        noise{i}=noise{i}+mean(cv{i},1);
        
        [latf,lonf]=...
            uv2latlon(num{i},real(noise{i}),imag(noise{i}),lat{i}(1),lon{i}(1),'forward');
        %[latb,lonb]=...
        %    uv2latlon(num{i},real(noise{i}),imag(noise{i}),lat{i}(end),lon{i}(end),'forward');
        lat1{i}=latf;
        lon1{i}=lonf;
        %lat1{i}=latb;
        %lon1{i}=lonb;
        %   [xf,yf,zf]=latlon2xyz(latf,lonf);
        %    [xb,yb,zb]=latlon2xyz(latb,lonb);
        %  [lat1{i},lon1{i}]=xyz2latlon((xf+xb)./2,(yf+yb)./2,(zf+zb)./2);
        %[lat1{i},lon1{i}]=xyz2latlon((xf+xb)./2,(yf+yb)./2,(zf+zb)./2);
        %this doesn't seem to help
        %     [lat1{i},lon1{i}]=...
        %         uv2latlon(num{i},real(noise{i}),imag(noise{i}),lat{i}(1),lon{i}(1));
        
        cv1{i}=latlon2uv(num{i},lat1{i},lon1{i});
        %cv1{i}=latlon2uv(num{i},lat1{i},lon1{i},'forward');
    end
end
