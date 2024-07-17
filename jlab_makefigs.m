function[varargout]=jlab_makefigs(namestr,str)
%JLAB_MAKEFIGS  Make figures for papers by J. M. Lilly.
%
%   JLAB_MAKEFIGS NAME makes all figures for the publication NAME. 
%   Please note, you may need to download addition datafiles or software 
%   for these to run correctly; see MAKEFIGS_NAME for details.
%
%   Supported publications are as follows.
%
%  'jlab_makefigs stokes':
%   Lilly, J. M., J. Feske, B. Fox-Kemper, and J. J. Early (2024). 
%       Integral theorems for the gradient of a vector field, with a fluid
%       dynamical application.  Proceedings of the Royal Society of London, 
%       Series A. 480 (2293): 20230550, 1–30.
%
%  'jlab_makefigs gulfcensus':
%   Lilly, J. M. and P. Perez-Brunius (2021b).  Extracting statistically 
%       significant eddy signals from large Lagrangian datasets using 
%       wavelet ridge analysis, with application to the Gulf of Mexico.
%       Nonlinear Processes in Geophysics, 28: 181–212.
%
%  'jlab_makefigs gulfdrifters':
%   Lilly, J. M. and P. Perez-Brunius (2021a). A gridded surface current
%       product for the Gulf of Mexico from consolidated  drifter
%       measurements.  Earth System Science Data, 13: 645–669.
%
%  'jlab_makefigs transfer':
%   Lilly, J. M. and S. Elipot (2021). A unifying perspective on transfer
%       function solutions to the unsteady Ekman problem.  Fluids, 6 (2): 
%       85, 1--36. 
%
%  'jlab_makefigs kinematics':
%   Lilly, J. M. (2018) Kinematics of a fluid ellipse in a linear flow. 
%       Fluids, 3 (1) 16: 1--29.
%
%  'jlab_makefigs matern':
%   Lilly, J. M., A. M. Sykulski, J. J. Early, and S. C. Olhede (2017).
%       Fractional Brownian motion, the Matern process, and stochastic 
%       modeling of turbulent dispersion.  Nonlinear Processes in
%       Geophysics, 24: 481--514.
%
%   'jlab_makefigs element':
%   Lilly, J. M.  (2017).  Element analysis: a wavelet-based method for
%       analyzing time-localized events in noisy time series.  Proceedings 
%       of the Royal Society of London, Series A, 473 (2200), 1--28.
%
%   'jlab_makefigs superfamily':
%   Lilly, J. M., and S. C. Olhede (2012). Generalized Morse wavelets as
%      a superfamily of analytic wavelets.  IEEE Transactions on Signal
%      Processing 60 (11), 6036--6041.
%
%   'jlab_makefigs multivariate':
%   Lilly, J. M., and S. C. Olhede (2012). Analysis of modulated 
%      multivariate oscillations. IEEE Transactions on Signal Processing, 
%      60 (2), 600--612.
%
%   'jlab_makefigs vortex': 
%   Lilly, J. M., R. K. Scott, and S. C. Olhede (2011). Extracting waves 
%      and vortices from Lagrangian trajectories. Geophysical Research 
%      Letters, 38, L23605, 1--5.
%
%   'jlab_makefigs trivariate':
%   Lilly, J. M. (2011). Modulated oscillations in three dimensions. IEEE
%      Transactions on Signal Processing, 59 (12), 5930--5943.
%
%   'jlab_makefigs analytic':
%   Lilly, J. M., and S. C. Olhede (2010). On the analytic wavelet 
%      transform. IEEE Transactions on Information Theory, 56 (8),
%      4135--4156.
%
%   'jlab_makefigs bandwidth':
%   Lilly, J. M., and S. C. Olhede (2010). Bivariate instantaneous 
%      frequency and bandwidth. IEEE Transactions on Signal Processing,
%      58 (2), 591--603.
%
%   'jlab_makefigs asilomar':
%   Lilly, J. M., and S. C. Olhede (2009). Wavelet ridge estimation of 
%      jointly modulated multivariate oscillations. IEEE 43rd Asilomar 
%      Conference on Signals, Systems, and Computers, 452--456. 
%
%   'jlab_makefigs morsies':
%   Lilly, J. M., and S. C. Olhede (2009). Higher-order properties of 
%      analytic wavelets. IEEE Transactions on Signal Processing, 57 (1),
%      146--160.
%   
%   'jlab_makefigs ridges':
%   Lilly, J. M., and J.-C. Gascard (2006). Wavelet ridge diagnosis of 
%      time-varying elliptical signals with application to an oceanic
%      eddy. Nonlinear Processes in Geophysics 13, 467--483.
%   
%   'jlab_makefigs all' makes the figures for all papers; this will take
%    a while, and requires the needed files and software to be installed.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2023 J.M. Lilly --- type 'help jlab_license' for details     
 

if nargin==1
    str='noprint';
end

jj=0;
jj=jj+1;names{jj}='ridges';
jj=jj+1;names{jj}='morsies';
jj=jj+1;names{jj}='asilomar';
jj=jj+1;names{jj}='bandwidth';
jj=jj+1;names{jj}='analytic';
jj=jj+1;names{jj}='trivariate';
jj=jj+1;names{jj}='vortex';
jj=jj+1;names{jj}='multivariate';
jj=jj+1;names{jj}='superfamily';
jj=jj+1;names{jj}='element';
jj=jj+1;names{jj}='matern';
jj=jj+1;names{jj}='kinematics';
jj=jj+1;names{jj}='transfer';
jj=jj+1;names{jj}='gulfdrifters';
jj=jj+1;names{jj}='gulfcensus';
jj=jj+1;names{jj}='stokes';

%cd(jlab_settings('dirnames.figures'))

dti=get(0,'defaultTextInterpreter');
dli=get(0,'defaultTextInterpreter');
set(0,'defaultTextInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')


if strcmpi(namestr(1:3),'all')
    for i=1:length(names)
        eval(['jlab_makefigs_' names{i} '(''' str ''');'])
    end
else
    if exist(['jlab_makefigs_' namestr])
        eval(['jlab_makefigs_' namestr '(''' str ''');'])
    else
        disp(['Sorry, formula for paper ''' namestr ''' not found.'])
    end
end

set(0,'defaultTextInterpreter',dti)
set(0,'defaultLegendInterpreter',dli)
%-----------------------------------------------------------------------
function[varargout]=jlab_makefigs_stokes(str) 
makefigs_stokes(str);
%-----------------------------------------------------------------------
function[varargout]=jlab_makefigs_gulfcensus(str) 
makefigs_gulfcensus(str);
%-----------------------------------------------------------------------
function[varargout]=jlab_makefigs_gulfdrifters(str) 
makefigs_gulfdrifters(str);
%-----------------------------------------------------------------------
function[varargout]=jlab_makefigs_transfer(str) 
makefigs_transfer(str);
%-----------------------------------------------------------------------
function[]=jlab_makefigs_kinematics(str) 
makefigs_kinematics(str);
%-----------------------------------------------------------------------
function[varargout]=jlab_makefigs_matern(str) 
makefigs_matern(str);
%-----------------------------------------------------------------------
function[varargout]=jlab_makefigs_element(str) 
makefigs_element(str);
%-----------------------------------------------------------------------
function[varargout]=jlab_makefigs_superfamily(str)
makefigs_superfamily(str);
%-----------------------------------------------------------------------
function[varargout]=jlab_makefigs_multivariate(str)
makefigs_multivariate(str);
%-----------------------------------------------------------------------
function[varargout]=jlab_makefigs_vortex(str)
makefigs_vortex(str);
%-----------------------------------------------------------------------
function[varargout]=jlab_makefigs_trivariate(str)
makefigs_trivariate(str);
%-----------------------------------------------------------------------
function[]=jlab_makefigs_analytic(str)
makefigs_analytic(str);
%-----------------------------------------------------------------------
function[]=jlab_makefigs_bandwidth(str)
makefigs_bandwidth(str);
%-----------------------------------------------------------------------
function[varargout]=jlab_makefigs_asilomar(str)
makefigs_asilomar(str);
%-----------------------------------------------------------------------
function[varargout]=jlab_makefigs_morsies(str)
makefigs_morsies(str);
%-----------------------------------------------------------------------
function[]=jlab_makefigs_ridges(str)
makefigs_ridges(str);

