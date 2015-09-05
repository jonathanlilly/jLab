function[teps]=maternedge(varargin)
%MATERNEDGE  Long-time cutoff edge for the Matern impulse response function.
%
%   TE=MATERNEDGE(ALPHA,H,EPSILON) determines a long-time cutoff for the
%   Green's function of a Matern process characterized by slope parameter 
%   ALPHA and range parameter H.
%
%   TE is the (approximately) the first time at which the ratio of the 
%   time-integrated magnitude of the Green's function, to the value it 
%   obtains at infinity, has risen to within EPSILON of unity.
%
%   TE is approximate because of discretization of time.  It will always be
%   greater than the actual time at which the EPSILON level is reached.  
%
%   See Sykulski, Lilly, Olhede, and Danioux (2013) for details.
%   
%   MATERNEDGE is a low-level function that is called by MATERNOISE.
%
%   See also MATERNSPEC, MATERNCOV, MATERNIMP, MATERNOISE.
%
%   'maternedge --t' runs some tests.
%
%   Usage: te=maternedge(alpha,h,epsilon);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    maternedge_test,return
end

alpha=varargin{1};
h=varargin{2};
epsilon=varargin{3};

arrayify(alpha,h,epsilon);
teps=zeros(size(alpha));

tnorm=[0:0.001:20]';
%tnorm=[0:0.01:20]';
for i=1:length(alpha)
    
    gint=gammainc(tnorm,alpha(i));  %This is integral divided by its total
    %Note Matlab defines the incomplete gamma function as a ratio to gamma
    ii=find(gint>1-epsilon(i),1,'first');
    %plot(tnorm,gint),hlines(1),hold on
    if ~isempty(ii)
        teps(i)=tnorm(ii)./h(i);
    else
        %disp('The Green''s function does not cross this threshold within 20 e-folding times.')
        teps(i)=inf;
    end 
end

if aresame(size(varargin{1}),size(varargin{2}))
    teps=reshape(teps,size(varargin{1},1),size(varargin{1},2));
end


function[]=maternedge_test

alpha=[1 1.5 2 3 4];
%alpha=[3/4 1 1.5 2 3 4];
h=[.01 .02 .05 .2 1];
[alpha,h]=meshgrid(alpha,h);
h=h.*alpha;
epsilon=1e-2;

vcolon(alpha,h);
te=maternedge(alpha,h,epsilon);

t=[1e-6:0.0001:20]';
tnorm=oprod(t,1./h);
g=zeros(size(tnorm));
gint=zeros(size(tnorm));
for i=1:length(alpha)
    d=sqrt(frac(h(i).^(2*alpha(i)-1),materncfun(alpha(i))));
    g(:,i)=frac(1,d)*maternimp(tnorm(:,i),alpha(i),h(i));
    gint(:,i)=cumsum(g(:,i)).*(tnorm(2,i)-tnorm(1,i));
    gint(:,i)=gint(:,i).*h(i).^alpha(i);
end

te2=zeros(size(te));
for i=1:length(alpha)
     ii=find(gint(:,i)>1-epsilon,1,'first');
     if ~isempty(ii)
         te2(i)=tnorm(ii,i);
     else
         te2(i)=nan;
     end
end

%abs(te2-te)./te

reporttest('MATERNEDGE matches explicit calculation to within one percent',allall(abs(te2-te)./te<.01))
