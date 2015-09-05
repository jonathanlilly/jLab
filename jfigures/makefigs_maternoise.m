function[]=makefigs_maternoise
%MAKEFIGS_MATERNOISE  Makes a sample figure for MATERNOISE.

N=1000; 
alpha=[0.6 1 1.5 2 3 4];
h=[.01 .02 .05 .2 1];
[alpha,h]=meshgrid(alpha,h);
h=h.*alpha;

rng(1);  %set seed
z=maternoise(1,N,10,alpha,h);
z=z./vrep(std(z,1,1),size(z,1),1);  %Set to unit std    
y=detrend(cumsum(z),'constant');
y=y./vrep(std(y,1,1),size(y,1),1);  %Set to unit std    

%This is just to make an offset
[xo,yo]=meshgrid(1:6,1:5);
zo=xo(:)'+1i*yo(:)';
zo=vrep(zo,length(z),1);

figure,
plot(y+zo*3),axis equal,
xlabel('Increasing \alpha \rightarrow'),ylabel('Increasing h \rightarrow')
set(gca,'xticklabel',[]),set(gca,'xticklabel',[])
title('Example of Matern Processes')
