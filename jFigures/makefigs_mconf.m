function[]=makefigs_mconf
%MAKEFIGS_MCONF  Makes sample figure for MCONF.

load bravo94
use bravo94
use bravo94.rcm

cv=cv(:,5);

figure
for i=1:2
    subplot(1,2,i)
    P=[3 8];
    psi=sleptap(length(cv),P(i));
    [f,spp,snn,spn]=mspec(cv,conj(cv),psi);
    h=plot(f,spp);linestyle -h h 2T
    [ra,rb]=mconf(2*P(i)-1,0.95);
    hold on,h1=plot(f,spp.*ra,'k');plot(f,spp.*rb,'k')
    [ra,rb]=mconf(2*P(i)-1,0.95,'log10');
    hold on,h2=plot(f,10.^(log10(spp)+ra),'r');plot(f,10.^(log10(spp)+rb),'r')
    ylog,axis tight,xlim([0 1.5])
    title(['P = ' int2str(P(i)) ' power spectrum with 95\% confidence intervals'])
    if i==2
        legend([h h1 h2],'Spectral estimate','Linear intervals','Log10 intervals')
    end
end