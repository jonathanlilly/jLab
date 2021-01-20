function[]=makefigs_seminfhaxby
%MAKEFIGS_SEMINFHAXBY  Makes a sample figure for SEMINFHAXBY.

figure
subplot(2,2,1),contourf(abs(peaks(200)),60),colormap(gca,'seminfhaxby'),%xlabel('seminf-haxby')
subplot(2,2,2),contourf(abs(peaks(200)),60),colormap(gca,'lansey'),%xlabel('lansey')
subplot(2,2,3),contourf(abs(peaks(200)),60),colormap(gca,'parula'),%xlabel('parula')
subplot(2,2,4),contourf(abs(peaks(200)),60),colormap(gca,'jet'),%xlabel('jet')
for i=1:4,subplot(2,2,i),axis equal,axis tight,caxis([0 8.2]),...
        noxlabels,noylabels,nocontours,ylim([-50 200])
end
h=packfig(2,2,'both');
str={'seminf-haxby','lansey','parula','jet'};
for i=1:4,axes(h(i))
    hc=colorbar('south'),
    %hc.AxisLocation='out';
    hc.Label.String=str{i};
end
set(gcf,'paperposition',[1 1 8 9.5 ])

%To print
if 0
    currentdir=pwd;
    cd([whichdir('jlab_license') '/figures'])
    print -dpng seminfhaxby
    crop seminfhaxby.png
    cd(currentdir)
end

