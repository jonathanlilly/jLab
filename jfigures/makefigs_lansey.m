function[]=makefigs_lansey
%MAKEFIGS_LANSEY  Makes a sample figure for LANSEY.

figure
subplot(1,3,1),contourf(peaks,30),colormap(gca,'lansey'),xlabel('lansey')
subplot(1,3,2),contourf(peaks,30),colormap(gca,'parula'),xlabel('parula')
subplot(1,3,3),contourf(peaks,30),colormap(gca,'jet'),xlabel('jet')
for i=1:3,subplot(1,3,i),axis equal,axis tight,caxis([-7 7]),...
        noxlabels,noylabels,nocontours,xtick([-100 100]),ytick([-100 100])
end
hc=packfig(1,3);
for i=1:3,axes(hc(i)),colorbar('southoutside'),end
set(gcf,'paperposition',[1 1 12 6])



% figure
% subplot(2,2,1),contourf(peaks(200),60),colormap(gca,'lansey'),%xlabel('lansey')
% subplot(2,2,2),contourf(peaks(200),60),colormap(gca,'parula'),%xlabel('parula')
% subplot(2,2,3),contourf(peaks(200),60),colormap(gca,'jet'),%xlabel('jet')
% subplot(2,2,4),contourf(peaks(200),60),colormap(gca,'seminfhaxby'),%xlabel('seminf-haxby')
% for i=1:4,subplot(2,2,i),axis equal,axis tight,caxis([-7 7]),...
%         noxlabels,noylabels,nocontours,ylim([-50 200])
% end
% h=packfig(2,2,'both');
% str={'lansey','parula','jet','seminf-haxby'};
% for i=1:4,axes(h(i))
%     hc=colorbar('south'),
%     %hc.AxisLocation='out';
%     hc.Label.String=str{i};
% end
% set(gcf,'paperposition',[1 1 8 9.5 ])


%To print
if 0
    currentdir=pwd;
    cd([whichdir('jlab_license') '/figures'])
    print -dpng lansey
    crop lansey.png
    cd(currentdir)
end

