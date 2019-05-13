function[]=makefigs_jfig
%MAKEFIGS_JFIG  Makes a sample figure for JFIG.

figure, jpcolor(peaks)
jfig equal ticksout yrev axis|[5 45 5 45] caxis|[-8 8] ...
          title|['Demonstration of JFIG'] colorbar|['eastoutside']