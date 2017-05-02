function[]=makefigs_stickvect
%MAKEFIGS_STICKVECT   Makes a sample figure for STICKVECT.

load bravo94
cv=vfilt(bravo94.rcm.cv,100);
figure
stickvect(yearfrac(bravo94.rcm.num),180,cv,300,-30);