function[]=makefigs_stickvect
%MAKEFIGS_STICKVECT   Makes a sample figure for STICKVECT.

load bravo94
cv=vfilt(bravo.rcm.cv,100);
figure
stickvect(bravo.rcm.yearf,180,cv,300,-30);