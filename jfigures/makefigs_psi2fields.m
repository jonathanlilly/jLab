function[]=makefigs_psi2fields
%MAKEFIGS_PSI2FIELDS  Makes a sample figure for PSI2FIELDS.

load qgsnapshot,use qgsnapshot
psi=qgsnapshot.psi;

[cv,zeta,N,S,P]=psi2fields(x(2)-x(1),psi);
x=x/1000;
y=y/1000;
figure
subplot(2,2,1),jpcolor(x,y,psi),title('Streamfunction')
subplot(2,2,2),jpcolor(x,y,zeta),title('Vorticity')
subplot(2,2,3),jpcolor(x,y,P),title('Okubo-Weiss')
subplot(2,2,4),ii=1:10:length(x);
quiver(x(ii),y(ii),real(cv(ii,ii)),imag(cv(ii,ii))),title('Velocity')

for i=1:4
    subplot(2,2,i)
    axis equal, axis tight,xtick([-3:3]),ytick([-3:1:3])
    fontsize 16 14 14 14
end
