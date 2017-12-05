%% Dist of r
clc, close all
data = importdata('dist.dat');


Z_cf = 2;       % Central field approximation
Z_vo = 27/16;   % Variationally optimized
r=linspace(0,6,500);
rho_cf=Z_cf^3*4*r.^2.*exp(-2*Z_cf*r);   % PDFs
rho_vo=Z_vo^3*4*r.^2.*exp(-2*Z_vo*r);

figure(1)
hold on
plot(r,rho_cf,'b')
plot(r,rho_vo,'r')
[h,x] = hist(data,1000);
bar(x,h/max(h))

xlabel('Distance to nucleus [a.u.]','interpreter','LaTeX')
legend('Sampled data','Z=2','Z=27/16')

trapz(x,h/max(h)) % Should be close to one