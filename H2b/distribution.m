clc, close all
data = importdata('dist.dat');
[h,x] = hist(data,1000);
bar(x,h/max(h))
hold on
Z_cf = 2;
Z_vo = 27/16;
r=linspace(0,6,50);
rho_cf=Z_cf^3*4*r.^2.*exp(-2*Z_cf*r);
rho_vo=Z_vo^3*4*r.^2.*exp(-2*Z_vo*r);
figure(1)
hold on
plot(r,rho_cf,'b')
plot(r,rho_vo,'r')
xlabel('Distance to nucleus [a.u.]')
legend('Sampled data','Z=2','Z=27/16')
trapz(x,h/max(h)) % Normalized p? 1.0840 for 1000 bins, 1.0818 for 2000 bins, 1.0909 for 500 bins