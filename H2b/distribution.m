%% Dist of r
clc, close all
data = importdata('dist.dat');


Z_cf = 2;       % Central field approximation
Z_vo = 27/16;   % Variationally optimized
r=linspace(0,6,500);
rho_cf=Z_cf^3*4*r.^2.*exp(-2*Z_cf*r);   % PDFs
rho_vo=Z_vo^3*4*r.^2.*exp(-2*Z_vo*r);

figure(1)
clf
hold on
histogram(data,100,'normalization','pdf','faceColor','y')
plot(r,rho_cf,'b')
plot(r,rho_vo,'r')


xlabel('Distance to nucleus [a.u.]','interpreter','LaTeX','fontSize',16)
ylabel('Probability density','interpreter','LaTeX','fontSize',16)
legend('Sampled data','Z=2','Z=27/16')

%% Dist of x
clc, close all
data = importdata('corr.dat');

figure(2)
hold on
[h,x] = hist(data,800);
h=h/trapz(x,h);
bar(x,h)


xlabel('$x=\cos\theta$','interpreter','LaTeX','fontSize',16)
ylabel('Probability density','interpreter','LaTeX','fontSize',16)

