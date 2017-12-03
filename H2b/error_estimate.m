clc, close all
data = importdata('phi.dat');
k=data(:,1);
phi=data(:,2);

hold on
plot(k,phi,'-')
plot([min(k) max(k)],[exp(-2) exp(-2)],'-r')
xlim([0 150])


s = 2*sum(phi) % statistical inefficiency