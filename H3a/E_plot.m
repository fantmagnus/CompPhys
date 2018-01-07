%% E
clc, close all
data1 = importdata('E.dat');
E = data1(:,1);
tau = data1(:,2); 
N = 18000;
hold on
plot(tau,cumsum(E)./(1:N)','r')
mean(E)
var(E)
%%
clc, clf 
data2 = importdata('M.dat'); 
M = data2(:,1); 
tau = data2(:,2); 
plot(tau, M)
mean(M)
%%
close all
dist = importdata('dist.dat'); 
min(dist)
max(dist)
histogram(dist,20)
%%
close all
clc
data3 = importdata('bin.dat'); 
x = data3(:,1); 
bin = data3(:,2);
bar(x,bin)
hold on
xx = linspace(-5,5);
plot(xx,(1/sqrt(2*pi)*exp(-xx.^2/2)).^2,'r-')