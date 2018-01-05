%% E
clc, close all
data1 = importdata('E.dat');
E = data1(:,1);

hold on
plot(E,'r')
mean(E)
%%
clc, clf 
data2 = importdata('M.dat'); 
M = data2(:,1); 
tau = data2(:,2); 
plot(tau, M)
mean(M)