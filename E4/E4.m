clc, close all
data1 = importdata('x_data.dat');
data2 = importdata('v_data.dat');

n=5; % Nbr of trajectorys
t=data1(:,1);
mean_x=data1(:,n+2);
mean_v=data2(:,n+2);
var_x=data1(:,n+3);
var_v=data2(:,n+3);
std_x=sqrt(var_x);
std_v=sqrt(var_v);

figure(1)
hold on
for i = 1:n
    plot(t,data1(:,i+1),'-')
end
plot(t,mean_x,'black-')
plot(t,mean_x+std_x,'cyan-')
plot(t,mean_x-std_x,'cyan-')
figure(2)
hold on
for i = 1:n
    plot(t,data2(:,i+1),'-')    
end
plot(t,mean_v,'black-')
plot(t,mean_v+std_v,'cyan-')
plot(t,mean_v-std_v,'cyan-')
%%
clc
close all
clear all
x_dist = importdata('x_dist.dat');
v_dist = importdata('v_dist.dat');
bins = 25; 
figure(1)
hold on 
histogram(x_dist(6,:), bins,'Normalization','probability' )
histogram(x_dist(5,:), bins,'Normalization','probability' )
histogram(x_dist(4,:), bins, 'Normalization','probability')
histogram(x_dist(3,:), bins ,'Normalization','probability' )
histogram(x_dist(2,:), bins,'Normalization','probability' )
histogram(x_dist(1,:), 100, 'Normalization','probability')
ylim([0 0.13])
xlim([-0.1 0.2])
legend('t=1 ms','t=0.3 ms','t=0.2 ms','t=0.15 ms','t=0.05 ms','t=0 ms')

figure(2)
hold on 
histogram(v_dist(6,:), bins,'Normalization','probability' )
histogram(v_dist(5,:), bins,'Normalization','probability' )
histogram(v_dist(4,:), bins, 'Normalization','probability')
histogram(v_dist(3,:), bins ,'Normalization','probability' )
histogram(v_dist(2,:), bins,'Normalization','probability' )
histogram(v_dist(1,:), 20, 'Normalization','probability')
ylim([0 0.13])
xlim([-2 2.1])
legend('t=1 ms','t=0.3 ms','t=0.2 ms','t=0.15 ms','t=0.05 ms','t=0 ms')
%%
