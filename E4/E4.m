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