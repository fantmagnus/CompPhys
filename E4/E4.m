clc, close all
x = importdata('x_data.dat');
v = importdata('v_data.dat');
tA=[0 0.05 0.15 0.2 0.3 1];
tB=[0 0.075 0.3 0.45 0.95 2.5];

n=5; % Nbr of trajectorys
t=x(:,1);
mean_x=x(:,n+2);
mean_v=v(:,n+2);
var_x=x(:,n+3);
var_v=v(:,n+3);
std_x=sqrt(var_x);
std_v=sqrt(var_v);

index_B=zeros(1,6);
index_A=zeros(1,6);
for i=1:6
   [a,I] = min(abs(t-tA(i)));
   index_B(i) = I;   
   [a,I] = min(abs(t-tA(i)));
   index_A(i) = I;
end


figure(1)
hold on
for i = 1:n
    plot(t,x(:,i+1),'-')
end
plot(t,mean_x,'black-')
plot(t,mean_x+std_x,'cyan-')
plot(t,mean_x-std_x,'cyan-')
plot(tB,mean_x(index_B),'b*')
xlabel('t [ms]')
ylabel('x [\mum]')
xlim([0 3])

figure(2)
hold on
for i = 1:n
    plot(t,v(:,i+1),'-')    
end
plot(t,mean_v,'black-')
plot(t,mean_v+std_v,'cyan-')
plot(t,mean_v-std_v,'cyan-')
plot(tB,mean_v(index_B),'b*')
xlabel('t [ms]')
ylabel('v [\mum/ms]')
xlim([0 3])
%%
clc
clear all
figure(3)
clf
figure(4)
clf
x_dist = importdata('x_dist.dat');
v_dist = importdata('v_dist.dat');
bins = 40;
N = 6;

figure(3)
hold on 
histogram(x_dist(1,:), 100, 'Normalization','pdf')
for i=2:N
    histogram(x_dist(i,:), bins,'Normalization','pdf' ) 
end
ylim([0 37])
xlim([-0.07 0.15])
legend('t=0 ms','t=0.05 ms','t=0.15 ms','t=0.2 ms','t=0.3 ms','t=9 ms')
%legend('t=0 ms','t=0.075 ms','t=0.3 ms','t=0.45 ms','t=0.9 ms','t=9 ms')
xlabel('x [\mum]')
title('p_x')


figure(4)
hold on 
histogram(v_dist(1,:), 10, 'Normalization','pdf')
for i=2:N
    histogram(v_dist(i,:), bins,'Normalization','pdf' )
end
ylim([0 2])
xlim([-2 2.1])
legend('t=0 ms','t=0.05 ms','t=0.15 ms','t=0.2 ms','t=0.3 ms','t=9 ms')
%legend('t=0 ms','t=0.075 ms','t=0.3 ms','t=0.45 ms','t=0.9 ms','t=9 ms')
xlabel('v [\mum/ms]')
title('p_v')
