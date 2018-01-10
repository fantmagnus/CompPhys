%% E
clc, %close all
data = importdata('E.dat');
E = data(:,1);
tau = data(:,2); 
N = 20000;




ylabel('Energy [$E_H$]','interpreter', 'latex','fontsize',16)
hold on
plot(tau,E,'linewidth',2)
plot(tau,cumsum(E)./(1:N)','r','linewidth',2)
plot([min(tau) max(tau)],[0.5 0.5], 'b-')
xlabel('Time [a.u.]','interpreter', 'latex','fontsize',16)
%title('$\Delta\tau = 0.01$ a.u., $\alpha = 5\cdot10^{-5}$','interpreter', 'latex','fontsize',16)
%legend({'$E_T$','Cumulative average of $E_T$'},'interpreter','latex','fontsize',13)
%ylim([-1 2])

E0 = mean(E)
Var_E = var(E) 

%% Independant runs 
% 2e4 steps
clc
E = [0.4943 0.5015 0.5078 0.5017 0.5064 0.5010 0.5084 0.5034 0.4948 0.5140 0.4944 0.4918  0.4990 0.4982 0.5006];
Var_E = [0.2220 0.1452 0.2228 0.1713 0.1073 0.1785 0.1124 0.0811 0.4123 0.2398 0.0507 0.0363 0.1027 0.1095 0.7343];
n = length(E);
mean_E = mean(E)
% 4e4 steps
E = [0.5000 0.5032 0.5017 0.5088 0.5029 0.5021 0.5016 0.4990 0.5014 0.5090];
Var_E = [0.6081 0.0330 0.1858 1.3942 0.3313 0.0713 0.1232 0.3517 0.1579 0.7937];
n = length(E);
mean_E = mean(E)

%% alpha
clc, close all
data1 = importdata('Ealpha2.dat');
E1 = data1(:,1);
data1 = importdata('Ealpha3.dat');
E2 = data1(:,1);
data1 = importdata('Ealpha4.dat');
E3 = data1(:,1);
data1 = importdata('Ealpha5.dat');
E4 = data1(:,1);
tau = data1(:,2); 
N = 10000;

hold on
%plot(tau,E1)
plot(tau,E2,'linewidth',2)
plot(tau,E3,'linewidth',2)
plot(tau,E4,'linewidth',2)
plot([min(tau) max(tau)],[0.5 0.5], 'b-')
legend({'$\alpha = 10^{-3}$','$\alpha = 10^{-4}$','$\alpha = 10^{-5}$'},'interpreter','latex')
xlabel('Time [a.u.]','interpreter', 'latex','fontsize',16)
ylabel('Energy [$E_H$]','interpreter', 'latex','fontsize',16)
title('$\Delta\tau = 0.01$ a.u.','interpreter', 'latex','fontsize',16)
ylim([-1 2])
xlim([0 50])
[mean(E2) mean(E3) mean(E4)]
[var(E2) var(E3) var(E4)]
%%
clc, clf 
data1 = importdata('M1.dat'); 
M1 = data1(:,1); 
tau = data1(:,2); 
data1 = importdata('M2.dat'); 
M2 = data1(:,1);
data1 = importdata('M3.dat'); 
M3 = data1(:,1);

subplot(3,1,1)
plot(tau, M1)
title('$\Delta\tau = 0.01$ a.u.','interpreter', 'latex','fontsize',16)

subplot(3,1,2)
plot(tau, M2)
ylabel('Number of walkers','interpreter', 'latex','fontsize',16)

subplot(3,1,3)
plot(tau, M3)
xlabel('Time [a.u.]','interpreter', 'latex','fontsize',16)

mean(M1)
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
bar(x,bin, 'y')
hold on
xx = linspace(-5,5,10000);
plot(xx,(1/sqrt(2*pi)*exp(-xx.^2/2)).^2,'r-')