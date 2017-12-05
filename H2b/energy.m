%% E(alpha)
clc, close all
data = importdata('E.dat');
alpha = data(:,1);
E = data(:,2);
var_E = 20*data(:,3);
s=sqrt(var_E);

errorbar(alpha,E,var_E,'r.','vertical')
%xlim([0.125 0.175])
xlabel('$\alpha$','interpreter','LaTeX')
ylabel('$E$ [a.u.]','interpreter','LaTeX')

%% E_L
clc, close all
data = importdata('E.dat');
E = data(1:20000,1);
plot(E)
xlabel('Iteration','interpreter','LaTeX')
ylabel('$E$ [a.u.]','interpreter','LaTeX')

%% E independant runs.
runs = 10;
N=1e7;
n_s=12;
N_tot=N*runs;
E=[-2.877649 -2.878007 -2.878397 -2.877272 -2.878002 -2.878062 -2.877571 -2.877951 -2.879088 -2.877711];
var_E=[0.113638 0.113446 0.113689 0.113625 0.113621 0.113856 0.113784 0.113645 0.113759 0.113587];
E_mean = mean(E)
E_variance = mean(var_E)*n_s/N_tot
E_std = sqrt(E_variance)

