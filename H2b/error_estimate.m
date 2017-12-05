%% Correlation function
clc, close all, clear all
data = importdata('phi.dat');
k=data(:,1);
phi=data(:,2);

figure(1)
hold on
plot(k,phi,'-')
plot([min(k) max(k)],[exp(-2) exp(-2)],'-r')
xlim([0 150])

xlabel('k','interpreter','LaTeX')
ylabel('Correlation function $\Phi_k$','interpreter','LaTeX')

dim = [.5 .01 .5 .5];
str = 'exp(-2)';
%annotation('textbox',dim,'String',str,'FitBoxToText','on');


s = 2*sum(phi) % statistical inefficiency


% s = 10.28 when using phi_k=s=exp

%% Block method
data = importdata('S.dat');
B = data(:,1);
S = data(:,2);
Blim=150;
figure(2)
plot(B(1:Blim),S(1:Blim),'.')
xlabel('$B$','interpreter','LaTeX')
ylabel('Statistical inefficiency $n_s$','interpreter','LaTeX')
