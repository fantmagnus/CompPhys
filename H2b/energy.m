%%
clc, close all
data = importdata('E.dat');
alpha = data(:,1);
E = data(:,2);
var_E = data(:,3);
s=sqrt(var_E);

errorbar(alpha,E,var_E,'r.','vertical')
%xlim([0.125 0.175])
xlabel('$\alpha$','interpreter','LaTeX')
ylabel('$E$ [a.u.]','interpreter','LaTeX')

%%
clc, close all
data = importdata('E.dat');
E = data(:,1);
var(E)