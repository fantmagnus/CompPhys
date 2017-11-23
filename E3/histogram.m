% load the data file
data = importdata('generated_points.dat');
clc, close all
%plot
figure(1);
hold on
hist(data,30)
s=sqrt(7e-5)

% title
title('N = 10 000, \sigma = 0.0084')

% labels
%xlabel('Time [(m/\kappa)^{1/2}]');
%ylabel('Energy [\kappa a_0^2]');
xlabel('\eta');
text()


% axis limits
%xlim([0,10]);
%ylim([0,33]);