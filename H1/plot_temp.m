% plot the energies
% Created by Martin Gren 2014-10-25.

% load the data file
data = importdata('temp.dat');
close all
%plot
figure;
plot(data(:,1),data(:,2),'-');

% labels
xlabel('Time / [dim. unit]');
ylabel('Energy / [dim. unit]');

% legend
legend('Temperature')

% axis limits
%xlim([0,0.1]);
%ylim([0,0.55]);
