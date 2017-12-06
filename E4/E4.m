data1 = importdata('x_data.dat');
data2 = importdata('v_data.dat');
n=5; % Nbr of trajectorys
figure(1)
hold all
for i = 1:n
    plot(data1(:,1),data1(:,i+1),'-')
    plot(data2(:,1),data2(:,i+1),'-')    
end