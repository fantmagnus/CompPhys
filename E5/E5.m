%% Task 1
close all, clc

% r
r_min = 0;
r_max = 10;
N = 100;
r = linspace(r_min,r_max,N+2);
r = r(2:end-1)';
h = r(2) - r(1);
Z=1;
n=@(r) 1/(4*pi)*Z^3*4*exp(-2*Z*r);
UN = 1;

% A
A = zeros(N,N);
for i = 2:N-1
   A(i,i) = -2/h^2;
   A(i,i-1) = 1/h^2;
   A(i,i+1) = 1/h^2;   
end
A(1,1) = -2/h^2;
A(N,N) = -2/h^2;
A(N,N-1) = 1/h^2;
A(1,2) = 1/h^2;

% b
b = zeros(N,1);
for i = 1:N
   b(i) = -4*pi*r(i)*n(r(i));
end
b(N) = b(N) - UN/h^2;

U = A\b;
plot(r,U./r,'cyan')
hold on
V_H = 1./r - (1+1./r).*exp(-2*r);
plot(r,V_H,'--r')
xlabel('r [a.u]')
ylabel('E [a.u]')
legend('FD solution', 'Hartree potential')

%% Task 2
% r
close all,clc
r_min = 0;
r_max = 10;
Z = 1;
N = 1000;
r = linspace(r_min,r_max,N+2);
r = r(2:end-1)';
h = r(2) - r(1);

% A
A = zeros(N,N);
for i = 2:N-1
   A(i,i) = 1/h^2 - 1/r(i);
   A(i,i-1) = -1/(2*h^2);
   A(i,i+1) = -1/(2*h^2);   
end
A(1,1) = 1/h^2 - 1/r(1);
A(N,N) = 1/h^2 - 1/r(N);
A(N,N-1) = -1/(2*h^2);
A(1,2) = -1/(2*h^2);

% epsilon
[F,lambda] = eig(A);
epsilon = lambda(1,1)
f = -F(:,1)/sqrt(trapz(r,F(:,1).^2));
R_0 = Z*2*h*f(1)/r(1)+f(2)/r(2);
plot([0;r],[R_0/sqrt(4*pi);f./(sqrt(4*pi)*r)])
hold on
f_anal = 1/sqrt(pi)*exp(-r);
plot(r,f_anal,'--r')
xlabel('r [a.u]')
ylabel('E [a.u]')
legend('FD solution', 'Analytic solution')

%% Task 3
clc, close all
% r
r_min = 0;
r_max = 10;
N = 1000;
r = linspace(r_min,r_max,N+2);
r = r(2:end-1)';
Z = 2;
tol = 1e-6;
% Initial wave-function
phi = 1/sqrt(4*pi)*Z^(3/2)*2*r.*exp(-Z*r);
%phi = phi/sqrt(trapz(r,phi.^2));
E1 = 0;
E2 = 3;

i = 0
while abs(E2-E1) > tol
    E1 = E2;
    VH = calc_pot(phi,r);
    
    A = hamiltonian(VH,r);    
    
    [F,lambda] = eig(A);
    
    f = F(:,1);
    phi = f./(sqrt(4*pi)*r);
    phi = phi / sqrt(trapz(r,4*pi*r.^2.*phi.^2));
    E2 = 2*lambda(1,1) - trapz(r,4*pi*r.^2.*VH.*abs(phi).^2);
    
   abs(E2-E1)
   i = i+1;
end
i
E = E2

clf,
n = abs(phi).^2;
rho = 4*pi*r.^2.*n;
plot(r,rho/trapz(r,rho),'b.')
hold on
Z = 2;
rho_cf = 4*r.^2*Z^3.*exp(-2*Z*r);
plot(r,rho_cf,'black')
Z = 27/16;
rho_cf = 4*r.^2*Z^3.*exp(-2*Z*r);
plot(r,rho_cf,'red')
xlim([0 5])