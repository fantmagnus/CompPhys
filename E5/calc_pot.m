function V = calc_pot(phi,r)
h = r(2) - r(1);
N=length(r);
n = abs(phi).^2;
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
   b(i) = -4*pi*r(i)*n(i);
end
b(N) = b(N) - UN/h^2;

U = A\b;
V = U./r;
end