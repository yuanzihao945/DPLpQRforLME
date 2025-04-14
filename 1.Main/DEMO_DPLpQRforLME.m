clearvars
clc
N = 300;
n = 4;
k = 10;
P = 6;
s = ceil(n * rand(N, 1));
x = randn(N, k);
z = randn(N, P);
beta0 = [randn(round(k/2), 1); zeros(k-round(k/2), 1)];
alpha0 = [randn(n, round(P/2)), zeros(n, P-round(P/2))];
tau = 0.5;
y = zeros(N, 1);
for i = 1 : N
    y(i) = 1.5 + x(i, :) * beta0 + z(i, :) * alpha0(s(i), :)' + randn(1);
end
p = 1.0;
[beta, alpha, lambda_beta, lambda_alpha] = DPLpQRforLME(s, y, x, z, p, tau);