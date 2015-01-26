alpha = 0.05
L = 100
lambda1 = 2.8889
tau = 10

syms x y

sol = double(solve(x + alpha * exp(-L * (x + alpha) / (tau * lambda1 * alpha))))