% General parameters of the problem
tol = 1e-10;
itmax = 100;

% General loop to compute for all mus

mu = 0:0.01:0.5;
n = size(mu, 2); %% n = length(mu) ???

L = zeros(3, n);

for i = 1:n
    s0 = ( mu(i)/(3*(1 - mu(i))) )^1/3; % Falta raiz cubica
    [xk, res, it] = fixed_pont_iter(s0, tol, itmax, @(x) G1(x, mu(i)));
    L(1, i) = -(1 - mu(i)) + xk(end);
    it
end

plot(mu, L(1,1:n));


function g1 = G1(s, mu)
    g1 = ( (mu*(1 + s)^2)/(3 - 2*mu + s*(3 - mu + s)) )^1/3;
end