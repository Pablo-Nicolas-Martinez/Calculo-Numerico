%% Parameters and preliminary declarations
clear all
close all

% General parameters of the problem
tol = 1e-10;
itmax = 100;

% General loop to compute for all mus

mu = 0:0.0001:0.5;
n = length(mu);

L1 = zeros(1, n);
L2 = zeros(1, n);
L3 = zeros(1, n);

%% Loop for the fixed point iterations
for i = 1:n
    % L1 point
    s0 = ( mu(i)/(3*(1 - mu(i))) )^(1/3);
    [xk1, res1, it1] = iteracio_simple(s0, tol, itmax, @(x) G1(x, mu(i)));
    L1(i) = -(1 - mu(i)) - xk1(end);
    
    % L2 point
    [xk2, res2, it2] = iteracio_simple(s0, tol, itmax, @(x) G2(x, mu(i)));
    L2(i) = -(1 - mu(i)) + xk2(end);

    % L3 point
    s0p = 1 - 7*mu(i)/12;
    [xk3, res3, it3] = iteracio_simple(s0p, tol, itmax, @(x) G3(x, mu(i)));
    L3(i) = mu(i) + xk3(end);
end

%% Plots and visualization

close all

% Plots for the Lagrange points with simple iteration
hold on
plot(mu, L1);
plot(mu, L2);
plot(mu, L3);
xlabel('\mu');
ylabel('x');
legend({'Lagrange 1', 'Lagrange 2', 'Lagrange 2'}, 'Location', 'Northeast')
title('Plot of the Lagrange points as a function of \mu by simple iteration');
figure;
hold off

% Plots for the Jacobi constant of equilibrium points for each mu
hold on
plot(mu, omega(L1, mu));
plot(mu, omega(L2, mu));
plot(mu, omega(L3, mu));
xlabel('\mu');
ylabel('C');
legend({'Lagrange 1', 'Lagrange 2', 'Lagrange 3'}, 'Location', 'Northwest')
title('Plot of the Jacobi constant as a function of \mu by simple iteration');
hold off

%% Auxiliary functions for the PCRTBP

% Functions for the L1 point
function g1 = G1(s, mu)
    g1 = ( (mu*(1 + s)^2)/(3 - 2*mu + s*(3 - mu + s)) )^(1/3);
end

% Functions for the L1 point
function g2 = G2(s, mu)
    g2 = ( (mu*(1 - s)^2)/(3 - 2*mu - s*(3 - mu - s)) )^(1/3);
end

% Functions for the L1 point
function g3 = G3(s, mu)
    g3 = ( ( (1 - mu)*(1 + s)^2 )/(1 + 2*mu + s*(2 + mu + s)) )^(1/3);
end

% Function for the omega plot
function res = omega(x, mu)
    res = x.*x/2 + (1 - mu)./abs(x - mu) + mu./abs(x - mu + 1) + mu.*(1 - mu)/2;
end