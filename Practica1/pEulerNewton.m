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
    [xk1, res1, it1] = newton(s0, tol, itmax, @(x) P1(x, mu(i)), @(x) dP1(x, mu(i)));
    L1(i) = -(1 - mu(i)) - xk1(end);
    
    % L2 point
    [xk2, res2, it2] = newton(s0, tol, itmax, @(x) P2(x, mu(i)), @(x) dP2(x, mu(i)));
    L2(i) = -(1 - mu(i)) + xk2(end);

    % L3 point
    s0p = 1 - 7*mu(i)/12;
    [xk3, res3, it3] = newton(s0p, tol, itmax, @(x) P3(x, mu(i)), @(x) dP3(x, mu(i)));
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
title('Plot of the Lagrange points as a function of \mu by Newton');
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
title('Plot of the Jacobi constant as a function of \mu by Newton');
hold off

%% Auxiliary functions for the PCRTBP

% Functions for the L1 point
function pol1 = P1(s, mu)
    pol1 = s^5 + (3 - mu)*s^4 + (3 - 2*mu)*s^3 - mu*s^2 - 2*mu*s - mu;
end

function dpol1 = dP1(s, mu)
    dpol1 = 5*s^4 + 4*(3 - mu)*s^3 + 3*(3 - 2*mu)*s^2 - 2*mu*s - 2*mu;
end

% Functions for the L2 point
function pol2 = P2(s, mu)
    pol2 = s^5 - (3 - mu)*s^4 + (3 - 2*mu)*s^3 - mu*s^2 + 2*mu*s - mu;
end

function dpol2 = dP2(s, mu)
    dpol2 = 5*s^4 - 4*(3 - mu)*s^3 + 3*(3 - 2*mu)*s^2 - 2*mu*s + 2*mu;
end

% Functions for the L3 point
function pol3 = P3(s, mu)
    pol3 = s^5 + (2 + mu)*s^4 + (1 + 2*mu)*s^3 - (1 - mu)*s^2 - 2*(1 - mu)*s - (1 - mu);
end

function dpol3 = dP3(s, mu)
    dpol3 = 5*s^4 + 4*(2 + mu)*s^3 + 3*(1 + 2*mu)*s^2 - 2*(1 - mu)*s - 2*(1 - mu);
end

% Function for the omega plot
function res = omega(x, mu)
    res = x.*x/2 + (1 - mu)./abs(x - mu) + mu./abs(x - mu + 1) + mu.*(1 - mu)/2;
end