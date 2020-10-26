%% Parameters and preliminary declarations

% General parameters of the problem
tol = 1e-10;
itmax = 100;

% General loop to compute for all mus

mu = 0:0.01:0.5;
n = size(mu, 2); %% n = length(mu) ???

L1 = zeros(3, n); % L1 point for all mus and different methods
L2 = zeros(3, n);
L3 = zeros(3, n);

%% Loop for the L1 iterations
for i = 1:n
    % Fixed point iteration
    s0 = ( mu(i)/(3*(1 - mu(i))) )^(1/3);
    [xk, res, it] = iteracio_simple(s0, tol, itmax, @(x) G1(x, mu(i)));
    L1(1, i) = -(1 - mu(i)) - xk(end);
    
    % Newton method
    [xk, res, it] = newton(s0, tol, itmax, @(x) P1(x, mu(i)), @(x) dP1(x, mu(i)));
    L1(2, i) = -(1 - mu(i)) - xk(end);
    
    % Secant method
    [xk, res, it] = secant(0.3, 0.8, tol, itmax, @(x) P1(x, mu(i)));
    L1(3, i) = -(1 - mu(i)) - xk(end);
end

%% Loop for the L2 iterations
for i = 1:n
    % Fixed point iteration
    s0 = ( mu(i)/(3*(1 - mu(i))) )^(1/3);
    [xk, res, it] = iteracio_simple(s0, tol, itmax, @(x) G2(x, mu(i)));
    L2(1, i) = -(1 - mu(i)) + xk(end);
    
    % Newton method
    [xk, res, it] = newton(s0, tol, itmax, @(x) P2(x, mu(i)), @(x) dP2(x, mu(i)));
    L2(2, i) = -(1 - mu(i)) + xk(end);
    
    % Secant method
    [xk, res, it] = secant(0.3, 0.8, tol, itmax, @(x) P2(x, mu(i)));
    L2(3, i) = -(1 - mu(i)) + xk(end);
end

%% Loop for the L3 iterations
for i = 1:n
    % Fixed point iteration
    s0 = 1 - 7*mu(i)/12;
    [xk, res, it] = iteracio_simple(s0, tol, itmax, @(x) G3(x, mu(i)));
    L3(1, i) = mu(i) + xk(end);
    
    % Newton method
    [xk, res, it] = newton(s0, tol, itmax, @(x) P3(x, mu(i)), @(x) dP3(x, mu(i)));
    L3(2, i) = mu(i) + xk(end);
    
    % Secant method
    [xk, res, it] = secant(0.3, 0.8, tol, itmax, @(x) P3(x, mu(i)));
    L3(3, i) = mu(i) + xk(end);
end

%% Plots and visualization

close all

% Plots for the L1 point with each method
hold on
plot(mu, L1(1,:));
plot(mu, L1(2,:));
plot(mu, L1(3,:));
xlabel('\mu');
ylabel('x_{L1}');
legend({'Simple', 'Newton', 'Tangent'}, 'Location', 'Northeast')
title('Plot of the L1 point as a function of \mu for different methods');
figure;
hold off

% Plots for the L2 point with each method
hold on
plot(mu, L2(1,:));
plot(mu, L2(2,:));
plot(mu, L2(3,:));
xlabel('\mu');
ylabel('x_{L2}');
legend({'Simple', 'Newton', 'Tangent'}, 'Location', 'Northwest')
title('Plot of the L2 point as a function of \mu for different methods');
figure;
hold off

% Plots for the L3 point with each method
hold on
plot(mu, L3(1,:));
plot(mu, L3(2,:));
plot(mu, L3(3,:));
xlabel('\mu');
ylabel('x_{L3}');
legend({'Simple', 'Newton', 'Tangent'}, 'Location', 'Northwest')
title('Plot of the L3 point as a function of \mu for different methods');
figure;
hold off

% Plots for all points with a simple iteration procedure
hold on
plot(mu, L1(1,:));
plot(mu, L2(1,:));
plot(mu, L3(1,:));
xlabel('\mu');
ylabel('x');
legend({'Lagrange 1', 'Lagrange 2', 'Lagrange3'}, 'Location', 'Northeast')
title('Plot of the Lagrange points as a function of \mu by simple iteration');
figure;
hold off

% Plots for the Jacobi constant of equilibrium points for each mu
hold on
plot(mu, omega(L1(1,:), mu));
plot(mu, omega(L2(1,:), mu));
plot(mu, omega(L3(1,:), mu));
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

function pol1 = P1(s, mu)
    pol1 = s^5 + (3 - mu)*s^4 + (3 - 2*mu)*s^3 - mu*s^2 - 2*mu*s - mu;
end

function dpol1 = dP1(s, mu)
    dpol1 = 5*s^4 + 4*(3 - mu)*s^3 + 3*(3 - 2*mu)*s^2 - 2*mu*s - 2*mu;
end

% Functions for the L1 point
function g2 = G2(s, mu)
    g2 = ( (mu*(1 - s)^2)/(3 - 2*mu - s*(3 - mu - s)) )^(1/3);
end

function pol2 = P2(s, mu)
    pol2 = s^5 - (3 - mu)*s^4 + (3 - 2*mu)*s^3 - mu*s^2 + 2*mu*s - mu;
end

function dpol2 = dP2(s, mu)
    dpol2 = 5*s^4 - 4*(3 - mu)*s^3 + 3*(3 - 2*mu)*s^2 - 2*mu*s + 2*mu;
end

% Functions for the L1 point
function g3 = G3(s, mu)
    g3 = ( ( (1 - mu)*(1 + s)^2 )/(1 + 2*mu + s*(2 + mu + s)) )^(1/3);
end

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