clear all
close all

%% Parameters for the problem
mu = 0.0121529;
C = [3.01, 3.0242, 3.08, 3.19, 3.3];

%% Stable points of the Lagrange solution
s0 = ( mu/(3*(1 - mu)) )^(1/3);
s0p = 1 - 7*mu/12;

[L1, ~, ~] = iteracio_simple(s0, 1e-10, 100, @(x) G1(x, mu));
[L2, ~, ~] = iteracio_simple(s0, 1e-10, 100, @(x) G2(x, mu));
[L3, ~, ~] = iteracio_simple(s0p, 1e-10, 100, @(x) G3(x, mu));

l1 = -(1 - mu) - L1(end);
l2 = -(1 - mu) + L2(end);
l3 = mu + L3(end);

%% Orbits for each Jacobi constant
X1 = contzvc([-1/2; sqrt(3)/2], 0.001, 13, 1, @(x) JacobDiff(x, mu, C(1)), @(x) DJacobDiff(x, mu));
X2 = contzvc([-1/2; sqrt(3)/2], 0.001, 13, 1, @(x) JacobDiff(x, mu, C(2)), @(x) DJacobDiff(x, mu));
X3 = contzvc([-1/2; sqrt(3)/2], 0.001, 13, 1, @(x) JacobDiff(x, mu, C(3)), @(x) DJacobDiff(x, mu));
X4o = contzvc([-1/2; sqrt(3)/2], 0.001, 13, 1, @(x) JacobDiff(x, mu, C(4)), @(x) DJacobDiff(x, mu));
X4i = contzvc([0; 0.3], 0.001, 13, 1, @(x) JacobDiff(x, mu, C(4)), @(x) DJacobDiff(x, mu));
X5o = contzvc([-1/2; sqrt(3)/2], 0.001, 13, 1, @(x) JacobDiff(x, mu, C(5)), @(x) DJacobDiff(x, mu));
X5i = contzvc([0; 0.3], 0.001, 13, 1, @(x) JacobDiff(x, mu, C(5)), @(x) DJacobDiff(x, mu));

%% Plots for the curves found

hold on

% Plot of the various curves
plot(X1(1,:), X1(2,:), 'm')
plot(X2(1,:), X2(2,:), 'b')
plot(X3(1,:), X3(2,:), 'g')
plot(X4o(1,:), X4o(2,:), 'r')
plot(X5o(1,:), X5o(2,:), 'k')

plot(X1(1,:), -X1(2,:), 'm')
plot(X4i(1,:), X4i(2,:), 'r')
plot(X5i(1,:), X5i(2,:), 'k')

% Plot of the complementary lines
plot([-1, -1/2, 0], [0, sqrt(3)/2, 0], '--k');
plot([-1, -1/2, 0], [0, -sqrt(3)/2, 0], '--k');

% Plot of the different points
fill(l1 + 0.01*cos(0:0.1:2*pi), 0.01*sin(0:0.1:2*pi), 'r');
fill(l2 + 0.01*cos(0:0.1:2*pi), 0.01*sin(0:0.1:2*pi), 'r');
fill(l3 + 0.01*cos(0:0.1:2*pi), 0.01*sin(0:0.1:2*pi), 'r');
fill(-1/2 + 0.01*cos(0:0.1:2*pi), sqrt(3)/2 + 0.01*sin(0:0.1:2*pi), 'r');
fill(-1/2 + 0.01*cos(0:0.1:2*pi), -sqrt(3)/2 + 0.01*sin(0:0.1:2*pi), 'r');
fill(0.03*cos(0:0.1:2*pi), 0.03*sin(0:0.1:2*pi), 'b');
fill(-1 + 0.02*cos(0:0.1:2*pi), 0.02*sin(0:0.1:2*pi), 'w');

grid on
xlabel('x');
ylabel('y');
legend({'C = 3.01', 'C = 3.0242', 'C = 3.08', 'C = 3.19', 'C = 3.33'}, 'Location', 'Southeast');
title('Curves of zero velocity for different values of Jacobi constant');

hold off

%% Functions for the script
function res = JacobDiff(x, mu, Cj)
    res = x(1)^2 + x(2)^2 + 2*(1 - mu)/sqrt((x(1) - mu)^2 + x(2)^2) + 2*mu/sqrt((x(1) - mu + 1)^2 + x(2)^2) + mu*(1 - mu) - Cj;
end

function sol = DJacobDiff(x, mu)
    pd1 = 2*x(1) - 2*(1 - mu)*(x(1) - mu)/((x(1) - mu)^2 + x(2)^2)^(3/2) - 2*mu*(x(1) - mu + 1)/((x(1) - mu + 1)^2 + x(2)^2)^(3/2);
    pd2 = 2*x(2) - 2*(1 - mu)*x(2)/((x(1) - mu)^2 + x(2)^2)^(3/2) - 2*mu*x(2)/((x(1) - mu + 1)^2 + x(2)^2)^(3/2);
    sol = [pd1, pd2];
end

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