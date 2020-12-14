%% Function to compute the accessible areas
clear all
close all

%% Parameters for the problem
m = 0.0121529; % Reduced mass
C = 3.08; % Jacobi constant

%% Computation of the values
X = contzvc([0.75; 0.75], 0.001, 12, 1, @(x) JacobDiff(x, m, C), @(x) DJacobDiff(x, m));

%% Plot and visualization
fill(X(1,:), X(2,:), 'r');
grid on;
xlabel('x');
ylabel('y');
title('Hill´s region´s complement for C = 3.08');

%% Functions for the problem
function res = JacobDiff(x, mu, Cj)
    res = x(1)^2 + x(2)^2 + 2*(1 - mu)/sqrt((x(1) - mu)^2 + x(2)^2) + 2*mu/sqrt((x(1) - mu + 1)^2 + x(2)^2) + mu*(1 - mu) - Cj;
end

function sol = DJacobDiff(x, mu)
    pd1 = 2*x(1) - 2*(1 - mu)*(x(1) - mu)/((x(1) - mu)^2 + x(2)^2)^(3/2) - 2*mu*(x(1) - mu + 1)/((x(1) - mu + 1)^2 + x(2)^2)^(3/2);
    pd2 = 2*x(2) - 2*(1 - mu)*x(2)/((x(1) - mu)^2 + x(2)^2)^(3/2) - 2*mu*x(2)/((x(1) - mu + 1)^2 + x(2)^2)^(3/2);
    sol = [pd1, pd2];
end