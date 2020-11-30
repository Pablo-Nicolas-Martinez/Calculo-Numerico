%% Parameters for the problem
mu = 0.0121529;
C = [3.01, 3.0242, 3.08, 3.19, 3.3];

%% Stable points of the Lagrange solution
%L1 = ;
%L2 = ;
%L3 = ;

%% Orbits for each Jacobi constant
X1 = contzvc([-1/2; sqrt(3)/2], 0.001, 12, 1, @(x) JacobDiff(x, mu, C(1)), @(x) DJacobDiff(x, mu));
X2 = contzvc([-1/2; sqrt(3)/2], 0.001, 12, 1, @(x) JacobDiff(x, mu, C(2)), @(x) DJacobDiff(x, mu));
X3 = contzvc([-1/2; sqrt(3)/2], 0.001, 12, 1, @(x) JacobDiff(x, mu, C(3)), @(x) DJacobDiff(x, mu));
X4 = contzvc([-1/2; sqrt(3)/2], 0.001, 12, 1, @(x) JacobDiff(x, mu, C(4)), @(x) DJacobDiff(x, mu));
X5 = contzvc([-1/2; sqrt(3)/2], 0.001, 12, 1, @(x) JacobDiff(x, mu, C(5)), @(x) DJacobDiff(x, mu));

%% Plots for the curves found
hold on
plot(X1(1,:), X1(2,:), 'y')
plot(X1(1,:), -X1(2,:), 'y')
plot(X2(1,:), X2(2,:), 'b')
plot(X3(1,:), X3(2,:), 'g')
plot(X4(1,:), X4(2,:), 'r')
plot(X5(1,:), X5(2,:), 'k')
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