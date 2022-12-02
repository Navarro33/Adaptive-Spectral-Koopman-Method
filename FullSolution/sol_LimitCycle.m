%%% Full Solutions Limit Cycle %%%
clear; close all; clc;

addpath('funcs');
addpath('models');

%% Parameters
m = 200;
T = 20;
tspan = linspace(0, T, m+1);
n = 200;
N = 9;
x0 = [sqrt(2)/2, -sqrt(2)/2];
r = [sqrt(2)/8, sqrt(2)/8];
frac = 0.2;
op = 1;

p.k = 1;

fmodel = @(t, x, p) LimitCycle(t, x, p);
fdynamics = @(X, Y, p) LimitCycleDynamics(X, Y, p);
%% Koopman
tic;
[xn, n_decomp] = ASK_Full_2D(fdynamics, p, n, tspan, N, x0, r, frac, op);
elapsed = toc;
%clc;
fprintf('--- Koopman completed in %.4f seconds --- \n', elapsed);

%% Runge-Kutta & Reference
tic;
xt = RungeKutta4(fmodel, x0, tspan, p);
elapsed2 = toc;
fprintf('--- Runge-Kutta completed in %.4f seconds --- \n', elapsed2);

xr = solve_limit_cycle(tspan, x0, p);
%% Errors
Koopman_errors = abs(xn - xr);
RungeKutta_errors = abs(xt - xr);

%% Display & Save
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name = 'results\LimitCycle\';

results = struct();
resuts.KoopmanTime = elapsed;
resuts.RungeKuttaTime = elapsed2;
resuts.KoomanSol = xn;
resuts.RungeKuttaSol = xt;
resuts.ReferenceSol = xr;
results.KoopmanErr = Koopman_errors;
results.RungeKuttaErr = RungeKutta_errors;
save(strcat(name, 'results', format1), 'results');

% performance x1
figure;
plot(tspan, xn(1, :), 'o', 'MarkerSize', 7);
hold on
plot(tspan, xt(1, :), '^', 'MarkerSize', 7);
plot(tspan, xr(1, :), '--', 'LineWidth', 2);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylim([-1.5, 1.5]);
ylabel('x1', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
% set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, strcat(name, 'performance_', 'x1', format2));
saveas(gcf, strcat(name, 'performance_', 'x1', format3));

% performance x2
figure;
plot(tspan, xn(2, :), 'o', 'MarkerSize', 7);
hold on
plot(tspan, xt(2, :), '^', 'MarkerSize', 7);
plot(tspan, xr(2, :), '--', 'LineWidth', 2);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylim([-1.5, 1.5]);
ylabel('x2', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(name, 'performance_', 'x2', format2));
saveas(gcf, strcat(name, 'performance_', 'x2', format3));

% error x1
figure;
semilogy(tspan, Koopman_errors(1, :), 'o', 'MarkerSize', 5);
hold on
semilogy(tspan, RungeKutta_errors(1, :), '^', 'MarkerSize', 5);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylim([1e-12, 1e-2]);
yticks(10.^(-12:2:-2));
ylabel('x1 Error', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(name, 'error_', 'x1', format2));
saveas(gcf, strcat(name, 'error_', 'x1', format3));

% error x2
figure;
semilogy(tspan, Koopman_errors(2, :), 'o', 'MarkerSize', 5);
hold on
semilogy(tspan, RungeKutta_errors(2, :), '^', 'MarkerSize', 5);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylim([1e-12, 1e-2]);
yticks(10.^(-12:2:-2));
ylabel('x2 Error', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(name, 'error_', 'x2', format2));
saveas(gcf, strcat(name, 'error_', 'x2', format3));
%%
% 3D plot
figure;
plot3(xn(1,:), xn(2, :), tspan, 'o', 'LineWidth', 2);
hold on
plot3(xr(1,:), xr(2, :), tspan, '-', 'LineWidth', 2);
hold off
xlabel('x1', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x2', 'Interpreter', 'latex', 'FontSize', 18);
zlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(name, '3D', format2));
saveas(gcf, strcat(name, '3D', format3));





