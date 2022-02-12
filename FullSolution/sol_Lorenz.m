%%% Full Solutions Lorenz %%%
clear; close all; clc;

addpath('funcs');
addpath('models');

%% Parameters
m = 2000;
T = 20;
tspan = linspace(0, T, m+1);
n = 2000;
T = 20;
N = 5;
x0 = [5,5,5];
r = [1,1,1];
frac = 0.75;
op = 1;

p.alpha = 10;
p.beta = 28;
p.gamma = 3;

fmodel = @(t, x, p) LorenzModel(t, x, p);
fdynamics = @(X, Y, Z, p) LorenzDynamics(X, Y, Z, p);

%% Koopman
tic;
[xn, n_decomp] = ASK_Full_3D(fdynamics, p, n, tspan, N, x0, r, frac, op);
elapsed = toc;
clc;
fprintf('--- Koopman completed in %.4f seconds --- \n', elapsed);

%% Runge-Kutta & Reference
tic;
xt = RungeKutta4(fmodel, x0, tspan, p);
elapsed2 = toc;
fprintf('--- Runge-Kutta completed in %.4f seconds --- \n', elapsed2);

M = 20000;
tspan2 = linspace(0, T, M+1);
xr = RungeKutta9(fmodel, x0, tspan2, p);

is_member = ismember(tspan2, tspan);
indexes = find(is_member);
xr = xr(:, indexes);

%% Errors
Koopman_errors = abs(xn - xr);
RungeKutta_errors = abs(xt - xr);

%% Display & Save
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name = 'results\Lorenz\';

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
plot(tspan(1:8:end), xn(1, 1:8:end), 'o', 'MarkerSize', 3.5);
hold on
plot(tspan(1:8:end), xt(1, 1:8:end), '^', 'MarkerSize', 3.5);
plot(tspan(1:8:end), xr(1, 1:8:end), '--', 'LineWidth', 2);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x1', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
% set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, strcat(name, 'performance_', 'x1', format2));
saveas(gcf, strcat(name, 'performance_', 'x1', format3));

% performance x2
figure;
plot(tspan(1:8:end), xn(2, 1:8:end), 'o', 'MarkerSize', 3.5);
hold on
plot(tspan(1:8:end), xt(2, 1:8:end), '^', 'MarkerSize', 3.5);
plot(tspan(1:8:end), xr(2, 1:8:end), '--', 'LineWidth', 2);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x2', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(name, 'performance_', 'x2', format2));
saveas(gcf, strcat(name, 'performance_', 'x2', format3));

% performance x3
figure;
plot(tspan(1:8:end), xn(3, 1:8:end), 'o', 'MarkerSize', 3.5);
hold on
plot(tspan(1:8:end), xt(3, 1:8:end), '^', 'MarkerSize', 3.5);
plot(tspan(1:8:end), xr(3, 1:8:end), '--', 'LineWidth', 2);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x3', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(name, 'performance_', 'x3', format2));
saveas(gcf, strcat(name, 'performance_', 'x3', format3));

% error x1
figure;
semilogy(tspan(1:4:end), Koopman_errors(1, 1:4:end), 'o', 'MarkerSize', 5);
hold on
semilogy(tspan(1:4:end), RungeKutta_errors(1, 1:4:end), '^', 'MarkerSize', 5);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x1 Error', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(name, 'error_', 'x1', format2));
saveas(gcf, strcat(name, 'error_', 'x1', format3));

% error x2
figure;
semilogy(tspan(1:4:end), Koopman_errors(2, 1:4:end), 'o', 'MarkerSize', 5);
hold on
semilogy(tspan(1:4:end), RungeKutta_errors(2, 1:4:end), '^', 'MarkerSize', 5);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x2 Error', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(name, 'error_', 'x2', format2));
saveas(gcf, strcat(name, 'error_', 'x2', format3));

% error x3
figure;
semilogy(tspan(1:4:end), Koopman_errors(3, 1:4:end), 'o', 'MarkerSize', 5);
hold on
semilogy(tspan(1:4:end), RungeKutta_errors(3, 1:4:end), '^', 'MarkerSize', 5);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x3 Error', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(name, 'error_', 'x3', format2));
saveas(gcf, strcat(name, 'error_', 'x3', format3));

% 3D plot
figure;
pp = plot3(xr(1,:), xr(2, :), xr(3, :), 'Color', 'r', 'LineWidth', 2);
direction = [0.25, 0.25, 1];
rotate(pp, direction, 30);
xlabel('x1', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x2', 'Interpreter', 'latex', 'FontSize', 18);
zlabel('x3', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(name, '3D', format2));
saveas(gcf, strcat(name, '3D', format3));



