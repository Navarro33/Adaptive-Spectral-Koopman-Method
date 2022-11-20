%%% Full Solutions Simple Pendulum %%%
clear; close all; clc;

addpath('funcs');
addpath('models');

%% Parameters
m = 1000;                           % number of time steps
T = 20;                             % terminal time
tspan = linspace(0, T, m+1);
n = 200;                            % number of check points
N = 7;                              % number of Gauss-Lobatto points
x0 = [- pi/4, pi/6];
r = [pi/8, pi/12];
frac = 0.2;
op = 1;

p = struct();
p.g = 9.8;
p.L = 9.8;

fmodel = @(t, x, p) SimplePendulum(t, x, p);
fdynamics = @(X, Y, p) SimplePendulumDynamics(X, Y, p);

%% Reference
M = 3000;
tspan2 = linspace(0, T, M+1);
x_Ref = RungeKutta9(fmodel, x0, tspan2, p);
is_member = ismember(tspan2, tspan);
indexes = find(is_member);
x_Ref = x_Ref(:, indexes);

%% ASK
tic;
[xnT, ~] = ASK_2D(fdynamics, p, n, T, N, x0, r, frac, op);
t_ASK = toc;
[x_ASK, n_decomp] = ASK_Full_2D(fdynamics, p, n, tspan, N, x0, r, frac, op);

%% Other algorithms

% Euler 
tic;
x_Euler = Euler(fmodel, x0, tspan, p);
t_Euler = toc;

% 4th order Runge-Kutta
tic;
x_RK4 = RungeKutta4(fmodel, x0, tspan, p);
t_RK4 = toc;

% 5-step Adams-Bashforth
tic;
x_AB5 = AdamsBashforth5(fmodel, x0, tspan, p);
t_AB5 = toc;

% 4-step Adams-Moulton 
tic;
x_AM4 = AdamsMoulton4(fmodel, x0, tspan, p);
t_AM4 = toc;

%% Errors
e_ASK = abs(x_ASK - x_Ref);
e_Euler = abs(x_Euler - x_Ref);
e_RK4 = abs(x_RK4 - x_Ref);
e_AB5 = abs(x_AB5 - x_Ref);
e_AM4 = abs(x_AM4 - x_Ref);

%% Results
format1 = '.png';
format2 = '.eps';
file = 'results\SimplePendulum\';

% time
fprintf('***** Running Time *****\n');
fprintf('ASK: %.6f \n', t_ASK);
fprintf('Euler: %.6f \n', t_Euler);
fprintf('RK4: %.6f \n', t_RK4);
fprintf('AB5: %.6f \n', t_AB5);
fprintf('AM4: %.6f \n', t_AM4);
fprintf('******************************\n\n');

% integration error
format short e
Algorithms = ["ASK"; "Euler"; "RK4"; "AB5"; "AM4"];
e1 = [e_ASK(1, end); e_Euler(1, end); e_RK4(1, end); e_AB5(1, end); e_AM4(1, end)];
e2 = [e_ASK(2, end); e_Euler(2, end); e_RK4(2, end); e_AB5(2, end); e_AM4(2, end)];
fprintf('********** Running Time **********\n\n');
tab = table(Algorithms, e1, e2);
disp(tab);
fprintf('******************************\n\n');
writetable(tab, strcat(file, 'SimplePendulum.csv'));

% performance x1
figure;
plot(tspan, x_Ref(1, :), '--', 'LineWidth', 2);
hold on
plot(tspan, x_ASK(1, :), 'o', 'MarkerSize', 7);
plot(tspan, x_Euler(1, :), '+', 'MarkerSize', 7);
plot(tspan, x_RK4(1, :), '^', 'MarkerSize', 7);
plot(tspan, x_AB5(1, :), 's', 'MarkerSize', 7);
plot(tspan, x_AM4(1, :), '*', 'MarkerSize', 7);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylim([-1.5, 1.5]);
ylabel('x1', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(file, 'performance_', 'x1', format1));
saveas(gcf, strcat(file, 'performance_', 'x1', format2));

% performance x2
figure;
plot(tspan, x_Ref(2, :), '--', 'LineWidth', 2);
hold on
plot(tspan, x_ASK(2, :), 'o', 'MarkerSize', 7);
plot(tspan, x_Euler(2, :), '+', 'MarkerSize', 7);
plot(tspan, x_RK4(2, :), '^', 'MarkerSize', 7);
plot(tspan, x_AB5(2, :), 's', 'MarkerSize', 7);
plot(tspan, x_AM4(2, :), '*', 'MarkerSize', 7);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylim([-1.5, 1.5]);
ylabel('x2', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(file, 'performance_', 'x2', format1));
saveas(gcf, strcat(file, 'performance_', 'x2', format2));

% error x1
figure;
semilogy(tspan, e_ASK(1, :), 'o', 'MarkerSize', 5);
hold on
semilogy(tspan, e_Euler(1, :), '+', 'MarkerSize', 5);
semilogy(tspan, e_RK4(1, :), '^', 'MarkerSize', 5);
semilogy(tspan, e_AB5(1, :), 's', 'MarkerSize', 5);
semilogy(tspan, e_AM4(1, :), '*', 'MarkerSize', 5);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylim([1e-12, 1e-2]);
yticks(10.^(-12:2:-2));
ylabel('x1 Error', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(file, 'error_', 'x1', format1));
saveas(gcf, strcat(file, 'error_', 'x1', format2));

% error x2
figure;
semilogy(tspan, e_ASK(2, :), 'o', 'MarkerSize', 5);
hold on
semilogy(tspan, e_Euler(2, :), '+', 'MarkerSize', 5);
semilogy(tspan, e_RK4(2, :), '^', 'MarkerSize', 5);
semilogy(tspan, e_AB5(2, :), 's', 'MarkerSize', 5);
semilogy(tspan, e_AM4(2, :), '*', 'MarkerSize', 5);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylim([1e-12, 1e-2]);
yticks(10.^(-12:2:-2));
ylabel('x2 Error', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(file, 'error_', 'x2', format1));
saveas(gcf, strcat(file, 'error_', 'x2', format2));

% 3D plot
figure;
plot3(x_ASK(1,:), x_ASK(2, :), tspan, 'o', 'LineWidth', 2);
hold on
plot3(x_Ref(1,:), x_Ref(2, :), tspan, '-', 'LineWidth', 2);
hold off
xlabel('x1', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x2', 'Interpreter', 'latex', 'FontSize', 18);
zlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(file, '3D', format1));
saveas(gcf, strcat(file, '3D', format2));
