%%% Full Solutions Lorenz %%%
clear; close all; clc;

addpath('funcs');
addpath('models');

%% Parameters
m = 2000;                          % number of time steps
T = 20;                             % terminal time
tspan = linspace(0, T, m+1);        
n = 2000;                           % number of check points
N = 5;                              % number of Gauss-Lobatto points
x0 = [5,5,5];
r = [1,1,1];
frac = 0.75;
op = 1;

p.alpha = 10;
p.beta = 28;
p.gamma = 3;

fmodel = @(t, x, p) LorenzModel(t, x, p);
fdynamics = @(X, Y, Z, p) LorenzDynamics(X, Y, Z, p);

%% Reference - Runge-Kutta 9
M = 20000;
tspan2 = linspace(0, T, M+1);
x_Ref = RungeKutta9(fmodel, x0, tspan2, p);
is_member = ismember(tspan2, tspan);
indexes = find(is_member);
x_Ref = x_Ref(:, indexes);

%% ASK
tic;
[xnT, ~] = ASK_3D(fdynamics, p, n, T, N, x0, r, frac, op);
t_ASK = toc;
[x_ASK, n_decomp] = ASK_Full_3D(fdynamics, p, n, tspan, N, x0, r, frac, op);
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
file = 'results\Lorenz\';

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
e3 = [e_ASK(3, end); e_Euler(3, end); e_RK4(3, end); e_AB5(3, end); e_AM4(3, end)];
fprintf('********** Running Time **********\n\n');
tab = table(Algorithms, e1, e2, e3);
disp(tab);
fprintf('******************************\n\n');
writetable(tab, strcat(file, 'Lorenz.csv'));

% performance x1
figure;
plot(tspan(1:8:end), x_Ref(1, 1:8:end), '--', 'MarkerSize', 3.5);
hold on
plot(tspan(1:8:end), x_ASK(1, 1:8:end), 'o', 'MarkerSize', 3.5);
plot(tspan(1:8:end), x_Euler(1, 1:8:end), '+', 'MarkerSize', 3.5);
plot(tspan(1:8:end), x_RK4(1, 1:8:end), '^', 'MarkerSize', 3.5);
plot(tspan(1:8:end), x_AB5(1, 1:8:end), 's', 'MarkerSize', 3.5);
plot(tspan(1:8:end), x_AM4(1, 1:8:end), '*', 'MarkerSize', 3.5);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x1', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(file, 'performance_', 'x1', format1));
saveas(gcf, strcat(file, 'performance_', 'x1', format2));

% performance x2
figure;
plot(tspan(1:8:end), x_Ref(2, 1:8:end), '--', 'MarkerSize', 3.5);
hold on
plot(tspan(1:8:end), x_ASK(2, 1:8:end), 'o', 'MarkerSize', 3.5);
plot(tspan(1:8:end), x_Euler(2, 1:8:end), '+', 'MarkerSize', 3.5);
plot(tspan(1:8:end), x_RK4(2, 1:8:end), '^', 'MarkerSize', 3.5);
plot(tspan(1:8:end), x_AB5(2, 1:8:end), 's', 'MarkerSize', 3.5);
plot(tspan(1:8:end), x_AM4(2, 1:8:end), '*', 'MarkerSize', 3.5);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x2', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(file, 'performance_', 'x2', format1));
saveas(gcf, strcat(file, 'performance_', 'x2', format2));

% performance x3
figure;
plot(tspan(1:8:end), x_Ref(3, 1:8:end), '--', 'MarkerSize', 3.5);
hold on
plot(tspan(1:8:end), x_ASK(3, 1:8:end), 'o', 'MarkerSize', 3.5);
plot(tspan(1:8:end), x_Euler(3, 1:8:end), '+', 'MarkerSize', 3.5);
plot(tspan(1:8:end), x_RK4(3, 1:8:end), '^', 'MarkerSize', 3.5);
plot(tspan(1:8:end), x_AB5(3, 1:8:end), 's', 'MarkerSize', 3.5);
plot(tspan(1:8:end), x_AM4(3, 1:8:end), '*', 'MarkerSize', 3.5);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x3', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(file, 'performance_', 'x3', format1));
saveas(gcf, strcat(file, 'performance_', 'x3', format2));


% error x1
figure;
semilogy(tspan(1:4:end), e_ASK(1, 1:4:end), 'o', 'MarkerSize', 5);
hold on
semilogy(tspan(1:4:end), e_Euler(1, 1:4:end), '+', 'MarkerSize', 5);
semilogy(tspan(1:4:end), e_RK4(1, 1:4:end), '^', 'MarkerSize', 5);
semilogy(tspan(1:4:end), e_AB5(1, 1:4:end), 's', 'MarkerSize', 5);
semilogy(tspan(1:4:end), e_AM4(1, 1:4:end), '*', 'MarkerSize', 5);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x1 Error', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(file, 'error_', 'x1', format1));
saveas(gcf, strcat(file, 'error_', 'x1', format2));

% error x2
figure;
semilogy(tspan(1:4:end), e_ASK(2, 1:4:end), 'o', 'MarkerSize', 5);
hold on
semilogy(tspan(1:4:end), e_Euler(2, 1:4:end), '+', 'MarkerSize', 5);
semilogy(tspan(1:4:end), e_RK4(2, 1:4:end), '^', 'MarkerSize', 5);
semilogy(tspan(1:4:end), e_AB5(2, 1:4:end), 's', 'MarkerSize', 5);
semilogy(tspan(1:4:end), e_AM4(2, 1:4:end), '*', 'MarkerSize', 5);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x2 Error', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(file, 'error_', 'x2', format1));
saveas(gcf, strcat(file, 'error_', 'x2', format2));

% error x3
figure;
semilogy(tspan(1:4:end), e_ASK(3, 1:4:end), 'o', 'MarkerSize', 5);
hold on
semilogy(tspan(1:4:end), e_Euler(3, 1:4:end), '+', 'MarkerSize', 5);
semilogy(tspan(1:4:end), e_RK4(3, 1:4:end), '^', 'MarkerSize', 5);
semilogy(tspan(1:4:end), e_AB5(3, 1:4:end), 's', 'MarkerSize', 5);
semilogy(tspan(1:4:end), e_AM4(3, 1:4:end), '*', 'MarkerSize', 5);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x3 Error', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(file, 'error_', 'x3', format1));
saveas(gcf, strcat(file, 'error_', 'x3', format2));

% 3D plot
figure;
pp = plot3(x_Ref(1,:), x_Ref(2, :), x_Ref(3, :), 'Color', 'r', 'LineWidth', 2);
direction = [0.25, 0.25, 1];
rotate(pp, direction, 30);
xlabel('x1', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('x2', 'Interpreter', 'latex', 'FontSize', 18);
zlabel('x3', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
saveas(gcf, strcat(file, '3D', format1));
saveas(gcf, strcat(file, '3D', format2));



