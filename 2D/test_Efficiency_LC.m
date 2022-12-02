%%% Efficiency Test - Limit Cycle %%%
clear; close all; clc;

addpath('funcs');
addpath('models');

%% Parameters
T = 10;
x0 = [sqrt(2)/2, - sqrt(2)/2];
p = struct();
p.k = 1;

fmodel = @(t, x, p) LimitCycle(t, x, p);
fdynamics = @(X, Y, p) LimitCycleDynamics(X, Y, p);

% reference
ref = solve_limit_cycle(T, x0, p)';


%% ASK
n = 60;
r = [sqrt(2)/10, sqrt(2)/10];
frac = 0.2;
op = 1;

N = 3:2:11;
err_ASK = zeros(2, length(N));
n_call_ASK = zeros(1, length(N));
for s = 1:length(N)
    [sol, n_decomp] = ASK_2D(fdynamics, p, n, T, N(s), x0, r, frac, op);
    err_ASK(:, s) = abs(sol - ref);
    
    n_call_ASK(s) = N(s)^2 * n_decomp;
end

%% Runge-Kutta 4th Order
m = 5;
n_RK4 = zeros(1, m);
err_RK4 = zeros(2, m);
n_call_RK4 = zeros(1, m);

n_RK4(1) = 123;
for s = 1:m
    n_RK4(s) = 2^(s - 1)*n_RK4(1);
    
    t = linspace(0, T, n_RK4(s)+1);
    sol = RungeKutta4(fmodel, x0, t, p);
    err_RK4(:, s) = abs(sol(:, end) - ref);
    
    n_call_RK4(s) = 4*n_RK4(s);
end

%% Runge-Kutta 9th Order
m = 4;
n_RK9 = zeros(1, m);
err_RK9 = zeros(2, m);
n_call_RK9 = zeros(1, m);

n_RK9(1) = 35;
for s = 1:m
    n_RK9(s) = 2^(s - 1)*n_RK9(1);
    
    t = linspace(0, T, n_RK9(s)+1);
    sol = RungeKutta9(fmodel, x0, t, p);
    err_RK9(:, s) = abs(sol(:, end) - ref);
    
    n_call_RK9(s) = 14*n_RK9(s);
end

%% Display
model_name = 'LC';
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name = 'results\Efficiency\';

results = struct();
results.err_ASK = err_ASK;
results.n_call_ASK = n_call_ASK;
results.err_RK4 = err_RK4;
results.n_call_RK4 = n_call_RK4;
results.err_RK9 = err_RK9;
results.n_call_RK9 = n_call_RK9;
save(strcat(name, 'efficiency_', model_name, format1), 'results');

figure;
loglog(n_call_ASK, sqrt((err_ASK(1,:).^2 + err_ASK(2,:).^2) / 2), 'o-', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
loglog(n_call_RK4, sqrt((err_RK4(1,:).^2 + err_RK4(2,:).^2) / 2), 's-', 'MarkerSize', 10, 'LineWidth', 2);
loglog(n_call_RK9, sqrt((err_RK9(1,:).^2 + err_RK9(2,:).^2) / 2), '^-', 'MarkerSize', 10, 'LineWidth', 2);
set(gca, 'fontsize', 18);
set(gca, 'linewidth', 2);
%legend('ASK', 'RK4', 'RK9');
xlabel('Number of function calls', 'fontsize', 18);
ylabel('Error', 'fontsize', 18);
xlim([3e2, 1e4]);
ylim([1e-17, 1e-2]);
yticks(10.^(-16:4:-2));

saveas(gcf, strcat(name, 'efficiency_', model_name, format2));
saveas(gcf, strcat(name, 'efficiency_', model_name, format3));


