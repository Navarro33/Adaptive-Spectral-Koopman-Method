%%% Running Time - Kraichnan-Orszag %%%
clear; close all; clc;

addpath('funcs');
addpath('models');


%% Parameters
T = 10;
x0 = [1, 2, -3];

p = struct();
p.k1 = 1;
p.k2 = 1;
p.k3 = 1;

fmodel = @(t, x, p) KraichnanOrszagModel_slow(t, x, p);
fdynamics = @(X, Y, Z, p) KraichnanOrszagDynamics_slow(X, Y, Z, p);

% reference
t = linspace(0, T, 10000+1);
ref = RungeKutta9(@(t, x, p) KraichnanOrszagModel(t, x, p), x0, t, p);

%% ASK
n = 200;
r = [0.3, 0.3, 0.3];
frac = 0.15;
op = 1;

N = 3:2:7;
err_ASK = zeros(3, length(N));
time_ASK = zeros(1, length(N));
n_loop = 10;
for s = 1:length(N)
    [sol, n_decomp] = ASK_3D(fdynamics, p, n, T, N(s), x0, r, frac, op);
    
    tic;
    for k = 1:n_loop
        [sol, n_decomp] = ASK_3D(fdynamics, p, n, T, N(s), x0, r, frac, op);
    end
    elapse = toc / n_loop;
    
    err_ASK(:, s) = abs(sol - ref(:,end));
    time_ASK(s) = elapse;
end

%% Runge-Kutta 4th Order
m = 5;
n_RK4 = zeros(1, m);
err_RK4 = zeros(3, m);
time_RK4 = zeros(1, m);

n_RK4(1) = 123;
for s = 1:m
    n_RK4(s) = 2^(s - 1)*n_RK4(1);
    
    tic;
    t = linspace(0, T, n_RK4(s)+1);
    sol = RungeKutta4(fmodel, x0, t, p);
    elapse = toc;
    
    err_RK4(:, s) = abs(sol(:, end) - ref(:, end));
    time_RK4(s) = elapse;
end

%% Runge-Kutta 9th Order
m = 4;
n_RK9 = zeros(1, m);
err_RK9 = zeros(3, m);
time_RK9 = zeros(1, m);

n_RK9(1) = 35;
for s = 1:m
    n_RK9(s) = 2^(s - 1)*n_RK9(1);
    
    tic;
    t = linspace(0, T, n_RK9(s)+1);
    sol = RungeKutta9(fmodel, x0, t, p);
    elapse = toc;
    
    err_RK9(:, s) = abs(sol(:, end) - ref(:, end));
    time_RK9(s) = elapse;
end

%% Euler
m = 5;
n_Euler = zeros(1, m);
err_Euler = zeros(3, m);
time_Euler = zeros(1, m);

n_Euler(1) = 1000;
for s = 1:m
    n_Euler(s) = 2^(s - 1)*n_Euler(1);
    
    tic;
    t = linspace(0, T, n_Euler(s)+1);
    sol = Euler(fmodel, x0, t, p);
    elapse = toc;
    
    err_Euler(:, s) = abs(sol(:, end) - ref(:, end));
    time_Euler(s) = elapse;
end


%% Adams-Bashforth 5 step
m = 5;
n_AB5 = zeros(1, m);
err_AB5 = zeros(3, m);
time_AB5 = zeros(1, m);

n_AB5(1) = 150;
for s = 1:m
    n_AB5(s) = 2^(s - 1)*n_AB5(1);
    
    tic;
    t = linspace(0, T, n_AB5(s)+1);
    sol =  AdamsBashforth5(fmodel, x0, t, p);
    elapse = toc;
    
    err_AB5(:, s) = abs(sol(:, end) - ref(:, end));
    time_AB5(s) = elapse;
end


%% Display
model_name = 'KO';
format1 = '.png';
format2 = '.eps';
name = 'results\Efficiency\';

figure;
loglog(time_ASK, sqrt((err_ASK(1,:).^2 + err_ASK(2,:).^2 + err_ASK(3,:).^2) / 3), 'o-', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
loglog(time_RK4, sqrt((err_RK4(1,:).^2 + err_RK4(2,:).^2 + + err_RK4(3,:).^2) / 3), '^-', 'MarkerSize', 10, 'LineWidth', 2);
loglog(time_RK9, sqrt((err_RK9(1,:).^2 + err_RK9(2,:).^2 + err_RK9(3,:).^2) / 3), 'd-', 'MarkerSize', 10, 'LineWidth', 2);
loglog(time_Euler, sqrt((err_Euler(1,:).^2 + err_Euler(2,:).^2 + err_Euler(3,:).^2) / 3), '+-', 'MarkerSize', 10, 'LineWidth', 2);
loglog(time_AB5, sqrt((err_AB5(1,:).^2 + err_AB5(2,:).^2 + + err_AB5(3,:).^2) / 3), 's-', 'MarkerSize', 10, 'LineWidth', 2);
set(gca, 'fontsize', 18);
set(gca, 'linewidth', 2);
%legend('ASK', 'RK4', 'RK9', 'Euler', 'AB5', 'Location', 'best');
xlabel('Time (s)', 'fontsize', 18);
ylabel('Error', 'fontsize', 18);
xlim([2e-1, 2e1]);
xticks(10.^(-1:1:1));
ylim([1e-11, 1e1]);
yticks(10.^(-10:4:0));

saveas(gcf, strcat(name, 'efficiency_time_', model_name, format1));
saveas(gcf, strcat(name, 'efficiency_time_', model_name, format2));



