%%% Efficiency Test %%%
clear; close all; clc;

addpath('funcs');
addpath('models');

%% Parameters
T = 10;
x0 = pi/4;
p = struct();
p.k = -0.5;

fmodel = @(t, x, p) CosineSquareModel_slow(t, x, p);
fdynamics = @(x, p) CosineSquareDynamics_slow(x, p);

%%
Ref = atan(p.k * T + tan(x0));

%% ASK
n = 50;
r = pi/20;
frac = 0.2;
op = 1;

N = 5:2:11;
err_ASK = zeros(1, length(N));
time_ASK = zeros(1, length(N));
n_loop = 50;
for s = 1:length(N)
    [sol, n_decomp] = ASK_1D(fdynamics, p, n, T, N(s), x0, r, frac, op);
    
    tic;
    for k = 1:n_loop
        [sol, n_decomp] = ASK_1D(fdynamics, p, n, T, N(s), x0, r, frac, op);
    end
    elapse = toc / n_loop;
    
    err_ASK(s) = abs(sol - Ref);
    time_ASK(s) = elapse;
end

%% Runge-Kutta 4th Order
m = 5;
n_RK4 = zeros(1, m);
err_RK4 = zeros(1, m);
time_RK4 = zeros(1, m);

n_RK4(1) = 50;
for s = 1:m
    n_RK4(s) = 2^(s - 1)*n_RK4(1);
    
    tic;
    t = linspace(0, T, n_RK4(s)+1);
    sol = RungeKutta4(fmodel, x0, t, p);
    elapse = toc;
    
    err_RK4(s) = abs(sol(end) - Ref);
    time_RK4(s) = elapse;
end


%% Runge-Kutta 9th Order
m = 4;
n_RK9 = zeros(1, m);
err_RK9 = zeros(1, m);
time_RK9 = zeros(1, m);

n_RK9(1) = 5;
for s = 1:m
    n_RK9(s) = 2^(s - 1)*n_RK9(1);
    
    tic;
    t = linspace(0, T, n_RK9(s)+1);
    sol = RungeKutta9(fmodel, x0, t, p);
    elapse = toc;
    
    err_RK9(s) = abs(sol(end) - Ref);
    time_RK9(s) = elapse;
end

%% Euler
m = 5;
n_Euler = zeros(1, m);
err_Euler = zeros(1, m);
time_Euler = zeros(1, m);

n_Euler(1) = 2000;
for s = 1:m
    n_Euler(s) = 2^(s - 1)*n_Euler(1);
    
    tic;
    t = linspace(0, T, n_Euler(s)+1);
    sol = Euler(fmodel, x0, t, p);
    elapse = toc;
    
    err_Euler(s) = abs(sol(end) - Ref);
    time_Euler(s) = elapse;
end

%% Adams-Bashforth 5 step
m = 5;
n_AB5 = zeros(1, m);
err_AB5 = zeros(1, m);
time_AB5 = zeros(1, m);

n_AB5(1) = 150;
for s = 1:m
    n_AB5(s) = 2^(s - 1)*n_AB5(1);
    
    tic;
    t = linspace(0, T, n_AB5(s)+1);
    sol = AdamsBashforth5(fmodel, x0, t, p);
    elapse = toc;
    
    err_AB5(s) = abs(sol(end) - Ref);
    time_AB5(s) = elapse;
end



%% Display
model_name = 'Cosine';
format1 = '.png';
format2 = '.eps';
name = 'results\Efficiency\';

figure;
loglog(time_ASK, err_ASK, 'o-', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
loglog(time_RK4, err_RK4, '^-', 'MarkerSize', 10, 'LineWidth', 2);
loglog(time_RK9, err_RK9, 'd-', 'MarkerSize', 10, 'LineWidth', 2);
loglog(time_Euler, err_Euler, '+-', 'MarkerSize', 10, 'LineWidth', 2);
loglog(time_AB5, err_AB5, 's-', 'MarkerSize', 10, 'LineWidth', 2);
set(gca, 'fontsize', 18);
set(gca, 'linewidth', 2);
%legend('ASK', 'RK4', 'RK9', 'Euler', 'AB5', 'Location', 'best');
xlabel('Time (s)', 'fontsize', 18);
ylabel('Error', 'fontsize', 18);
xlim([1e-2, 1e2]);
ylim([1e-15, 1e-2]);

saveas(gcf, strcat(name, 'efficiency_time_', model_name, format1));
saveas(gcf, strcat(name, 'efficiency_time_', model_name, format2));

