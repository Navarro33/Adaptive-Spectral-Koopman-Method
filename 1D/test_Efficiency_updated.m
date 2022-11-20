%%% Efficiency Test %%%
clear; close all; clc;

addpath('funcs');
addpath('models');

%% Parameters
T = 10;
x0 = pi/4;
p = struct();
p.k = -0.5;

fmodel = @(t, x, p) CosineSquareModel(t, x, p);
fdynamics = @(x, p) CosineSquareDynamics(x, p);

%%
Ref = atan(p.k * T + tan(x0));

%% ASK
n = 50;
r = pi/20;
frac = 0.2;
op = 1;

N = 5:2:11;
err_ASK = zeros(1, length(N));
n_call_ASK = zeros(1, length(N));
for s = 1:length(N)
    [sol, n_decomp] = ASK_1D(fdynamics, p, n, T, N(s), x0, r, frac, op);
    err_ASK(s) = abs(sol - Ref);
    
    n_call_ASK(s) = N(s) * n_decomp;
end

%% Runge-Kutta 4th Order
m = 5;
n_RK4 = zeros(1, m);
err_RK4 = zeros(1, m);
n_call_RK4 = zeros(1, m);

n_RK4(1) = 20;
for s = 1:m
    n_RK4(s) = 2^(s - 1)*n_RK4(1);
    
    t = linspace(0, T, n_RK4(s)+1);
    sol = RungeKutta4(fmodel, x0, t, p);
    err_RK4(s) = abs(sol(end) - Ref);
    
    n_call_RK4(s) = 4*n_RK4(s);
end


%% Runge-Kutta 9th Order
m = 4;
n_RK9 = zeros(1, m);
err_RK9 = zeros(1, m);
n_call_RK9 = zeros(1, m);

n_RK9(1) = 5;
for s = 1:m
    n_RK9(s) = 2^(s - 1)*n_RK9(1);
    
    t = linspace(0, T, n_RK9(s)+1);
    sol = RungeKutta9(fmodel, x0, t, p);
    err_RK9(s) = abs(sol(end) - Ref);
    
    n_call_RK9(s) = 14*n_RK9(s);
end

%% Euler
m = 4;
n_Euler = zeros(1, m);
err_Euler = zeros(1, m);
n_call_Euler = zeros(1, m);

n_Euler(1) = 1000;
for s = 1:m
    n_Euler(s) = 2^(s - 1)*n_Euler(1);
    
    t = linspace(0, T, n_Euler(s)+1);
    sol = Euler(fmodel, x0, t, p);
    err_Euler(s) = abs(sol(end) - Ref);
    
    n_call_Euler(s) = n_Euler(s);
end

%% Adams-Bashforth 5 step
m = 5;
n_AB5 = zeros(1, m);
err_AB5 = zeros(1, m);
n_call_AB5 = zeros(1, m);

n_AB5(1) = 50;
for s = 1:m
    n_AB5(s) = 2^(s - 1)*n_AB5(1);
    
    t = linspace(0, T, n_AB5(s)+1);
    sol = AdamsBashforth5(fmodel, x0, t, p);
    err_AB5(s) = abs(sol(end) - Ref);
    
    n_call_AB5(s) = 10 + 5*(n_AB5(s)-5);
end



%% Display
model_name = 'Cosine';
format1 = '.png';
format2 = '.eps';
name = 'results\Efficiency\';

figure;
loglog(n_call_ASK, err_ASK, 'o-', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
loglog(n_call_RK4, err_RK4, '^-', 'MarkerSize', 10, 'LineWidth', 2);
loglog(n_call_RK9, err_RK9, 'd-', 'MarkerSize', 10, 'LineWidth', 2);
loglog(n_call_Euler, err_Euler, '+-', 'MarkerSize', 10, 'LineWidth', 2);
loglog(n_call_AB5, err_AB5, 's-', 'MarkerSize', 10, 'LineWidth', 2);
set(gca, 'fontsize', 18);
set(gca, 'linewidth', 2);
%legend('ASK', 'RK4', 'RK9', 'Euler', 'AB5', 'Location', 'best');
xlabel('Number of function calls', 'fontsize', 18);
ylabel('Error', 'fontsize', 18);
xlim([40, 1e4]);
ylim([1e-15, 5e-2]);

saveas(gcf, strcat(name, 'efficiency_', model_name, format1));
saveas(gcf, strcat(name, 'efficiency_', model_name, format2));

