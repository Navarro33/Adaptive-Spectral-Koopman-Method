%%% Test on Check Points %%%
clear; close all; clc;

addpath('funcs');
addpath('models');

%% Slow Manifold
n = 200:200:2000;
T = 20;
N = 3;
x0 = [1,1];
r = [0.1, 0.1];
frac = 0.2;
op = 1;

p = struct();
p.alpha = -0.05;
p.beta = -1;

fmodel = @(t, x, p) SlowManifold(t, x, p);
fdynamics = @(X, Y, p) SlowManifoldDynamics(X, Y, p);

% test
running_time = zeros(1, length(n));
error_terminal = zeros(2, length(n));

xtT = [exp(p.alpha*T); (-2*p.alpha / (p.beta - 2*p.alpha)) * exp(p.beta*T) + (p.beta / (p.beta - 2*p.alpha)) * exp(2*p.alpha*T)];
for i = 1:length(n)
    tic;
    [xnT, n_decomp] = ASK_2D(fdynamics, p, n(i), T, N, x0, r, frac, op);
    running_time(i) = toc;
    
    error_terminal(:, i) = abs(xtT - xnT);
end

% display & store
model_name = 'SM';
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name1 = strcat('results\CheckPoints\time_', model_name, '_n_', num2str(T));
name2 = strcat('results\CheckPoints\error_', model_name, '_n_', num2str(T));

save(strcat(name1, format1), 'running_time');
save(strcat(name2, format1), 'error_terminal');

figure;
plot(n, running_time, 'o-', 'linewidth', 2, 'markersize', 10);
xlim([100, 2100]);
xlabel('n', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Running Time (s)', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name1, format2));
saveas(gcf, strcat(name1, format3));

NameArray = {'Marker'};
ValueArray = {'o'; 's'};
figure;
p = semilogy(n, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([100, 2100]);
xlabel('n', 'interpreter', 'latex', 'fontsize', 18);
ylim([5e-16, 1e-14]);
ylabel('Error', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name2, format2));
saveas(gcf, strcat(name2, format3));

disp('Slow manifold completed!');

%% Lotka-Volterra
n = 200:200:2000;
T = 20;
N = 5;
x0 = [10,5];
r = [0.5, 0.5];
frac = 0.5;
op = 1;

p = struct();
p.alpha = 1.1;
p.beta = 0.4;
p.gamma = 0.4;
p.delta = 0.1;

fmodel = @(t, x, p) LotkaVolterraModel(t, x, p);
fdynamics = @(X, Y, p) LotkaVolterraDynamics(X, Y, p);

% test
running_time = zeros(1, length(n));
error_terminal = zeros(2, length(n));

t2 = linspace(0, T, 20000+1);
xt = RungeKutta9(fmodel, x0, t2, p);
for i = 1:length(n)
    tic;
    [xnT, n_decomp] = ASK_2D(fdynamics, p, n(i), T, N, x0, r, frac, op);
    running_time(i) = toc;
    
    error_terminal(:, i) = abs(xt(:, end) - xnT);
end

% display & store
model_name = 'LV';
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name1 = strcat('results\CheckPoints\time_', model_name, '_n_', num2str(T));
name2 = strcat('results\CheckPoints\error_', model_name, '_n_', num2str(T));

save(strcat(name1, format1), 'running_time');
save(strcat(name2, format1), 'error_terminal');

figure;
plot(n, running_time, 'o-', 'linewidth', 2, 'markersize', 10);
xlim([100, 2100]);
xlabel('n', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Running Time (s)', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name1, format2));
saveas(gcf, strcat(name1, format3));

NameArray = {'Marker'};
ValueArray = {'o'; 's'};
figure;
p = semilogy(n, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([100, 2100]);
xlabel('n', 'interpreter', 'latex', 'fontsize', 18);
ylim([1e-8, 1e-6]);
ylabel('Error', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name2, format2));
saveas(gcf, strcat(name2, format3));

disp('Lotka-Volterrra model completed!');

%% Simple pendulum
clear;
n = 200:200:2000;
T = 20;
N = 7;
x0 = [- pi/4, pi/6];
r = [pi/8, pi/12];
%r = [pi/8, pi/8];
frac = 0.2;
op = 1;

p = struct();
p.g = 9.8;
p.L = 9.8;

fmodel = @(t, x, p) SimplePendulum(t, x, p);
fdynamics = @(X, Y, p) SimplePendulumDynamics(X, Y, p);

% test
running_time = zeros(1, length(n));
error_terminal = zeros(2, length(n));

t2 = linspace(0, T, 20000+1);
xt = RungeKutta9(fmodel, x0, t2, p);
for i = 1:length(n)
    tic;
    [xnT, n_decomp] = ASK_2D(fdynamics, p, n(i), T, N, x0, r, frac, op);
    running_time(i) = toc;
    
    error_terminal(:, i) = abs(xt(:, end) - xnT);
end

% display & store
model_name = 'SP';
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name1 = strcat('results\CheckPoints\time_', model_name, '_n_', num2str(T));
name2 = strcat('results\CheckPoints\error_', model_name, '_n_', num2str(T));

save(strcat(name1, format1), 'running_time');
save(strcat(name2, format1), 'error_terminal');

figure;
plot(n, running_time, 'o-', 'linewidth', 2, 'markersize', 10);
xlim([100, 2100]);
xlabel('n', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Running Time (s)', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name1, format2));
saveas(gcf, strcat(name1, format3));

NameArray = {'Marker'};
ValueArray = {'o'; 's'};
figure;
p = semilogy(n, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([100, 2100]);
xlabel('n', 'interpreter', 'latex', 'fontsize', 18);
ylim([0.8e-8, 2.8e-8]);
ylabel('Error', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name2, format2));
saveas(gcf, strcat(name2, format3));

disp('Simple pendulum model completed!');

%% Limit Cycle
clear;
n = 200:200:2000;
T = 20;
N = 7;
x0 = [sqrt(2)/2, -sqrt(2)/2];
r = [sqrt(2)/6, sqrt(2)/6];
frac = 0.2;
op = 1;

p = struct();
p.k = 1;

fmodel = @(t, x, p) LimitCycle(t, x, p);
fdynamics = @(X, Y, p) LimitCycleDynamics(X, Y, p);

% test
running_time = zeros(1, length(n));
error_terminal = zeros(2, length(n));

xtT = solve_limit_cycle(T, x0, p)'; %closed-form
for i = 1:length(n)
    tic;
    [xnT, n_decomp] = ASK_2D(fdynamics, p, n(i), T, N, x0, r, frac, op);
    running_time(i) = toc;
    
    error_terminal(:, i) = abs(xtT - xnT);
end

% display & store
model_name = 'LC';
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name1 = strcat('results\CheckPoints\time_', model_name, '_n_', num2str(T));
name2 = strcat('results\CheckPoints\error_', model_name, '_n_', num2str(T));

save(strcat(name1, format1), 'running_time');
save(strcat(name2, format1), 'error_terminal');

figure;
plot(n, running_time, 'o-', 'linewidth', 2, 'markersize', 10);
xlim([100, 2100]);
xlabel('n', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Running Time (s)', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name1, format2));
saveas(gcf, strcat(name1, format3));

NameArray = {'Marker'};
ValueArray = {'o'; 's'};
figure;
p = semilogy(n, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([100, 2100]);
xlabel('n', 'interpreter', 'latex', 'fontsize', 18);
ylim([1.5e-7, 1.2e-6]);
ylabel('Error', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name2, format2));
saveas(gcf, strcat(name2, format3));

disp('Limit cycle model completed!');
