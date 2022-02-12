%%% Test on Gauss-Lobatto Points %%%
clear; close all; clc;

addpath('funcs');
addpath('models');

%% Slow Manifold
n = 200;
T = 20;
N = [3,5,7,9,11];
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
running_time = zeros(1, length(N));
error_terminal = zeros(2, length(N));

xtT = [exp(p.alpha*T); (-2*p.alpha / (p.beta - 2*p.alpha)) * exp(p.beta*T) + (p.beta / (p.beta - 2*p.alpha)) * exp(2*p.alpha*T)];
for i = 1:length(N)
    tic;
    [xnT, n_decomp] = ASK_2D(fdynamics, p, n, T, N(i), x0, r, frac, op);
    running_time(i) = toc;
    
    error_terminal(:, i) = abs(xtT - xnT);
end

% display & store
model_name = 'SM';
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name1 = strcat('results\GaussLobatto\time_', model_name, '_N_', num2str(T));
name2 = strcat('results\GaussLobatto\error_', model_name, '_N_', num2str(T));

save(strcat(name1, format1), 'running_time');
save(strcat(name2, format1), 'error_terminal');

figure;
plot(N, running_time, 'o-', 'linewidth', 2, 'markersize', 10);
xlim([2, 12]);
xlabel('N', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Running Time (s)', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name1, format2));
saveas(gcf, strcat(name1, format3));

NameArray = {'Marker'};
ValueArray = {'o'; 's'};
figure;
p = semilogy(N-1, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([1, 11]);
xlabel('N', 'interpreter', 'latex', 'fontsize', 18);
ylim([1e-16, 1e-9]);
ylabel('Error', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name2, format2));
saveas(gcf, strcat(name2, format3));

disp('Slow manifold completed!');

%% Lotka-Volterra
n = 200;
T = 20;
N = [3,5,7,9,11];
x0 = [10,5];
r = [1.5, 1.5];
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
running_time = zeros(1, length(N));
error_terminal = zeros(2, length(N));

t2 = linspace(0, T, 20000);
xt = RungeKutta9(fmodel, x0, t2, p);
for i = 1:length(N)
    tic;
    [xnT, n_decomp] = ASK_2D(fdynamics, p, n, T, N(i), x0, r, frac, op);
    running_time(i) = toc;
    
    error_terminal(:, i) = abs(xt(:, end) - xnT);
end

% display & store
model_name = 'LV';
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name1 = strcat('results\GaussLobatto\time_', model_name, '_N_', num2str(T));
name2 = strcat('results\GaussLobatto\error_', model_name, '_N_', num2str(T));

save(strcat(name1, format1), 'running_time');
save(strcat(name2, format1), 'error_terminal');

figure;
plot(N, running_time, 'o-', 'linewidth', 2, 'markersize', 10);
xlim([2, 12]);
xlabel('N', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Running Time (s)', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name1, format2));
saveas(gcf, strcat(name1, format3));

NameArray = {'Marker'};
ValueArray = {'o'; 's'};
figure;
p = semilogy(N-1, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([1, 11]);
xlabel('N', 'interpreter', 'latex', 'fontsize', 18);
yticks([10.^(-10:2:1)]);
ylim([1e-11, 1e0]);
ylabel('Error', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name2, format2));
saveas(gcf, strcat(name2, format3));

disp('Lotka-Volterrra model completed!');

%% Simple pendulum
clear;
n = 200;
T = 20;
N = [3,5,7,9,11];
x0 = [- pi/4, pi/6];
r = [pi/8, pi/12];
frac = 0.2;
op = 1;

p = struct();
p.g = 9.8;
p.L = 9.8;

fmodel = @(t, x, p) SimplePendulum(t, x, p);
fdynamics = @(X, Y, p) SimplePendulumDynamics(X, Y, p);

% test
running_time = zeros(1, length(N));
error_terminal = zeros(2, length(N));

t2 = linspace(0, T, 20000);
xt = RungeKutta9(fmodel, x0, t2, p);
for i = 1:length(N)
    tic;
    [xnT, n_decomp] = ASK_2D(fdynamics, p, n, T, N(i), x0, r, frac, op);
    running_time(i) = toc;
    
    error_terminal(:, i) = abs(xt(:, end) - xnT);
end

% display & store
model_name = 'SP';
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name1 = strcat('results\GaussLobatto\time_', model_name, '_N_', num2str(T));
name2 = strcat('results\GaussLobatto\error_', model_name, '_N_', num2str(T));

save(strcat(name1, format1), 'running_time');
save(strcat(name2, format1), 'error_terminal');

figure;
plot(N, running_time, 'o-', 'linewidth', 2, 'markersize', 10);
xlim([2, 12]);
xlabel('N', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Running Time (s)', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name1, format2));
saveas(gcf, strcat(name1, format3));

NameArray = {'Marker'};
ValueArray = {'o'; 's'};
figure;
p = semilogy(N-1, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([1, 11]);
xlabel('N', 'interpreter', 'latex', 'fontsize', 18);
yticks([10.^(-14:2:1)]);
ylim([1e-15, 1e0]);
ylabel('Error', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name2, format2));
saveas(gcf, strcat(name2, format3));

disp('Simple pendulum model completed!');

%% Limit Cycle
clear;
n = 200;
T = 20;
N = [3,5,7,9,11];
x0 = [sqrt(2)/2, -sqrt(2)/2];
r = [sqrt(2)/6, sqrt(2)/6];
frac = 0.2;
op = 1;

p = struct();
p.k = 1;

fmodel = @(t, x, p) LimitCycle(t, x, p);
fdynamics = @(X, Y, p) LimitCycleDynamics(X, Y, p);

% test
running_time = zeros(1, length(N));
error_terminal = zeros(2, length(N));

xtT = [sin(T+pi/4); sin(T-pi/4)]; %closed-form
for i = 1:length(N)
    tic;
    [xnT, n_decomp] = ASK_2D(fdynamics, p, n, T, N(i), x0, r, frac, op);
    running_time(i) = toc;
    
    error_terminal(:, i) = abs(xtT - xnT);
end

% display & store
model_name = 'LC';
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name1 = strcat('results\GaussLobatto\time_', model_name, '_N_', num2str(T));
name2 = strcat('results\GaussLobatto\error_', model_name, '_N_', num2str(T));

save(strcat(name1, format1), 'running_time');
save(strcat(name2, format1), 'error_terminal');

figure;
plot(N, running_time, 'o-', 'linewidth', 2, 'markersize', 10);
xlim([2, 12]);
xlabel('N', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Running Time (s)', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name1, format2));
saveas(gcf, strcat(name1, format3));

NameArray = {'Marker'};
ValueArray = {'o'; 's'};
figure;
p = semilogy(N-1, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([1, 11]);
xlabel('N', 'interpreter', 'latex', 'fontsize', 18);
yticks([10.^(-10:2:-2)]);
ylim([1e-10, 1e-2]);
ylabel('Error', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name2, format2));
saveas(gcf, strcat(name2, format3));

disp('Limit cycle model completed!');
