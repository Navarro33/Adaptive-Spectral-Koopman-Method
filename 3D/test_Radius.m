 %%% Test on Radius %%%
clear; close all; clc;

addpath('funcs');
addpath('models');

%% Kraichnan-Orszag model
n = 400;
T = 20;
% n = 200;
% T = 10;
N = 3;
x0 = [1, 2, -3];
R = 0.1:0.1:1;
frac = 0.15;
op = 1;

p = struct();
p.k1 = 1;
p.k2 = 1;
p.k3 = 1;

fmodel = @(t, x, p) KraichnanOrszagModel(t, x, p);
fdynamics = @(X, Y, Z, p) KraichnanOrszagDynamics(X, Y, Z, p);

% test
running_time = zeros(1, length(R));
error_terminal = zeros(3, length(R));

t2 = linspace(0, T, 20000);
xt = RungeKutta9(fmodel, x0, t2, p);
for i = 1:length(R)
    r = [R(i), R(i), R(i)];
    tic;
    [xnT, n_decomp] = ASK_3D(fdynamics, p, n, T, N, x0, r, frac, op);
    running_time(i) = toc;
    
    error_terminal(:, i) = abs(xt(:, end) - xnT);
end

% display & store
model_name = 'KO';
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name1 = strcat('results\Radius\time_', model_name, '_r_', num2str(T));
name2 = strcat('results\Radius\error_', model_name, '_r_', num2str(T));

save(strcat(name1, format1), 'running_time');
save(strcat(name2, format1), 'error_terminal');

figure;
plot(R, running_time, 'o-', 'linewidth', 2, 'markersize', 10);
xlim([0, 1.1]);
xlabel('r', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Running Time (s)', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name1, format2));
saveas(gcf, strcat(name1, format3));

NameArray = {'Marker'};
ValueArray = {'o'; 's'; '^'};
figure;
p = semilogy(R, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([0, 1.1]);
xlabel('r', 'interpreter', 'latex', 'fontsize', 18);
ylim([1e-6, 1e1]);
ylabel('Error', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name2, format2));
saveas(gcf, strcat(name2, format3));

disp('Kraichnan-Orszag model completed!');

%% Lorenz model
clear;
n = 500;
T = 10;
N = 5;
x0 = [5,5,5];
R = 0.5:0.1:1.5;
frac = 0.5;
op = 1;

p = struct();
p.alpha = 10;
p.beta = 28;
p.gamma = 3;

fmodel = @(t, x, p) LorenzModel(t, x, p);
fdynamics = @(X, Y, Z, p) LorenzDynamics(X, Y, Z, p);

% test
running_time = zeros(1, length(R));
error_terminal = zeros(3, length(R));

t2 = linspace(0, T, 20000);
xt = RungeKutta9(fmodel, x0, t2, p);
for i = 1:length(R)
    r = [R(i), R(i), R(i)];
    tic;
    [xnT, n_decomp] = ASK_3D(fdynamics, p, n, T, N, x0, r, frac, op);
    running_time(i) = toc;
    
    error_terminal(:, i) = abs(xt(:, end) - xnT);
end

% display & store
model_name = 'LO';
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name1 = strcat('results\Radius\time_', model_name, '_r_', num2str(T));
name2 = strcat('results\Radius\error_', model_name, '_r_', num2str(T));

save(strcat(name1, format1), 'running_time');
save(strcat(name2, format1), 'error_terminal');

figure;
plot(R, running_time, 'o-', 'linewidth', 2, 'markersize', 10);
xlim([0.4, 1.6]);
xlabel('r', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Running Time (s)', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name1, format2));
saveas(gcf, strcat(name1, format3));

NameArray = {'Marker'};
ValueArray = {'o'; 's'; '^'};
figure;
p = semilogy(R, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([0.4, 1.6]);
xlabel('r', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Error', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name2, format2));
saveas(gcf, strcat(name2, format3));

disp('Lorenz model completed!');
