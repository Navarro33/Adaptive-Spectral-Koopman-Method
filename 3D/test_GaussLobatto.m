%%% Test on Gauss-Lobatto Points %%%
clear; close all; clc;

addpath('funcs');
addpath('models');

%% Kraichnan-Orszag model
n = 400;
T = 20;
% n = 200;
% T = 10;
N = [3,5,7,9,11];
x0 = [1, 2, -3];
r = [1, 1, 1];
frac = 0.15;
op = 1;

p = struct();
p.k1 = 1;
p.k2 = 1;
p.k3 = 1;

fmodel = @(t, x, p) KraichnanOrszagModel(t, x, p);
fdynamics = @(X, Y, Z, p) KraichnanOrszagDynamics(X, Y, Z, p);

% test
running_time = zeros(1, length(N));
error_terminal = zeros(3, length(N));

t2 = linspace(0, T, 20000+1);
xt = RungeKutta9(fmodel, x0, t2, p);
for i = 1:length(N)
    tic;
    [xnT, n_decomp] = ASK_3D(fdynamics, p, n, T, N(i), x0, r, frac, op);
    running_time(i) = toc;
    
    error_terminal(:, i) = abs(xt(:, end) - xnT);
end

% display & store
model_name = 'KO';
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
ValueArray = {'o'; 's'; '^'};
figure;
p = semilogy(N(1:4)-1, error_terminal(:, 1:4)', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([1, 9]);
xlabel('N', 'interpreter', 'latex', 'fontsize', 18);
yticks(10.^(-8:2:0));
ylim([1e-8, 5e0]);
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
N = [3,5,7,9,11];
x0 = [5,5,5];
r = [4,4,4];
frac = 0.5;
op = 1;

p = struct();
p.alpha = 10;
p.beta = 28;
p.gamma = 3;

fmodel = @(t, x, p) LorenzModel(t, x, p);
fdynamics = @(X, Y, Z, p) LorenzDynamics(X, Y, Z, p);

% test
running_time = zeros(1, length(N));
error_terminal = zeros(3, length(N));

t2 = linspace(0, T, 20000+1);
xt = RungeKutta9(fmodel, x0, t2, p);
for i = 1:length(N)
    tic;
    [xnT, n_decomp] = ASK_3D(fdynamics, p, n, T, N(i), x0, r, frac, op);
    running_time(i) = toc;
    
    error_terminal(:, i) = abs(xt(:, end) - xnT);
end
%%
% display & store
model_name = 'LO';
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
ValueArray = {'o'; 's'; '^'};
figure;
p = semilogy(N(1:4)-1, error_terminal(:, 1:4)', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([1, 9]);
xlabel('N', 'interpreter', 'latex', 'fontsize', 18);
yticks(10.^(-8:2:0));
ylim([1e-8, 5e1]);
ylabel('Error', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name2, format2));
saveas(gcf, strcat(name2, format3));

disp('Lorenz model completed!');
