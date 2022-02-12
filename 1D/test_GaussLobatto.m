%%% Test on Gauss-Lobatto Points %%%
clear; close all; clc;

addpath('funcs');
addpath('models');

%% Linear
n = 200;
T = 20;
N = [3,5,7,9,11];
x0 = 1;
r = 0.2;
frac = 0.2;
op = 1;

p = struct();
p.k = -0.5;

fdynamics = @(x, p) LinearDynamics(x, p);

% test
running_time = zeros(1, length(N));
error_terminal = zeros(1, length(N));

xtT = exp( p.k * T);
for i = 1:length(N)
    tic;
    [xnT, n_decomp] = ASK_1D(fdynamics, p, n, T, N(i), x0, r, frac, op);
    running_time(i) = toc;
    
    error_terminal(i) = abs(xtT - xnT);
end

% display & store
model_name = 'Linear';
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
ValueArray = {'o'};
figure;
p = semilogy(N-1, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([1, 11]);
xlabel('N', 'interpreter', 'latex', 'fontsize', 18);
yticks([10.^(-17:1:-15)]);
ylim([5e-18, 1e-15]);
ylabel('Error', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name2, format2));
saveas(gcf, strcat(name2, format3));

disp('Linear model completed!');

%% Cubic
clear;
n = 200;
T = 20;
N = [3,5,7,9,11];
x0 = 1;
r = 0.2;
frac = 0.2;
op = 1;

p = struct();
p.k = 10;

fdynamics = @(x, p) CubicDynamics(x, p);

% test
running_time = zeros(1, length(N));
error_terminal = zeros(1, length(N));

xtT = nthroot(p.k * T + 1, 3);
for i = 1:length(N)
    tic;
    [xnT, n_decomp] = ASK_1D(fdynamics, p, n, T, N(i), x0, r, frac, op);
    running_time(i) = toc;
    
    error_terminal(i) = abs(xtT - xnT);
end

% display & store
model_name = 'Cubic';
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
ValueArray = {'o'};
figure;
p = semilogy(N-1, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([1, 11]);
xlabel('N', 'interpreter', 'latex', 'fontsize', 18);
yticks([10.^(-10:2:-4)]);
ylim([1e-11, 1e-3]);
ylabel('Error', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name2, format2));
saveas(gcf, strcat(name2, format3));

disp('Cubic model completed!');

%% Cosine Square
clear;
n = 200;
T = 20;
N = [3,5,7,9,11];
x0 = pi/4;
r = pi/20;
frac = 0.2;
op = 1;

p = struct();
p.k = -0.5;

fdynamics = @(x, p) CosineSquareDynamics(x, p);

% test
running_time = zeros(1, length(N));
error_terminal = zeros(1, length(N));

xtT = atan(p.k * T + 1);
for i = 1:length(N)
    tic;
    [xnT, n_decomp] = ASK_1D(fdynamics, p, n, T, N(i), x0, r, frac, op);
    running_time(i) = toc;
    
    error_terminal(i) = abs(xtT - xnT);
end

% display & store
model_name = 'Cosine';
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
ValueArray = {'o'};
figure;
p = semilogy(N-1, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([1, 11]);
xlabel('N', 'interpreter', 'latex', 'fontsize', 18);
yticks([10.^(-8:2:-2)]);
ylim([1e-9, 1e-1]);
ylabel('Error', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name2, format2));
saveas(gcf, strcat(name2, format3));

disp('Cosine Square model completed!');
