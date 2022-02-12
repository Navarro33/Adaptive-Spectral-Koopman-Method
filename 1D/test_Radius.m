%%% Test on Radius %%%
clear; close all; clc;

addpath('funcs');
addpath('models');

%% Linear
n = 200;
T = 20;
N = 3;
x0 = 1;
r = 0.05:0.05:0.5;
frac = 0.2;
op = 1;

p = struct();
p.k = -0.5;

fdynamics = @(x, p) LinearDynamics(x, p);

% test
running_time = zeros(1, length(r));
error_terminal = zeros(1, length(r));

xtT = exp( p.k * T);
for i = 1:length(r)
    tic;
    [xnT, n_decomp] = ASK_1D(fdynamics, p, n, T, N, x0, r(i), frac, op);
    running_time(i) = toc;
    
    error_terminal(i) = abs(xtT - xnT);
end

% display & store
model_name = 'Linear';
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name1 = strcat('results\Radius\time_', model_name, '_r_', num2str(T));
name2 = strcat('results\Radius\error_', model_name, '_r_', num2str(T));

save(strcat(name1, format1), 'running_time');
save(strcat(name2, format1), 'error_terminal');

figure;
plot(r, running_time, 'o-', 'linewidth', 2, 'markersize', 10);
xlim([0, 0.55]);
xlabel('r', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Running Time (s)', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name1, format2));
saveas(gcf, strcat(name1, format3));

NameArray = {'Marker'};
ValueArray = {'o'};
figure;
p = semilogy(r, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([0, 0.55]);
xlabel('r', 'interpreter', 'latex', 'fontsize', 18);
ylim([1e-19, 1e-15]);
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
N = 7;
x0 = 1;
r = 0.05:0.05:0.5;
frac = 0.2;
op = 1;

p = struct();
p.k = 10;

fdynamics = @(x, p) CubicDynamics(x, p);

% test
running_time = zeros(1, length(r));
error_terminal = zeros(1, length(r));

xtT = nthroot(p.k * T + 1, 3);
for i = 1:length(r)
    tic;
    [xnT, n_decomp] = ASK_1D(fdynamics, p, n, T, N, x0, r(i), frac, op);
    running_time(i) = toc;
    
    error_terminal(i) = abs(xtT - xnT);
end

% display & store
model_name = 'Cubic';
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name1 = strcat('results\Radius\time_', model_name, '_r_', num2str(T));
name2 = strcat('results\Radius\error_', model_name, '_r_', num2str(T));

save(strcat(name1, format1), 'running_time');
save(strcat(name2, format1), 'error_terminal');

figure;
plot(r, running_time, 'o-', 'linewidth', 2, 'markersize', 10);
xlim([0, 0.55]);
xlabel('r', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Running Time (s)', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name1, format2));
saveas(gcf, strcat(name1, format3));

NameArray = {'Marker'};
ValueArray = {'o'};
figure;
p = semilogy(r, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([0, 0.55]);
xlabel('r', 'interpreter', 'latex', 'fontsize', 18);
ylim([1e-9, 1e-5]);
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
N = 9;
x0 = pi/4;
r = 0.01:0.01:0.1;
frac = 0.2;
op = 1;

p = struct();
p.k = -0.5;

fdynamics = @(x, p) CosineSquareDynamics(x, p);

% test
running_time = zeros(1, length(r));
error_terminal = zeros(1, length(r));

xtT = atan(p.k * T + 1);
for i = 1:length(r)
    tic;
    [xnT, n_decomp] = ASK_1D(fdynamics, p, n, T, N, x0, r(i), frac, op);
    running_time(i) = toc;
    
    error_terminal(i) = abs(xtT - xnT);
end

% display & store
model_name = 'Cosine';
format1 = '.mat';
format2 = '.png';
format3 = '.eps';
name1 = strcat('results\Radius\time_', model_name, '_r_', num2str(T));
name2 = strcat('results\Radius\error_', model_name, '_r_', num2str(T));

save(strcat(name1, format1), 'running_time');
save(strcat(name2, format1), 'error_terminal');

figure;
plot(r, running_time, 'o-', 'linewidth', 2, 'markersize', 10);
xlim([0, 0.11]);
xlabel('r', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Running Time (s)', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name1, format2));
saveas(gcf, strcat(name1, format3));

NameArray = {'Marker'};
ValueArray = {'o'};
figure;
p = semilogy(r, error_terminal', '-', 'linewidth', 2, 'markersize', 10);
set(p, NameArray, ValueArray);
xlim([0, 0.11]);
xlabel('r', 'interpreter', 'latex', 'fontsize', 18);
yticks(10.^(-15:2:-7));
ylim([1e-15, 1e-7]);
ylabel('Error', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 16);
saveas(gcf, strcat(name2, format2));
saveas(gcf, strcat(name2, format3));

disp('Cosine Square model completed!');
