%%% Test for Other Methods %%%
clear; close all; clc;

addpath('funcs');
addpath('models');

%% Parameters
m = 200;
T = 20;
tspan = linspace(0, T, m+1);
x0 = [sqrt(2)/2, -sqrt(2)/2];

p.k = 1;

fmodel = @(t, x, p) LimitCycle(t, x, p);
xr = [sin(tspan + pi/4); sin(tspan - pi/4)];    % reference
xt_Euler = Euler(fmodel, x0, tspan, p);
xt_AB5 = AdamsBashforth5(fmodel, x0, tspan, p);
xt_AM4 = AdamsMoulton4(fmodel, x0, tspan, p);

errors_Euler = abs(xt_Euler - xr);
errors_AB5 = abs(xt_AB5 - xr);
errors_AM4 = abs(xt_AM4 - xr);


% % performance x1
% figure;
% plot(tspan, xt_Euler(1, :), 'o', 'MarkerSize', 7);
% hold on
% plot(tspan, xt_AB5(1, :), '^', 'MarkerSize', 7);
% plot(tspan, xt_AM4(1, :), 's', 'MarkerSize', 7);
% plot(tspan, xr(1, :), '--', 'LineWidth', 2);
% hold off
% xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
% ylim([-1.5, 1.5]);
% ylabel('x1', 'Interpreter', 'latex', 'FontSize', 18);
% set(gca, 'LineWidth', 2);
% set(gca, 'FontSize', 16);
% 
% 
% % performance x2
% figure;
% plot(tspan, xt_Euler(2, :), 'o', 'MarkerSize', 7);
% hold on
% plot(tspan, xt_AB5(2, :), '^', 'MarkerSize', 7);
% plot(tspan, xt_AM4(2, :), 's', 'MarkerSize', 7);
% plot(tspan, xr(2, :), '--', 'LineWidth', 2);
% hold off
% xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
% ylim([-1.5, 1.5]);
% ylabel('x1', 'Interpreter', 'latex', 'FontSize', 18);
% set(gca, 'LineWidth', 2);
% set(gca, 'FontSize', 16);

% error x1
figure;
semilogy(tspan, errors_Euler(1, :), 'o', 'MarkerSize', 5);
hold on
semilogy(tspan, errors_AB5(1, :), '^', 'MarkerSize', 5);
semilogy(tspan, errors_AM4(1, :), 's', 'MarkerSize', 5);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylim([1e-12, 1e-2]);
yticks(10.^(-12:2:-2));
ylabel('x1 Error', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);

% error x2
figure;
semilogy(tspan, errors_Euler(2, :), 'o', 'MarkerSize', 5);
hold on
semilogy(tspan, errors_AB5(2, :), '^', 'MarkerSize', 5);
semilogy(tspan, errors_AM4(2, :), 's', 'MarkerSize', 5);
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 18);
ylim([1e-12, 1e-2]);
yticks(10.^(-12:2:-2));
ylabel('x2 Error', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 16);
