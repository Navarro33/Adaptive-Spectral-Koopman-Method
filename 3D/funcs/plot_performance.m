function plot_performance(t, yt, yn, y_label, legend_first, file)

plot(t, yt, 'r-', 'LineWidth', 2);
hold on
plot(t, yn, 'b:', 'LineWidth', 2);
hold off

xlabel('t');
ylabel(y_label);
legend(legend_first, 'Koopman');

if nargin == 6
    saveas(gcf, file);
end



