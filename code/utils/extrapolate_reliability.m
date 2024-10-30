function target_sample_size = extrapolate_reliability(sample_sizes, ...
                                                      reliabilities, ...
                                                      target_reliability, ... 0.9
                                                      max_sample_size, ... 
                                                      param_init ... [-10 10 1]
                                                      )

x = sample_sizes;
y = reliabilities;

power_func = @(x, params) params(1)./(x + params(2)) + params(3);
cost_func = @(params) mean((y - power_func(x, params)).^2);

options = optimset('MaxIter',500*3, 'MaxFunEval', 500*3);
est_params = fminsearch(cost_func, param_init, options);
x_test = sample_sizes(1):max_sample_size;
predicted_reliability = power_func(x_test, est_params);

figure
hold on
plot([x_test(1), x_test(end)], [target_reliability, target_reliability], '--r', 'LineWidth',2)
plot(x_test, round(predicted_reliability,2), 'b', 'LineWidth',2)
scatter(x, y, 50, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
grid on
box off
ylim([0, 1])
% xticks(x_test(1):10:x_test(end))
ylabel('Split-Half Reliability')
xlabel('Number of Subjects')

target_sample_size = x_test(find(round(predicted_reliability,2) >= target_reliability, 1));