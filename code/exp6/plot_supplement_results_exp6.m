%% initialization and settings
[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
wav_path = fullfile(file_path, sprintf('../../figures/%s/wavs/', exp_name)); % get path to raw waveforms
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)
save_path = fullfile(file_path, sprintf('../../figures/%s/', exp_name)); % get path to save figure
if ~exist(save_path, 'dir'); mkdir(save_path); end

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

% set colors
rep_color = [152 9 15];
alt_color = rep_color + (255 - rep_color) * 0.4;
clrs = [alt_color; rep_color] / 255;

%% plot results figure

% re-analyze data by splitting into first half versus second half of block
temp_block_position = repmat(1:2, 20, 1);
block_position = repmat(temp_block_position(:), 8, 1);
d_primes = nan(20,length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)
    responses = results(i).responses;
    responses_to_use = responses(:, 2:end);    
    responses_to_use(:, 2) = block_position;
    responses_to_use(sum(responses,2)==0,:) = []; % remove trials with no data                            

    [G, ID1, ID2] = findgroups(responses_to_use(:,2),responses_to_use(:,3)); % get the group info
    group_info = [ID1 ID2];
    hit_rates = splitapply(@compute_hit_rate, responses_to_use(:, 4), responses_to_use(:, 5), G);  % get hit rates and false alarm rates for all conditions (SNR by onset time)      
    mean_false_alarm_rates = splitapply(@compute_false_alarm_rate, responses_to_use(:, 4), responses_to_use(:, 5), responses_to_use(:, 2)); % compute false alarms over all "conditions" since they only apply with foregrounds present
    temp_false_alarm_rates = repmat(mean_false_alarm_rates, 1, 10)';
    d_primes(:, i) = norminv(hit_rates) - norminv(temp_false_alarm_rates(:));
    keep_subs(i) = results(i).pass;
end

d_primes = d_primes(:, logical(keep_subs));
n_subs = size(d_primes,2);

% fit elbow function
x = linspace(250,2500,10);
y_1 = mean(d_primes(1:10, :),2)';
y_2 = mean(d_primes(11:20, :),2)';
est_params_1 = fit_elbow_func(x,y_1);
est_params_2 = fit_elbow_func(x,y_2);
x_test = linspace(250, 2500, 1000);
y_pred_1 = elbow_func(x_test, est_params_1);
y_pred_2 = elbow_func(x_test, est_params_2);

% initialize figure
y_scale = 3 / 3.25; % adjust for size difference in y dimension since not plotting elbow point CI
figure('units','inches','innerposition',[0 0 9 11*y_scale])
hold on
timevec = linspace(250,2500,10);

% plot dummy variables for legend
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(1,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:)); % style of performance data
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(2,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(2,:)); % style of performance data

% plot shaded errorbars (SEM)
m = mean(d_primes(1:10, :),2);
s = std(d_primes(1:10, :),0,2)/sqrt(n_subs);
ub = m + s;
lb = m - s;
x_patch = [timevec, flip(timevec)];
y_patch = [ub(:,1)', flip(lb(:,1)')];
ptch = patch(x_patch, y_patch, clrs(1,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

m = mean(d_primes(11:20, :),2);
s = std(d_primes(11:20, :),0,2)/sqrt(n_subs);
ub = m + s;
lb = m - s;
x_patch = [timevec, flip(timevec)];
y_patch = [ub(:,1)', flip(lb(:,1)')];
ptch = patch(x_patch, y_patch, clrs(2,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

% plot elbow function fit to data
p1 = plot(x_test, y_pred_1, '--', 'LineWidth', 3, 'Color', clrs(1,:));
p2 = plot(x_test, y_pred_2, '--', 'LineWidth', 3, 'Color', clrs(2,:));

% plot human results
scatter(timevec, mean(d_primes(1:10, :),2), 125, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)
scatter(timevec, mean(d_primes(11:20, :),2), 125, 'MarkerEdgeColor','w', 'MarkerFaceColor', clrs(2,:), 'LineWidth', 0.75)

% plot deltas
delta_x = 2600;
plot([delta_x, delta_x], [y_pred_1(1) y_pred_1(end)], 'LineWidth', 2, 'Color', clrs(1,:))
plot([delta_x - 25, delta_x + 5], [y_pred_1(end) y_pred_1(end)], 'LineWidth', 2, 'Color', clrs(1,:))
plot([delta_x - 25, delta_x + 5], [y_pred_1(1) y_pred_1(1)], 'LineWidth', 2, 'Color', clrs(1,:))
delta_x = 2650;
plot([delta_x, delta_x], [y_pred_2(1) y_pred_2(end)], 'LineWidth', 2, 'Color', clrs(2,:))
plot([delta_x - 25, delta_x + 5], [y_pred_2(end) y_pred_2(end)], 'LineWidth', 2, 'Color', clrs(2,:))
plot([delta_x - 25, delta_x + 5], [y_pred_2(1) y_pred_2(1)], 'LineWidth', 2, 'Color', clrs(2,:))

% add subject numbers
text(75, 2.9, sprintf('n = %d', n_subs), ...
     'FontSize', 25, ...
     'FontAngle','italic', ...
     'Color','k')

% stylize figure and add labels
box off
ylabel('Detection Performance (d'')', 'Position', [-300, 1.5])
xlabel('Foreground Onset Time (ms)')
ylim([0 3])
yticks(0:0.5:3)
ytickformat('%.1f')
xlim([0 2750])
xticks((timevec))
set(gca,'linewidth',2)
set(gca, 'TickDir', 'both')
set(gca,'FontSize',25)
set(gca,'XTickLabelRotation',45)

% improve look of legend
[l, objh]  = legend({'    First Half', '    Second Half'},'Location','southeast');
objh(3).XData = [objh(3).XData(1), objh(3).XData(2)*1.18];
objh(3).Color = [objh(3).Color 0.5];
objh(4).XData = mean(objh(3).XData);
objh(4).LineWidth = 1;
objh(5).XData = [objh(5).XData(1), objh(5).XData(2)*1.18];
objh(5).Color = [objh(5).Color 0.5];
objh(6).XData = mean(objh(5).XData);
objh(6).LineWidth = 1;
l.Position(1) = l.Position(1)*1.075;
l.Position(2) = l.Position(2)*1.2;
l.EdgeColor = 'none';
text(2110, 0.5, 'Position in Block', ...
     'FontSize', 22.5, ...
     'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', ...
     'Color',ones(3,1)*0)

% save figure as svg
set(gca,'Position', [(0.145), (0.145/y_scale), (1 - 0.165), (1 - 0.195/y_scale)]) % remove space around edges
save_name = fullfile(save_path, sprintf('%s_supplement_results.svg', exp_name));
saveas(gcf, save_name)