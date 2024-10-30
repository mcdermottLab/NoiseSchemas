%% initialization and settings
[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
wav_path = fullfile(file_path, sprintf('../../figures/%s/wavs/', exp_name)); % get path to raw waveforms
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)
comparison_exp = 'exp6';
matched_permutation_test_file = fullfile(file_path, sprintf('../../results/%s/matched_permutation_test_%s_vs_%s.mat', exp_name, exp_name, comparison_exp));
load(matched_permutation_test_file)
indices = [6, 8, 10]; % corresponds to timepoints 1500, 2000, 2500 ms
test_indices = setdiff(1:10, indices);
save_path = fullfile(file_path, sprintf('../../figures/%s/', exp_name)); % get path to save figure
if ~exist(save_path, 'dir'); mkdir(save_path); end

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

% set colors
rep_color = [152 9 15];
nonrep_color = [153 153 153]; % gray to match backgrounds
clrs = [rep_color; nonrep_color] / 255;

% settings for waveforms...note must match saved waveforms in wav_path
nstims = 5;
stims_bg = [123 99 32 40 22]; % stim IDs for backgrounds
stims_fg = [32 22 40 123 99]; % stim IDs for foregrounds
outliers = [2 2 1 1 2]; % presence of outlier
times = [5 1 2 3 9]; % foreground onset times
snrs = [1 1 1 1 1]; % SNRs
trials = [1 2 3 319 320]; % trial numbers

% waveform plot settings
clip_thresh = 0.25; % value at which to clip waveforms to make them look a bit nicer
y_positions = flip((1:(nstims+1)))*clip_thresh*2; % vector of y-positions of waveforms
skip_position = y_positions(end-2); % position that is skipped and replaced with ellipsis
y_positions(end-2) = []; % actually remove the skipped position

% stimulus settings (taken from stimulus creation script)
seglen = 500; % of the foreground (in MS)
bg_length = 3.25; % length of background in seconds
audio_sr = 20000; % sample rate of all the audio
SNRvec = -8; % set of possible SNRs
onset_times = linspace(0.25,2.50,10); % set of possible positions of the outlier

%% plot diagram

% initialize diagram figure
figure('units','inches','innerposition',[0 0 7 8])
ylim([0.25, 4])
xlim([-0.5*audio_sr, bg_length*audio_sr + 0.1*audio_sr])
axis off
hold on
% loop over trials
for n = 1:nstims
    % load background waveform    
    bacgkround_file = sprintf('background_%d.wav', stims_bg(n));
    background = audioread(fullfile(wav_path, bacgkround_file));
    background(abs(background) > clip_thresh) = clip_thresh;

    % load foreground waveform
    foreground_file = sprintf('foreground_%d.wav', stims_fg(n));
    foreground = audioread(fullfile(wav_path, foreground_file));

    % get proper onset_time
    t = times(n);
    ind1 = onset_times(t)*audio_sr; % convert the ms timing to a sample

    % set SNR of foreground relative to background
    snr_ind = snrs(n);
    clip1 = background(ind1:ind1+(audio_sr*seglen/1000)-1); % extract corresponding segment of background            
    setsnr = SNRvec(snr_ind); % get SNR level
    c = (rms(clip1)/rms(foreground))*10^(setsnr/20); % compute the scaling factor to get the proper SNR
    outlier_switch = outliers(n);
    if outlier_switch == 1 % then no outlier
        c = 0; % will kill the outlier
    else % otherwise keep the right SNR level
    end
    foreground = c*foreground; % actually set proper SNR
    foreground(abs(foreground) > clip_thresh) = clip_thresh; % apply clip threshold

    % plot stimulus
    plot(background + y_positions(n), 'LineWidth', 1.5, 'Color', 0.6*ones(3,1))
    if outlier_switch==2 % then also plot the foreground
        temp_timevec = 1:length(background);
        fg_timevec = temp_timevec(ind1:ind1+(audio_sr*seglen/1000)-1);
        plot(fg_timevec, foreground + y_positions(n), 'LineWidth', 1.5, 'Color', 'k')
    end
    x_offset = 0;
    if n <= 3; x_offset = 0.025*audio_sr; end
    % add trial number label
    text(-0.25*audio_sr + x_offset, y_positions(n), sprintf('%d', trials(n)), 'HorizontalAlignment','center', 'FontSize', 18, 'VerticalAlignment','middle')

end

% add more text to figure
text(-0.25*audio_sr + 0.025*audio_sr, skip_position, '. . .', 'HorizontalAlignment','center', 'FontSize', 24, 'Rotation', 90, 'VerticalAlignment','baseline')
text(length(background)/2, skip_position, '. . .', 'HorizontalAlignment','center', 'FontSize', 24, 'Rotation', 90, 'VerticalAlignment','baseline')
text(-0.25*audio_sr + 0.025*audio_sr, 3.4, 'trial', 'HorizontalAlignment','center', 'FontSize', 18, 'VerticalAlignment','middle')
text(length(background)/2, 3.4, 'stimulus', 'HorizontalAlignment','center', 'FontSize', 18, 'VerticalAlignment','middle')
text(mean([-0.5*audio_sr, bg_length*audio_sr + 0.1*audio_sr]), 3.75, sprintf('Did you hear\none sound or two sounds?'), ...
    'HorizontalAlignment','center', 'FontSize', 24, 'VerticalAlignment','middle', 'FontWeight','bold', 'FontAngle','italic', 'Color', clrs(2,:))

% save diagram figure
set(gca,'LooseInset',get(gca,'TightInset')); % remove space around edges
save_name = fullfile(save_path, sprintf('%s_diagram.pdf', exp_name));
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');

%% plot results figure

% get data
d_primes = nan(length(results(1).d_primes),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)
    d_primes(:,i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end

d_primes = d_primes(:, logical(keep_subs));
n_subs = size(d_primes,2);

d_primes_nonrep = d_primes;
n_nonrep = n_subs;

prev_exp_name = sprintf('exp%d', str2num(exp_name(4:end))-1);
prev_results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', prev_exp_name)); % get path to results
load(prev_results_file)

d_primes_rep = nan(length(results(1).d_primes),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)    
    d_primes_rep(:,i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end

d_primes_rep = d_primes_rep(:, logical(keep_subs));
n_rep = size(d_primes_rep,2);

% fit elbow function
x = linspace(250,2500,10);
y_rep = mean(d_primes_rep,2)';
y_nonrep = mean(d_primes_nonrep,2)';
est_params_rep = fit_elbow_func(x,y_rep);
est_params_nonrep = fit_elbow_func(x,y_nonrep);
x_test = linspace(250, 2500, 1000);
y_pred_rep = elbow_func(x_test, est_params_rep);
y_pred_nonrep = elbow_func(x_test, est_params_nonrep);

% initialize figure
y_scale = 3 / 3.25; % adjust for size difference in y dimension since not plotting elbow point CI
figure('units','inches','innerposition',[0 0 9 11*y_scale])
hold on
timevec = linspace(250,2500,10);

% plot dummy variables for legend
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(1,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:)); % style of performance data
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(2,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(2,:)); % style of performance data

% plot shaded errorbars (SEM)
m = mean(d_primes_rep,2);
s = std(d_primes_rep,0,2)/sqrt(n_rep);
ub = m + s;
lb = m - s;
x_patch = [timevec, flip(timevec)];
y_patch = [ub(:,1)', flip(lb(:,1)')];
ptch = patch(x_patch, y_patch, clrs(1,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

m = mean(d_primes_nonrep,2);
s = std(d_primes_nonrep,0,2)/sqrt(n_nonrep);
ub = m + s;
lb = m - s;
x_patch = [timevec, flip(timevec)];
y_patch = [ub(:,1)', flip(lb(:,1)')];
ptch = patch(x_patch, y_patch, clrs(2,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

% plot elbow function fit to data
p1 = plot(x_test, y_pred_rep, '--', 'LineWidth', 3, 'Color', clrs(1,:));
p2 = plot(x_test, y_pred_nonrep, '--', 'LineWidth', 3, 'Color', clrs(2,:));

% plot human results
scatter(timevec, mean(d_primes_rep,2), 125, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)
scatter(timevec, mean(d_primes_nonrep,2), 125, 'MarkerEdgeColor','w', 'MarkerFaceColor', clrs(2,:), 'LineWidth', 0.75)

% plot deltas
delta_x = 2600;
plot([delta_x, delta_x], [y_pred_rep(1) y_pred_rep(end)], 'LineWidth', 2, 'Color', clrs(1,:))
plot([delta_x - 25, delta_x + 5], [y_pred_rep(end) y_pred_rep(end)], 'LineWidth', 2, 'Color', clrs(1,:))
plot([delta_x - 25, delta_x + 5], [y_pred_rep(1) y_pred_rep(1)], 'LineWidth', 2, 'Color', clrs(1,:))
plot([delta_x, delta_x], [y_pred_nonrep(1) y_pred_nonrep(end)], 'LineWidth', 2, 'Color', clrs(2,:))
plot([delta_x - 25, delta_x + 5], [y_pred_nonrep(end) y_pred_nonrep(end)], 'LineWidth', 2, 'Color', clrs(2,:))
plot([delta_x - 25, delta_x + 5], [y_pred_nonrep(1) y_pred_nonrep(1)], 'LineWidth', 2, 'Color', clrs(2,:))

% add subject numbers
text(75, 2.9, sprintf('n = %d', n_rep), ...
     'FontSize', 30, ...
     'FontAngle','italic', ...
     'Color',clrs(1,:))

text(75, 2.75, sprintf('n = %d', n_nonrep), ...
     'FontSize', 30, ...
     'FontAngle','italic', ...
     'Color',clrs(2,:))

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
set(gca,'FontSize',30)
set(gca,'XTickLabelRotation',45)

% improve look of legend
[l, objh]  = legend({'  Repeated (Exp. 6)', '  Non-Repeated (Exp. 7)'},'Location','northeast');
objh(3).XData = [objh(3).XData(1), objh(3).XData(2)*1.18];
objh(3).Color = [objh(3).Color 0.5];
objh(4).XData = mean(objh(3).XData);
objh(4).LineWidth = 1;
objh(5).XData = [objh(5).XData(1), objh(5).XData(2)*1.18];
objh(5).Color = [objh(5).Color 0.5];
objh(6).XData = mean(objh(5).XData);
objh(6).LineWidth = 1;
l.EdgeColor = 'none';
text(1035, 2.9, 'Background Type', ...
     'FontSize', 28, ...
     'FontWeight', 'bold', ...
     'HorizontalAlignment', 'left', ...
     'Color',ones(3,1)*0)

% reposition
set(gca,'Position', [(0.145), (0.165/y_scale), (1 - 0.165), (1 - 0.195/y_scale)]) % remove space around edges
l.Position(2) = l.Position(2)*0.955;

% save figure as svg
save_name = fullfile(save_path, sprintf('%s_results.pdf', exp_name));
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');

%% plot matched performance inset

d_primes_rep = d_primes_rep(:, matched_subs);
d_primes_nonrep = d_primes_nonrep(:, matched_subs_comparison);

% fit elbow function
x = linspace(250,2500,10);
x = x(test_indices);
y_rep = mean(d_primes_rep(test_indices,:),2,'omitnan')';
y_nonrep = mean(d_primes_nonrep(test_indices,:),2,'omitnan')';
est_params_rep = fit_elbow_func(x,y_rep);
est_params_nonrep = fit_elbow_func(x,y_nonrep);
x_test = linspace(250, 2500, 1000);
y_pred_rep = elbow_func(x_test, est_params_rep);
y_pred_nonrep = elbow_func(x_test, est_params_nonrep);

% initialize figure
y_scale = 1.2 / 3.25; % adjust for size difference in y dimension since not plotting elbow point CI
figure('units','inches','innerposition',[0 0 9 11*y_scale])
hold on
timevec = linspace(250,2500,10);

% plot shaded errorbars (SEM)
m = mean(d_primes_rep(test_indices,:),2,'omitnan');
s = std(d_primes_rep(test_indices,:),0,2,'omitnan')/sqrt(n_subs);
ub = m + s;
lb = m - s;
x_patch = [timevec(test_indices), flip(timevec(test_indices))];
y_patch = [ub(:,1)', flip(lb(:,1)')];
ptch = patch(x_patch, y_patch, clrs(1,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

m = mean(d_primes_nonrep(test_indices,:),2,'omitnan');
s = std(d_primes_nonrep(test_indices,:),0,2,'omitnan')/sqrt(n_subs);
ub = m + s;
lb = m - s;
x_patch = [timevec(test_indices), flip(timevec(test_indices))];
y_patch = [ub(:,1)', flip(lb(:,1)')];
ptch = patch(x_patch, y_patch, clrs(2,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

% plot elbow function fit to data
p1 = plot(x_test, y_pred_rep, '--', 'LineWidth', 3, 'Color', clrs(1,:));
p2 = plot(x_test, y_pred_nonrep, '--', 'LineWidth', 3, 'Color', clrs(2,:));

% plot human results
scatter(timevec(test_indices), mean(d_primes_rep(test_indices,:),2,'omitnan'), 125, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)
scatter(timevec(test_indices), mean(d_primes_nonrep(test_indices,:),2,'omitnan'), 125, 'MarkerEdgeColor','w', 'MarkerFaceColor', clrs(2,:), 'LineWidth', 0.75)

% plot deltas
delta_x = 2650;
plot([delta_x, delta_x], [y_pred_rep(1) y_pred_rep(end)], 'LineWidth', 2, 'Color', clrs(1,:))
plot([delta_x - 25, delta_x + 5], [y_pred_rep(end) y_pred_rep(end)], 'LineWidth', 2, 'Color', clrs(1,:))
plot([delta_x - 25, delta_x + 5], [y_pred_rep(1) y_pred_rep(1)], 'LineWidth', 2, 'Color', clrs(1,:))
delta_x = 2600;
plot([delta_x, delta_x], [y_pred_nonrep(1) y_pred_nonrep(end)], 'LineWidth', 2, 'Color', clrs(2,:))
plot([delta_x - 25, delta_x + 5], [y_pred_nonrep(end) y_pred_nonrep(end)], 'LineWidth', 2, 'Color', clrs(2,:))
plot([delta_x - 25, delta_x + 5], [y_pred_nonrep(1) y_pred_nonrep(1)], 'LineWidth', 2, 'Color', clrs(2,:))

% add text
text(1750, 1.6, sprintf('Performance\nmatched'), ...
     'FontSize', 32, ...
     'FontAngle', 'italic', ...
     'HorizontalAlignment', 'center', ...
     'Color',ones(3,1)*0)

% stylize figure and add labels
box off
ylabel('d''', 'Position', [-300, 1.75])
xlabel('Foreground Onset Time (ms)')
ylim([1.25 2.1])
yticks([1.25, 1.5, 2])
yticklabels({'0', '1.5', '2.0'})
% ytickformat('%.1f')
xlim([0 2750])
xticks((timevec(test_indices)))
set(gca,'linewidth',2)
set(gca, 'TickDir', 'both')
set(gca,'FontSize',30)
set(gca,'XTickLabelRotation',45)

% trim edges
set(gca,'Position', [(0.145), (0.165/y_scale), (1 - 0.165), (1 - 0.195/y_scale)]) % remove space around edges

% add axis break
ax1 = gca;
ax2 = axes('Position', ax1.Position);
px=[-50 50];
py1=[1.275 1.325]+0.05;
height=0.04;
py2=py1+height;
hold on
fill([px flip(px)],[py1 flip(py2)],'w','EdgeColor','w');
plot(px,py1,'k','LineWidth',2)
plot(px,py2,'k','LineWidth',2)
ax2.YLim = ax1.YLim;
ax2.XLim = ax1.XLim;
ax2.Position(1) = ax2.Position(1)*0.945;
axis off
set(gcf, 'InvertHardCopy', 'off'); 
set(gcf, 'Color', 'w');

% save figure as svg
save_name = fullfile(save_path, sprintf('%s_matched_results.pdf', exp_name));
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');