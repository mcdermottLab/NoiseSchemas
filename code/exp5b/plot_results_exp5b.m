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

% set random number generator
rng(1)

% set colors
clrs = [219, 133, 55; 239, 206, 89]/255;

% settings for waveforms...note must match saved waveforms in wav_path
nstims = 5;
stims = [123 99 32 40 22]; % stim IDs (bg/fg pairing)
outliers = [1 2 2 1 2]; % presence of outlier
times = [3 1 8 5 4]; % foreground onset times
snrs = [1 1 1 1 1]; % SNRs
ints = [2 1 2 1 1]; % interrupter types
trials = [1 2 3 319 320]; % trial numbers

% waveform plot settings
wav_scale = 2; % scales amplitude of everything for nicer plotting
clip_thresh = 0.25; % value at which to clip waveforms to make them look a bit nicer
y_positions = flip((1:(nstims+1)))*clip_thresh*2; % vector of y-positions of waveforms
skip_position = y_positions(end-2); % position that is skipped and replaced with ellipsis
y_positions(end-2) = []; % actually remove the skipped position

% stimulus settings (taken from stimulus creation script)
audio_sr = 20000; % sample rate of all the audio
seglen = 500; % of the foreground (in MS)
bg_length = 5; % length of background in seconds
int_length = 1.5; % length of noise in (s)
int_start = ((bg_length - int_length)/2)*audio_sr + 1; % starting position of interrupter
int_end = int_start + int_length*audio_sr - 1; % ending position of interrupter
SNRvec = -2; % set of possible SNRs
noise_dB = 12; % level of noise relative to background
desired_rms = .05 / (10^(noise_dB/20)); % sets noise to 0.05 rms and rest of texture lower
onset_times = [linspace(0.25,1,4), linspace(3.5,4.25,4)]; % set of possible positions of the outlier

d_primes = nan(length(results(1).d_primes),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)    
    d_primes(:,i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end

d_primes = d_primes(:, logical(keep_subs));
n_subs = size(d_primes,2);

% initialize figure
x_scale = 4.5 / 2.75;
y_scale = 3 / 3.25;
figure('units','inches','innerposition',[0 0 6.3*x_scale 11*y_scale])
hold on
timevec = [linspace(250,1000,4), linspace(3500,4250,4)]; % scale so plot matches previous figure scale
plotting_timevec = [timevec, timevec];
interrupter_times = [1750, 3250];

% plot dummy stuff for legend styling
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(1,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:)); % style of performance data
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(2,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(2,:)); % style of performance data

% plot shaded error bars & average performance
inds = [1:4, 9:12];
m = reshape(mean(d_primes(inds,:),2),4,2);
s = reshape(std(d_primes(inds,:),0,2),4,2)/sqrt(n_subs);
ub = m + s;
lb = m - s;
x_patch = [timevec(1:4), flip(timevec(1:4))];
y_patch1 = [ub(:,1)', flip(lb(:,1)')];
y_patch2 = [ub(:,2)', flip(lb(:,2)')];
ptch = patch(x_patch, y_patch1, clrs(1,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';
ptch = patch(x_patch, y_patch2, clrs(2,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

inds = [5:8, 13:16];
m = reshape(mean(d_primes(inds,:),2),4,2);
s = reshape(std(d_primes(inds,:),0,2),4,2)/sqrt(n_subs);
ub = m + s;
lb = m - s;
x_patch = [timevec(5:8), flip(timevec(5:8))];
y_patch1 = [ub(:,1)', flip(lb(:,1)')];
y_patch2 = [ub(:,2)', flip(lb(:,2)')];
ptch = patch(x_patch, y_patch1, clrs(1,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';
ptch = patch(x_patch, y_patch2, clrs(2,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

% plot human results
scatter(timevec, mean(d_primes(1:8,:),2), 125, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)
scatter(timevec, mean(d_primes(9:16,:),2), 125, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(2,:), 'LineWidth', 0.75)

% plot interrupter as shaded region
patch([interrupter_times(1), interrupter_times(1), interrupter_times(2), interrupter_times(2)], [0, 3, 3, 0], ones(1,3)*0.85, ...
    'FaceColor', ones(1,3)*0.85, 'FaceAlpha', 0.25, 'EdgeColor', 'None')

% text(mean(interrupter_times), 1.5, '1500 ms interrupter', ...
%      'FontSize', 30, ...
%      'Rotation', 90, ...
%      'Color',ones(3,1)*0.6, ...
%      'HorizontalAlignment','center')

% add subject numbers
text(75, 2.9, sprintf('n = %d', n_subs), ...
     'FontSize', 30, ...
     'FontAngle','italic', ...
     'Color',ones(3,1)*0)

% stylize figure and add labels
box off
ylim([0 3])
yticks(0:0.5:3)
ytickformat('%.1f')
ylabel('Detection Performance (d'')', 'Position', [-450, 1.5])
xlim([0 4500])
xticks(timevec)
xlabel('Foreground Onset Time (ms)')
set(gca,'linewidth',2)
set(gca, 'TickDir', 'both')
set(gca,'FontSize',30)
set(gca,'XTickLabelRotation',45)

% add legend
[l, objh]  = legend({'  Silence', '  White Noise'},'Location','southwest');
objh(3).XData = [objh(3).XData(1), objh(3).XData(2)*1.18];
objh(3).Color = [objh(3).Color 0.5];
objh(4).XData = mean(objh(3).XData);
objh(4).LineWidth = 1;
objh(5).XData = [objh(5).XData(1), objh(5).XData(2)*1.18];
objh(5).Color = [objh(5).Color 0.5];
objh(6).XData = mean(objh(5).XData);
objh(6).LineWidth = 1;
l.Position(1) = l.Position(1)*1.075;
l.Position(2) = l.Position(2)*1.05;
l.EdgeColor = 'none';

% legend title
text(100, 0.5, 'Interrupter Type', ...
     'FontSize', 28, ...
     'FontWeight', 'bold', ...
     'HorizontalAlignment', 'left', ...
     'Color',ones(3,1)*0)

% save figure as svg
set(gca,'Position', [(0.22/x_scale), (0.165/y_scale), (1 - 0.23/x_scale), (1 - 0.195/y_scale)]) % remove space around edges
save_name = fullfile(save_path, sprintf('%s_results.pdf', exp_name));
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');