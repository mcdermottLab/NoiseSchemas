%% adaptation account

rng(10)

close all
left_lim = -10;
right_lim = 20;
n_samples = 100;
offset = 0;
slope = 2;
y = -normcdf(linspace(left_lim,right_lim,n_samples), offset, slope);
noise = randn(size(y))/10;
signals(1,:) = smooth(y+noise, 5);

left_lim = -10;
right_lim = 20;
n_samples = 100;
offset = 7;
slope = 2;
max_scale = 1;
y = normpdf(linspace(left_lim,right_lim,n_samples), offset, slope);
y = (max_scale * y / max(y)) - 0.9;
noise = randn(size(y))/10;
signals(2,:) = smooth(y+noise, 5);

left_lim = -10;
right_lim = 20;
n_samples = 100;
offset = 4;
slope = 5;
y = -normcdf(linspace(left_lim,right_lim,n_samples), offset, slope);
noise = randn(size(y))/10;
signals(3,:) = smooth(y+noise, 5);

orig_signals = repmat(zeros(size(y)), 3, 1);
orig_signals = orig_signals + randn(size(orig_signals))/10;
orig_signals(2,:) = NaN;


f = 1/2;
dur = 1;
sr = 1000;
t = linspace(0, dur, dur*sr);
y = sin(2*pi*f*t);

[background, ~] = audioread('~/Documents/om/projects/HearingInNoise/hearing_in_noise/figures/exp1/wavs/background_32.wav');
background = background / abs(max(background));
[foreground, ~] = audioread('~/Documents/om/projects/HearingInNoise/hearing_in_noise/figures/exp1/wavs/foreground_32.wav');
foreground = foreground / rms(foreground) * 10^(-2/20)*rms(background);
background = background+4;
foreground = foreground+4;

figure('units','inches','innerposition',[0 0 7 3])
hold on

x_end = 0;
filt_dur = 24;
x_shift = 5;
x_ax_length = (filt_dur+1)*4 + 2*x_shift;
y_shift = 0.1;
gap = 50;
signal_length = 150;
foreground_length = signal_length / 3.25 * 0.5;
foreground_position = 1.25;
clrs = [0 72 158; 0 121 250; 0 194 249] / 255;

x_start = x_end + 1;
x_end = x_start + signal_length - 1;
plot(linspace(x_start, x_end, length(background)), background, 'Color', 0.6*ones(1,3))
foreground_start = x_start + signal_length / 3.25 * foreground_position;
foreground_end = foreground_start + foreground_length;
plot(linspace(foreground_start, foreground_end, length(foreground)), foreground, 'Color', 'k')

x_start = x_end + gap;
x_end = x_start + filt_dur;
y_pos = 6;
x_ax_start = x_start - x_shift;
x_ax_end = x_ax_start + x_ax_length;
plot(linspace(x_start, x_end, sr),y + y_pos + y_shift, 'LineWidth', 4, 'Color', clrs(1,:))
plot([x_ax_start,x_ax_end], [y_pos y_pos], 'k', 'LineWidth', 2)
plot([x_ax_start, x_ax_start], [y_pos, y_pos+1+(2*y_shift)], 'k', 'LineWidth', 2)

x_start = x_end + 1;
x_end = x_start + filt_dur;
y_pos = 4;
plot(linspace(x_start, x_end, sr),y + y_pos + y_shift, 'LineWidth', 4, 'Color', clrs(2,:))
plot([x_ax_start,x_ax_end], [y_pos y_pos], 'k', 'LineWidth', 2)
plot([x_ax_start, x_ax_start], [y_pos, y_pos+1+(2*y_shift)], 'k', 'LineWidth', 2)

x_start = x_end + 1;
x_end = x_start + filt_dur;
x_start = x_end + 1;
x_end = x_start + filt_dur;
y_pos = 1;
plot(linspace(x_start, x_end, sr),y + y_pos + y_shift, 'LineWidth', 4, 'Color', clrs(3,:))
plot([x_ax_start,x_ax_end], [y_pos y_pos], 'k', 'LineWidth', 2)
plot([x_ax_start, x_ax_start], [y_pos, y_pos+1+(2*y_shift)], 'k', 'LineWidth', 2)
text(mean([x_ax_start, x_ax_end]), 3, '. . .', 'HorizontalAlignment','center', 'FontSize', 20, 'Rotation', 90, 'FontWeight','bold')


x_start = x_end + gap;
x_end = x_start + signal_length - 1;

plot(linspace(x_start, x_end, length(orig_signals)), smooth(orig_signals(1,:), 5) + 2, 'LineWidth', 4, 'Color', 0.85*ones(1,3))    
plot(linspace(x_start, x_end, length(orig_signals)), signals(1,:) + 2, 'LineWidth', 4, 'Color', clrs(3,:))

plot(linspace(x_start, x_end, length(orig_signals)), smooth(orig_signals(2,:), 5) + 5, 'LineWidth', 4, 'Color', 0.85*ones(1,3))    
plot(linspace(x_start, x_end, length(orig_signals)), signals(2,:) + 5, 'LineWidth', 4, 'Color', clrs(2,:))

plot(linspace(x_start, x_end, length(orig_signals)), smooth(orig_signals(3,:), 5) + 7, 'LineWidth', 4, 'Color', 0.85*ones(1,3))    
plot(linspace(x_start, x_end, length(orig_signals)), signals(3,:) + 7, 'LineWidth', 4, 'Color', clrs(1,:))

text(mean([x_start, x_end]), 3, '. . .', 'HorizontalAlignment','center', 'FontSize', 20, 'Rotation', 90, 'FontWeight','bold')

xlim([-x_shift, x_end+x_shift])
axis off
ax = gca;
ax.Position = [0 -0.1 1 1.2];

save_path = '~/Documents/om/projects/HearingInNoise/hearing_in_noise/figures/fig1';
save_name = fullfile(save_path, 'adaptation_diagram.pdf');
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');

%% fixed filter account

rng(10)

n_samples = 100;
y = -1*ones(1, n_samples);
noise = randn(size(y))/10;
signals(1,:) = smooth(y+noise, 5);

left_lim = -10;
right_lim = 20;
n_samples = 100;
offset = 7;
slope = 2;
max_scale = 1;
y = normpdf(linspace(left_lim,right_lim,n_samples), offset, slope);
y = (max_scale * y / max(y)) - 0.9;
noise = randn(size(y))/10;
signals(2,:) = smooth(y+noise, 5);

y = -0.75*ones(1, n_samples);
noise = randn(size(y))/10;
signals(3,:) = smooth(y+noise, 5);

orig_signals = repmat(zeros(size(y)), 3, 1);
orig_signals = orig_signals + randn(size(orig_signals))/10;
orig_signals(2,:) = NaN;

f = 1/2;
dur = 1;
sr = 1000;
t = linspace(0, dur, dur*sr);
y = sin(2*pi*f*t);

figure('units','inches','innerposition',[0 0 7 3])
hold on

x_end = 0;
filt_dur = 24;
x_shift = 5;
x_ax_length = (filt_dur+1)*4 + 2*x_shift;
y_shift = 0.1;
gap = 50;
signal_length = 150;
foreground_length = signal_length / 3.25 * 0.5;
foreground_position = 1.25;

x_start = x_end + 1;
x_end = x_start + signal_length - 1;
plot(linspace(x_start, x_end, length(background)), background, 'Color', 0.6*ones(1,3))
foreground_start = x_start + signal_length / 3.25 * foreground_position;
foreground_end = foreground_start + foreground_length;
plot(linspace(foreground_start, foreground_end, length(foreground)), foreground, 'Color', 'k')

x_start = x_end + gap;
x_end = x_start + filt_dur;
y_pos = 6;
x_ax_start = x_start - x_shift;
x_ax_end = x_ax_start + x_ax_length;
plot(linspace(x_start, x_end, sr),y + y_pos + y_shift, 'LineWidth', 4, 'Color', clrs(1,:))
plot([x_ax_start,x_ax_end], [y_pos y_pos], 'k', 'LineWidth', 2)
plot([x_ax_start, x_ax_start], [y_pos, y_pos+1+(2*y_shift)], 'k', 'LineWidth', 2)

x_start = x_end + 1;
x_end = x_start + filt_dur;
y_pos = 4;
plot(linspace(x_start, x_end, sr),y + y_pos + y_shift, 'LineWidth', 4, 'Color', clrs(2,:))
plot([x_ax_start,x_ax_end], [y_pos y_pos], 'k', 'LineWidth', 2)
plot([x_ax_start, x_ax_start], [y_pos, y_pos+1+(2*y_shift)], 'k', 'LineWidth', 2)

x_start = x_end + 1;
x_end = x_start + filt_dur;
x_start = x_end + 1;
x_end = x_start + filt_dur;
y_pos = 1;
plot(linspace(x_start, x_end, sr),y + y_pos + y_shift, 'LineWidth', 4, 'Color', clrs(3,:))
plot([x_ax_start,x_ax_end], [y_pos y_pos], 'k', 'LineWidth', 2)
plot([x_ax_start, x_ax_start], [y_pos, y_pos+1+(2*y_shift)], 'k', 'LineWidth', 2)
text(mean([x_ax_start, x_ax_end]), 3, '. . .', 'HorizontalAlignment','center', 'FontSize', 20, 'Rotation', 90,'FontWeight','bold')

x_start = x_end + gap;
x_end = x_start + signal_length - 1;

plot(linspace(x_start, x_end, length(orig_signals)), smooth(orig_signals(1,:), 5) + 2, 'LineWidth', 4, 'Color', 0.85*ones(1,3))    
plot(linspace(x_start, x_end, length(orig_signals)), signals(1,:) + 2, 'LineWidth', 4, 'Color', clrs(3,:))

plot(linspace(x_start, x_end, length(orig_signals)), smooth(orig_signals(2,:), 5) + 5, 'LineWidth', 4, 'Color', 0.85*ones(1,3))    
plot(linspace(x_start, x_end, length(orig_signals)), signals(2,:) + 5, 'LineWidth', 4, 'Color', clrs(2,:))

plot(linspace(x_start, x_end, length(orig_signals)), smooth(orig_signals(3,:), 5) + 7, 'LineWidth', 4, 'Color', 0.85*ones(1,3))    
plot(linspace(x_start, x_end, length(orig_signals)), signals(3,:) + 7, 'LineWidth', 4, 'Color', clrs(1,:))

text(mean([x_start, x_end]), 3, '. . .', 'HorizontalAlignment','center', 'FontSize', 20, 'Rotation', 90,'FontWeight','bold')

xlim([-x_shift, x_end+x_shift])
axis off
ax = gca;
ax.Position = [0 -0.1 1 1.2];

save_path = '~/Documents/om/projects/HearingInNoise/hearing_in_noise/figures/fig1';
save_name = fullfile(save_path, 'fixed_diagram.pdf');
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');

%% noise modeling

rng(10)

n_samples = 100;
y = 0*ones(1, n_samples);
noise = randn(size(y))/10;
signals(1,:) = smooth(y+noise, 5);

left_lim = -10;
right_lim = 20;
n_samples = 100;
offset = 7;
slope = 2;
max_scale = 1;
y = normpdf(linspace(left_lim,right_lim,n_samples), offset, slope);
y = (max_scale * y / max(y)) - 0.9;
noise = randn(size(y))/10;
signals(2,:) = smooth(y+noise, 5);

y = -0.25*ones(1, n_samples);
noise = randn(size(y))/10;
signals(3,:) = smooth(y+noise, 5);

f = 1/2;
dur = 1;
sr = 1000;
t = linspace(0, dur, dur*sr);
y = sin(2*pi*f*t);

figure('units','inches','innerposition',[0 0 7 3])
hold on

x_end = 0;
filt_dur = 24;
x_shift = 5;
x_ax_length = (filt_dur+1)*4 + 2*x_shift;
y_shift = 0.1;
gap = 50;
signal_length = 150;
foreground_length = signal_length / 3.25 * 0.5;
foreground_position = 1.25;

x_start = x_end + 1;
x_end = x_start + signal_length - 1;
plot(linspace(x_start, x_end, length(background)), background, 'Color', 0.6*ones(1,3))
foreground_start = x_start + signal_length / 3.25 * foreground_position;
foreground_end = foreground_start + foreground_length;
plot(linspace(foreground_start, foreground_end, length(foreground)), foreground, 'Color', 'k')

x_start = x_end + gap;
x_end = x_start + filt_dur;
y_pos = 6;
x_ax_start = x_start - x_shift;
x_ax_end = x_ax_start + x_ax_length;
plot(linspace(x_start, x_end, sr),y + y_pos + y_shift, 'LineWidth', 4, 'Color', clrs(1,:))
plot([x_ax_start,x_ax_end], [y_pos y_pos], 'k', 'LineWidth', 2)
plot([x_ax_start, x_ax_start], [y_pos, y_pos+1+(2*y_shift)], 'k', 'LineWidth', 2)

x_start = x_end + 1;
x_end = x_start + filt_dur;
y_pos = 4;
plot(linspace(x_start, x_end, sr),y + y_pos + y_shift, 'LineWidth', 4, 'Color', clrs(2,:))
plot([x_ax_start,x_ax_end], [y_pos y_pos], 'k', 'LineWidth', 2)
plot([x_ax_start, x_ax_start], [y_pos, y_pos+1+(2*y_shift)], 'k', 'LineWidth', 2)

x_start = x_end + 1;
x_end = x_start + filt_dur;
x_start = x_end + 1;
x_end = x_start + filt_dur;
y_pos = 1;
plot(linspace(x_start, x_end, sr),y + y_pos + y_shift, 'LineWidth', 4, 'Color', clrs(3,:))
plot([x_ax_start,x_ax_end], [y_pos y_pos], 'k', 'LineWidth', 2)
plot([x_ax_start, x_ax_start], [y_pos, y_pos+1+(2*y_shift)], 'k', 'LineWidth', 2)
text(mean([x_ax_start, x_ax_end]), 3, '. . .', 'HorizontalAlignment','center', 'FontSize', 20, 'Rotation', 90,'FontWeight','bold')

x_start = x_end + gap;
x_end = x_start + signal_length - 1;

plot(linspace(x_start, x_end, length(orig_signals)), signals(1,:) + 2, 'LineWidth', 4, 'Color', clrs(3,:))
plot(linspace(x_start, x_end, length(orig_signals)), signals(2,:) + 5, 'LineWidth', 4, 'Color', clrs(2,:))
plot(linspace(x_start, x_end, length(orig_signals)), signals(3,:) + 7, 'LineWidth', 4, 'Color', clrs(1,:))

text(mean([x_start, x_end]), 3, '. . .', 'HorizontalAlignment','center', 'FontSize', 20, 'Rotation', 90,'FontWeight','bold')
xlim([-x_shift, x_end+x_shift])
ylim([0, 8])
axis off
ax = gca;
ax.Position = [0 -0.1 1 1.2];

save_path = '~/Documents/om/projects/HearingInNoise/hearing_in_noise/figures/fig1';
save_name = fullfile(save_path, 'noisemodeling_diagram.pdf');
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');